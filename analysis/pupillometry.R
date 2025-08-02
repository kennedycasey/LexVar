library(dtw)
library(glue)
library(furrr)
library(slider)

# full exposure sentence timings (switching to word-level for now)
exposure_frames <- list(
  "1" = 0:1450,
  "2" = 3450:4800,
  "3" = 8425:10000,
  "4" = 11950:15000
)

# set a reproducible seed for all future workers
options(future.rng.onMisuse = "error")
set.seed(20250801)

# parallelize
plan(multisession)
furrr_opts <- furrr_options(seed = TRUE)

if (reprocess.pupil) {
  
  # get window times for each word exposure
  exposure.label.times <- read_csv("data/metadata/exposure-label-times.csv") %>%
    pivot_longer(
      cols = starts_with("exp"),
      names_to = c("exposure", "bound"),
      names_pattern = "exp(\\d)_(onset|offset)",
      values_to = "frame"
    ) %>%
    pivot_wider(names_from = "bound", values_from = "frame") %>%
    mutate(exposure = as.character(exposure))
  
  outdir <- "data/synchrony"
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  # create metadata
  files <- list.files("data/processed/pupillometry", recursive = TRUE, full.names = TRUE)
  trial_meta <- tibble(
    file = files,
    participant = str_extract(file, "LexVar\\d+"),
    fname = basename(file),
    condition = str_extract(fname, "^[a-z]+"),
    word = str_extract(fname, "_[a-z]+_") %>% str_remove_all("_"),
    order = str_extract(fname, "(?<=_)[0-9]+(?=\\.csv$)")
  )
  
  trial_types <- trial_meta %>% distinct(condition, word)
  # loop over all trials
  # for each comibination of condition and word ...
  all_sync <- future_pmap_dfr(
    trial_types, 
    function(condition, word) {
      outdir <- "data/synchrony"
      
      # get relevant files for each participant
      files   <- trial_meta %>%
        filter(condition == !!condition, word == !!word) %>%
        pull(file)
      
      # require 2+ participants
      if (length(files) < 2) return(tibble())
      # extract pids
      participant_ids <- str_extract(files, "LexVar\\d+")
      
      # read in files and tag exposures
      raw_list <- map(files, ~ {
        df <- read_csv(.x, show_col_types = FALSE) %>%
          mutate(
            participant = str_extract(.x, "LexVar\\d+"),
            condition = !!condition,
            word = !!word
            )
        
      df$exposure <- NA_character_
      
      # for the relevant exposures ...
      for (exp in c("1", "4")) {
        timing <- exposure.label.times %>%
          filter(word == !!word, condition == !!condition, exposure == exp)
        
        if (nrow(timing) == 1) {
          frame_window <- timing$onset:timing$offset
          df$exposure[df$frame %in% frame_window] <- exp
        }
      }
      
      df %>% filter(!is.na(exposure))
    })
    
    names(raw_list) <- participant_ids
    
    exposures <- c("1", "4")
    
    map_dfr(exposures, function(exp) {
      timing <- exposure.label.times %>%
        filter(word == !!word, condition == !!condition, exposure == exp)
      
      if (nrow(timing) != 1) return(tibble())  # skip if not found
      
      frame_window <- timing$onset:timing$offset
      
      # extract usable segments
      segs <- raw_list %>%
        map(~ {
          df_exp <- filter(.x, exposure == exp) %>%
            arrange(frame)
          
          if (nrow(df_exp) < 2) return(NULL)
          
          df_exp <- df_exp %>%
            complete(
              frame = frame_window,
              fill = list(
                participant = unique(df_exp$participant),
                condition   = unique(df_exp$condition),
                word        = unique(df_exp$word)
              )
            )
          
          # remove segments with <50% data
          if (mean(is.na(df_exp$pupil_norm)) > 0.5) return(NULL)
          
          df_exp <- df_exp %>% na.omit()
          
          # duration normalization - 100 time points per word
          pts <- slice(df_exp, round(seq(1, n(), length.out = min(n(), 100))))$pupil_norm
          
          if (sum(!is.na(pts)) < 2) return(NULL)
          
          pts
        }) %>%
        compact()
      
      ids <- names(segs)
      n   <- length(segs)
      
      if (n < 2) return(tibble())
      
      # for all unique participant pairs ...
      # calculate dtw
      pairs <- combn(seq_len(n), 2, simplify = FALSE)
      dtw_dist <- function(x, y) {
        tryCatch(
          dtw(x, y, distance.only = TRUE)$distance,
          error = function(e) NA_real_
        )
      }
      
      dists <- future_map_dbl(
        pairs,
        function(idx) dtw_dist(segs[[idx[1]]], segs[[idx[2]]]),
        .options = furrr_opts
      )
      dist_mat <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
      for (k in seq_along(pairs)) {
        i <- pairs[[k]][1]; j <- pairs[[k]][2]
        dist_mat[i, j] <- dist_mat[j, i] <- dists[k]
      }
      # convert distance to synchrony
      max_d    <- max(dist_mat, na.rm = TRUE)
      sync_mat <- max_d - dist_mat
      sync_mean <- rowMeans(sync_mat, na.rm = TRUE)
      tibble(
        participant = ids,
        condition   = condition,
        word        = word,
        exposure    = exp,
        synchrony   = sync_mean
      )
    }) %>%
      { # save per-trial result
        outfile <- glue("{outdir}/sync_{condition}_{word}.csv")
        write_csv(., outfile)
      }
  }, 
  .options = furrr_opts)
  write_csv(all_sync, "data/pupil_synchrony_by_trial.csv")
}

# analysis ----------------------------------------------------------------
library(glmmTMB)
library(lme4) 
library(lmerTest)
# library(brms) # Uncomment if using Bayesian model
# library(ggeffects) # Uncomment if plotting predictions

# 1. Load Data
sync_data <- read_csv("data/pupil_synchrony_by_trial.csv") %>%
  group_by(participant, condition, word, exposure) %>%
  summarise(mean_synchrony = mean(synchrony, na.rm = TRUE), .groups = "drop")

test_data <- read_csv("data/test-props.csv")

# 2. Prepare for Join
# Rename columns to ensure consistency
# Adjust these if your column names differ
sync_data <- sync_data %>%
  mutate(pid = gsub("ex", "", participant)) %>%
  group_by(pid, word, condition) %>%
  summarize(mean_synchrony = mean(mean_synchrony)) %>%
  ungroup() %>%
  mutate(sync_z = datawizard::standardize(mean_synchrony))

test_data <- test_data %>%
  transmute(pid = pid,
            word = base_word, 
            condition = condition, 
            prop = prop) %>%
  group_by(pid, word, condition) %>%
  summarize(mean_accuracy = mean(prop)) %>%
  ungroup() %>%
  mutate(acc_z = datawizard::standardize(mean_accuracy))

# 3. Join Datasets
combined_data <- left_join(test_data, sync_data,
                           by = c("pid", "word", "condition")) %>%
  ungroup() %>%
  filter(!is.na(sync_z), !is.na(acc_z)) %>%
  mutate(condition = factor(condition, 
                            levels = c("stable", "unnaturalv", "naturalv"))) 

sync_model <- lm(sync_z ~ condition, data = combined_data)
summary(sync_model)

#outliers <- performance::check_outliers(sync_model)
#outlier_ids <- which(outliers)
#combined_data_clean <- combined_data[-outlier_ids, ]

#sync_model <- lmer(sync_z ~ condition + (1 | pid), data = combined_data)
#summary(sync_model)

ggplot(combined_data, aes(x = condition, y = sync_z, fill = condition, color = condition)) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  geom_boxplot(outlier.shape = NA, color = "black", fill = NA) +
  scale_color_manual(values = condition.colors) + 
  scale_fill_manual(values = condition.colors) + 
  scale_x_discrete(labels = c("Atypical variation", "Consistent", "Typical variation")) + 
  labs(
    title = "Pupil Synchrony Across Conditions",
    x = "Condition",
    y = "Mean Pupil Synchrony"
  ) +
  theme_minimal() + 
  theme(legend.position = "none")

ggplot(combined_data, aes(x = sync_z, y = acc_z, color = condition, fill = condition)) +
  geom_point(shape = 21, alpha = 0.1) + 
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~ condition) +
  scale_color_manual(values = condition.colors) + 
  scale_fill_manual(values = condition.colors) + 
  labs(
    x = "Mean Pupil Synchrony",
    y = "Mean Accuracy (Proportion Target Looks)",
    title = "Relationship Between Pupil Synchrony and Accuracy",
    subtitle = "Faceted by Condition"
  ) +
  theme_test() + 
  theme(legend.position = "none")

m_null <- lmer(acc_z ~ condition + (1 | pid), data = combined_data)
m_sync <- lmer(acc_z ~ sync_z + condition + (1 | pid), data = combined_data)
anova(m_null, m_sync)