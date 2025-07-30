library(tidyverse)
library(dtw)
library(glue)
library(furrr)

plan(multisession)
handlers(global = TRUE)

outdir <- "data/synchrony"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# Step 1: Create trial meta (one row per condition/word/order)
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
# Step 2: Parallel loop over all trials
all_sync <- future_pmap_dfr(trial_types, function(condition, word) {
  outdir <- "data/synchrony"
  files   <- trial_meta %>%
    filter(condition == !!condition, word == !!word) %>%
    pull(file)
  
  if (length(files) < 2) return(tibble())
  # Extract participant IDs from path
  participant_ids <- str_extract(files, "LexVar\\d+")
  
  # Read & tag exposures
  raw_list <- map(files, ~ read_csv(.x, show_col_types = FALSE) %>%
                    mutate(exposure = case_when(
                      frame %in%    0:1450   ~ "1",
                      frame %in% 3450:4800   ~ "2",
                      frame %in% 8425:10000  ~ "3",
                      frame %in% 11950:15000  ~ "4",
                      TRUE                   ~ NA_character_
                    )) %>%
                    filter(!is.na(exposure)))
  names(raw_list) <- participant_ids
  
  exposures <- c("1", "2", "3", "4")
  
  map_dfr(exposures, function(exp) {
    segs <- raw_list %>%
      keep(~ nrow(filter(.x, exposure == exp)) >= 2) %>%
      map(~ {
        df_exp <- filter(.x, exposure == exp) %>% arrange(frame)
        pts <- slice(df_exp, round(seq(1, n(), length.out = min(n(), 250))))$pupil_norm
        pts
      })
    ids <- names(segs)
    n   <- length(segs)
    
    if (n < 2) return(tibble())
    
    # All pairs (i, j) with i < j
    pairs <- combn(seq_len(n), 2, simplify = FALSE)
    dtw_dist <- function(x, y) dtw(x, y, distance.only = TRUE)$distance
    dists <- future_map_dbl(
      pairs,
      function(idx) dtw_dist(segs[[idx[1]]], segs[[idx[2]]]),
      .options = furrr_options(seed = TRUE)
    )
    dist_mat <- matrix(NA_real_, n, n, dimnames = list(ids, ids))
    for (k in seq_along(pairs)) {
      i <- pairs[[k]][1]; j <- pairs[[k]][2]
      dist_mat[i, j] <- dist_mat[j, i] <- dists[k]
    }
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
    { # Save per-trial result
      outfile <- glue("{outdir}/sync_{condition}_{word}.csv")
      write_csv(., outfile)
      .
    }
})

# Step 3: Save combined results
write_csv(all_sync, "data/pupil_synchrony_by_trial.csv")

# analysis ----------------------------------------------------------------
library(glmmTMB)
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
                            levels = c("unnaturalv", "stable", "naturalv"))) 

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

sync_model <- lmer(sync_z ~ condition + (1 | pid), data = combined_data)
summary(sync_model)

outliers <- performance::check_outliers(sync_model)
outlier_ids <- which(outliers >= 0.5)
combined_data_clean <- combined_data[-outlier_ids, ]

sync_model_clean <- lmer(sync_z ~ condition + (1 | pid), data = combined_data_clean)
summary(sync_model_clean)

ggplot(combined_data_clean, aes(x = condition, y = sync_z, fill = condition, color = condition)) +
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

ggplot(combined_data_clean, aes(x = sync_z, y = acc_z, color = condition, fill = condition)) +
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

m_null <- lmer(acc_z ~ (1 | pid), data = combined_data_clean)
m_sync <- lmer(acc_z ~ sync_z + (1 | pid), data = combined_data_clean)
anova(m_null, m_sync)
