library(tidyverse)
library(readxl)

orders <- read_csv("data/metadata/orders.csv")

path <- "data/raw/"

log <- read_xlsx("data/metadata/ANON-ptcp-log.xlsx")

label.times <- read_csv("data/metadata/label-times.csv")

create_frames <- function(start_time, end_time, frame_duration = 1) {
  seq(start_time, end_time, by = frame_duration)
}

test_props <- tibble()
frames <- tibble()

included.files <- list.files(path, ".txt")
usable <- tibble()

window_max_length = 2000
window_start_delay = 300

for (i in included.files) {
  
  curr_file <- read.table(paste0(path, i), header = TRUE, sep = "\t", fill = TRUE, quote = "")
  
  id <- str_extract(unique(curr_file$RECORDING_SESSION_LABEL), "[0-9][0-9][0-9]")
  
  order <- filter(log, str_detect(ptcp_id, id))$order
  
  age <- filter(log, str_detect(ptcp_id, id))$age
  
  age_m <- filter(log, str_detect(ptcp_id, id))$age_m
  
  age_d <- filter(log, str_detect(ptcp_id, id))$age_d
  
  status <- filter(log, str_detect(ptcp_id, id))$usable
  
  if (status != "n") {
    
    exists <- file.exists(paste0("data/processed/", id, "-test.csv"))
    
    if (reprocess | !exists) {
    
    curr_file <- curr_file %>%
      transmute(pid = RECORDING_SESSION_LABEL,
                id = id,
                order = order,
                age = age,
                age_m = age_m,
                age_d = age_d,
                trial_index = TRIAL_INDEX,
                trial_label = TRIAL_LABEL,
                trial_start = TRIAL_START_TIME,
                center = center,
                aoi = case_when(
                  CURRENT_FIX_X <= 680 & CURRENT_FIX_X >= 100 ~ "left",
                  CURRENT_FIX_X >= 760 & CURRENT_FIX_X <= 1340 ~ "right",
                  TRUE ~ NA), 
                x = CURRENT_FIX_X,
                fix_start = CURRENT_FIX_START,
                fix_end = CURRENT_FIX_END,
                duration = CURRENT_FIX_DURATION,
                pupil = CURRENT_FIX_PUPIL,
                video = str_extract(videofile, "([^/]+)$"),
                img = str_extract(imgfile, "([^/]+)$"),
                wav = str_extract(wavfilename, "([^/]+)$")) %>%
      filter(wav != "ITI") %>%
      distinct()
    
    curr_order <- orders %>%
      filter(Order == order) %>%
      select(-Order) %>%
      pivot_longer(everything()) %>%
      filter(!str_detect(name, "object"))
    
    conditions <- curr_order %>%
      slice_head(n = 6) %>%
      transmute(base_word = substr(name, 1, 3), 
                condition = value)
    
    target_sides <- curr_order %>%
      filter(str_detect(name, "targetside")) %>%
      transmute(trial = as.numeric(str_remove(name, "targetside_")), 
                target_side = value)
    
    test <- curr_file %>%
      filter(trial_index %in% c(20, 23, 26, 29, 45, 48, 51, 54, 70, 73, 76, 79)) %>%
      group_by(trial_index) %>%
      mutate(trial = cur_group_id()) %>%
      left_join(target_sides, by = "trial") %>%
      mutate(base_word = substr(wav, 1, 3)) %>%
      left_join(conditions, by = "base_word") %>%
      left_join(label.times, by = "base_word") %>%
      filter(aoi %in% c("left", "right")) %>%
      rowwise() %>%
      mutate(frame = list(create_frames(fix_start, fix_end))) %>%
      unnest(cols = frame) %>%
      mutate(frame = round(frame - first_label, 0), 
             aoi = toupper(substr(aoi, 1, 1)), 
             fix = ifelse(aoi == target_side, "target", "distracter")) %>%
      ungroup() %>%
      mutate(
        first_start = window_start_delay,
        first_end = window_max_length,
        second_start = round(second_label - first_label, 0) + window_start_delay,
        second_end = round(second_label - first_label, 0) + window_max_length,
        window = case_when(
          frame >= first_start & frame <= first_end ~ "first",
          frame >= second_start & frame <= second_end ~ "second",
          TRUE ~ NA_character_
        ),
      postnaming = ifelse(frame >= window_start_delay & frame <= second_label-first_label + window_max_length, "y", "n"))
    
    trial.conditions <- test %>%
      select(trial, condition) %>%
      distinct()
    
    prenaming_window_fix <- test %>%
      filter(frame <= 0) %>%
      pull(trial) %>%
      unique()
    
    first_window_fix <- test %>%
      filter(window == "first") %>%
      group_by(trial) %>%
      summarize(total = n()) %>%
      filter(total >= (window_max_length-window_start_delay)/3) %>%
      pull(trial)
    
    second_window_fix <- test %>%
      filter(window == "second") %>%
      group_by(trial) %>%      
      summarize(total = n()) %>%
      filter(total >= (window_max_length-window_start_delay)/3) %>%
      pull(trial)

    usable.trials.nos <- c(intersect(prenaming_window_fix, first_window_fix),
                           second_window_fix) %>%
      unique()
    
    usable.trials <- tibble(trial = usable.trials.nos) %>%
      left_join(trial.conditions, by = "trial") %>%
      mutate(id = id) %>%
      group_by(id, condition) %>%
      mutate(total_trials_per_block = n()) %>%
      filter(total_trials_per_block >= 2) %>%
      group_by(id) %>%
      mutate(total_usable_blocks = length(unique(condition))) %>%
      filter(total_usable_blocks >= 2) %>%
      select(id, condition, trial)
    
    if (nrow(usable.trials) > 0) {
      train <- curr_file %>%
        filter(trial_index %in% c(2, 5, 8, 11, 33, 36, 39, 42, 58, 61, 64, 67)) %>%
        group_by(trial_index) %>%
        mutate(trial = cur_group_id()) %>%
        mutate(base_word = substr(wav, 1, 3)) %>%
        left_join(conditions, by = "base_word") %>%
        rowwise() %>%
        mutate(frame = list(create_frames(fix_start, fix_end))) %>%
        unnest(cols = frame) %>%
        ungroup()

      write_csv(usable.trials, paste0("data/processed/usable/", id, ".csv"))
      
      write_csv(test, paste0("data/processed/", id, "-test.csv"))
      write_csv(train, paste0("data/processed/", id, "-train.csv"))
      
      # cluster <- test %>%
      #   arrange(frame) %>%
      #   group_by(pid, id, order, age_m, trial_index, trial_label, trial_start, 
      #            center, trial, target_side, base_word, condition, first_label, 
      #            second_label, window, postnaming) %>%
      #   group_modify(~ {
      #     full_frame_range <- data.frame(frame = -3152:7471)
      #     full_data <- full_frame_range %>%
      #       left_join(.x, by = "frame")
      #     
      #     full_data_filled <- full_data %>%
      #       tidyr::fill(names(.), .direction = "down") %>%
      #       distinct()
      #     
      #     full_data_filled <- full_data_filled %>%
      #       mutate(
      #         aoi = ifelse(is.na(aoi), NA, aoi),
      #         trackloss = ifelse(is.na(aoi), TRUE, FALSE)
      #       )
      #     return(full_data_filled)
      #   }) %>%
      #   ungroup()
      # 
      # write.csv(cluster, paste0("data/processed/", id, "-cluster.csv"), row.names = FALSE)
    }
    }
  }
}

counts <- c()
for (i in list.files("data/processed/usable/", ".csv")) {
  curr <- read_csv(paste0("data/processed/usable/", i)) %>%
    nrow()

  counts <- c(counts, curr)
}

# pupillometry ------------------------------------------------------------
library(stinepack)
library(zoo)

processed.path <- "data/processed/"
for (i in list.files(processed.path, "train.csv")) {
  curr.id <- str_extract(i, "[0-9][0-9][0-9]")

  # create dirs for participant data
  if (!dir.exists(paste0(processed.path, "pupillometry/LexVar", curr.id))) {
    dir.create(paste0(processed.path, "pupillometry/LexVar", curr.id))
  }

  curr.file <- read_csv(paste0(processed.path, i)) %>%
    group_by(base_word, trial_index) %>%
    mutate(order = cur_group_id()) %>%
    ungroup() %>%
    group_by(base_word) %>%
    mutate(word.order = ifelse(order == min(order), 1, 2)) %>%
    ungroup()

  for (j in unique(curr.file$order)) {
    curr.trial <- curr.file %>%
      filter(order == j)
    
    # convert raw pupil size to % change
    curr.trial <- curr.trial %>%
        mutate(pupil_norm = 100 * (pupil - mean(pupil, na.rm = TRUE)) / mean(pupil, na.rm = TRUE))

    # set window size for artifact detection
    window_size <- 50
    
    pupil_range <- rollapply(curr.trial$pupil_norm, width = window_size, FUN = function(x) max(x) - min(x),
                             fill = NA, align = "center")
    
    # drop blink artifacts (>15% change in 0.05-s window)
    artifact_mask <- pupil_range > 15
    
    curr.trial <- curr.trial %>%
      mutate(pupil_norm = pupil_norm,
             artifact = artifact_mask) %>%
      filter(!artifact)
    
    curr.word <- unique(curr.trial$base_word)
    curr.word.order <- unique(curr.trial$word.order)
    curr.condition <- unique(curr.trial$condition)
    
    # interpolate gaps <100ms (stineman)
    curr.trial <- curr.trial %>%
      mutate(
        pupil_norm_interp = {
          if (sum(!is.na(pupil_norm)) < 2) {
            rep(NA_real_, length(pupil_norm))
          } else {
            stinterp(
              x = which(!is.na(pupil_norm)),
              y = pupil_norm[!is.na(pupil_norm)],
              xout = seq_along(pupil_norm)
            )$y
          }
        }
      )

    curr.trial %>%
      transmute(frame = row_number(),
                pupil_norm = pupil_norm) %>%
      write_csv(paste0("data/processed/pupillometry/LexVar", curr.id, "/",
                       curr.condition, "_", curr.word, "_", curr.word.order, ".csv"))
  }
}