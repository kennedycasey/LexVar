library(tidyverse)

# Set path to folder with .txt files
path <- "../data/raw/"  # <- replace with actual path
files <- list.files(path, pattern = "\\.txt$", full.names = TRUE)

gaze_data <- tibble()
# Read and combine all files
for (i in files) {
  curr_file <- read.table(i, header = TRUE, sep = "\t", fill = TRUE, quote = "") 
  
  print(unique(curr_file$center))
    
    curr_file <- curr_file %>% 
      transmute(id = RECORDING_SESSION_LABEL, 
              trial = TRIAL_INDEX, 
              x = CURRENT_FIX_X, 
              y = CURRENT_FIX_Y, 
              center = as.character(center)) %>%
    filter(trial %in% c(20, 23, 26, 29, 45, 48, 51, 54, 70, 73, 76, 79))
  
  gaze_data <- bind_rows(gaze_data, curr_file)
}

gaze_data_clean <- gaze_data %>%
  drop_na(x, y)

# Plot distribution of X values
gaze_data_clean <- gaze_data_clean %>%
  mutate(aoi = case_when(
  x >= 100 & x <= 680 ~ "left", 
  x >= 760 & x <= 1340 ~ "right", 
  x < 100 | x > 1340 ~ "off"
))

ggplot(gaze_data_clean, aes(x = x, fill = aoi)) +
  geom_histogram(binwidth = 20) +
  labs(title = "Distribution of Gaze X Positions",
       x = "X Position (pixels)", y = "Count") +
  theme_minimal()

# Plot distribution of Y values
ggplot(gaze_data_clean, aes(x = y)) +
  geom_histogram(binwidth = 20, fill = "darkorange", color = "black") +
  labs(title = "Distribution of Gaze Y Positions",
       x = "Y Position (pixels)", y = "Count") +
  theme_minimal()
