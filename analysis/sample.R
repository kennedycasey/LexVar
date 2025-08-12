library(tidyverse)
library(readxl)
library(kableExtra)

log <- read_xlsx("../../2_data/LexVar-participant-log.xlsx") %>%
  mutate(dob = as.Date(dob, format="%d-%m-%y"), 
         dot = as.Date(dot, format="%d-%m-%y"), 
         age_d = as.numeric(dot - dob),
         age = age_d/30.4375)

n.tested <- nrow(log)
n.usable <- log %>%
  filter(usable_checked == "y") %>%
  nrow()

usable <- log %>%
  filter(usable_checked == "y")

usable %>%
  summarize(mean = round(mean(age), 1),
            median = round(median(age), 1),
            min = round(min(age), 1),  
            max = round(max(age), 1))

usable %>%
  mutate(total = n()) %>%
  group_by(sex, total) %>%
  summarize(n = n(), 
            prop = n/total) %>%
  distinct()

usable %>%
  mutate(ethnicity = ifelse(str_detect(ethnicity, "/"), "mixed", ethnicity), 
         total = n()) %>%
  group_by(ethnicity, total) %>%
  summarize(n = n(), 
            prop = n/total) %>%
  distinct() %>%
  mutate(across(where(is.numeric), ~ round(., 3) * 100))

log %>%
  filter(!is.na(usable_checked)) %>%
  mutate(age = floor(age), 
         usable_checked = factor(usable_checked,
                                 levels = c("n", "y", "not tested"), 
                                 labels = c("excluded", 
                                            "included", 
                                            "not tested"))) %>%
  group_by(age, usable_checked) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  ggplot(aes(x = age, y = n, color = usable_checked, fill = usable_checked)) +
  geom_bar(stat = "identity", alpha = 0.6, size = 1.5) + 
  geom_text(aes(label = n), color = "black", 
            position = position_stack(vjust = 0.9), size = 5) + 
  geom_text(aes(label = paste("Mean:", round(mean(usable$age), 1)), x = 34, y = 9), 
            color = "#5CAB87", size = 8) + 
  geom_text(aes(label = paste("Median:", round(median(usable$age), 1)), x = 34, y = 8.5), 
            color = "#5CAB87", size = 8) +
  scale_x_continuous(limits = c(17, 36), 
                     breaks = c(18, 21, 24, 27, 30, 33, 36)) +
  scale_color_manual(values = c("gray", "#5CAB87")) + 
  scale_fill_manual(values = c("gray", "#5CAB87")) + 
  labs(x = "Age (months)", y = "# of participants", 
       color = "", fill = "") +
  theme_classic(base_size = 20) + 
  theme(legend.position = "none")
ggsave("plots/age-dist-PREDICTED.png")

log %>%
  filter(!is.na(usable_checked)) %>%
  mutate(age = floor(age), 
         group = case_when(
           age <= 27 ~ "younger",
           age > 27 ~ "older"
         ),
         group = factor(group, levels = c("younger", "older")),
         usable_checked = factor(usable_checked,
                                 levels = c("n", "y", "not tested"), 
                                 labels = c("excluded", 
                                            "included", 
                                            "not tested"))) %>%
  group_by(group, usable_checked) %>%
  summarize(n = n()) %>%
  ungroup() %>%
  ggplot(aes(x = usable_checked, y = n, color = usable_checked, fill = usable_checked)) +
  facet_wrap(.~group) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.6, size = 1.5) + 
  geom_text(aes(label = n), color = "black", 
            position = position_stack(vjust = 0.9), size = 5) + 
  labs(x = "status", y = "# of participants", 
       color = "", fill = "") +
  theme_test(base_size = 15) + 
  theme(legend.position = "none", text = element_text(face = "bold"))
ggsave("plots/age-breakdown.png")