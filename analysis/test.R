library(tidyverse)
library(readxl)
library(tidybayes)
library(brms) # run bayes lm 
library(marginaleffects) # get posteriors 
library(ggeffects) # graph
library(easystats) # easystats packages # bayesttestR 
library(bayesplot) # graph trace plots
library(ggdist)

reprocess = FALSE
rewrite = FALSE

if (reprocess) {
  source("process.R")
}

path <- "processed/"

log <- read_xlsx("../../2_data/LexVar-participant-log.xlsx")

condition.colors <- c("stable" = "#ef5039", 
                      "unnaturalv" = "#04a48c", 
                      "naturalv" = "#3f5eaa")


if (rewrite) {
  test_props <- tibble()
  frames <- tibble()
  
  included.files <- list.files(path, "-test.csv")
  
  for (i in included.files) {
    
    curr_file <- read_csv(paste0(path, i))
    
    id <- str_extract(unique(curr_file$pid), "[0-9][0-9][0-9]")
    
    status <- filter(log, str_detect(ptcp_id, id))$usable_checked
    
    #status <- filter(log, str_detect(ptcp_id, id))$usable_strict
    
    if (status != "n" & !is.na(status)) {
      # read in usable trials from pre-processing stage
      usable_trials <- read_csv(paste0("usable/", id, ".csv")) %>%
        pull(trial)
      
      if (length(usable_trials) > 0) {
        
      curr_test_props <- curr_file %>%
        select(pid, age_m, order, trial, x, pupil, base_word:window) %>%
        filter(trial %in% usable_trials & !is.na(window)) %>%
        group_by(pid, trial, condition, fix, base_word) %>%
        summarize(total = n()) %>%
        pivot_wider(names_from = "fix", values_from = "total", 
                    id_cols = c(pid, trial, condition, base_word))
      
      if ("target" %in% colnames(curr_test_props)) {
        curr_test_props <- curr_test_props %>%
          mutate(target = ifelse(is.na(target), 0, target))
      } else {
        curr_test_props <- curr_test_props %>%
          mutate(target = 0)
      }
      
      if ("distracter" %in% colnames(curr_test_props)) {
        curr_test_props <- curr_test_props %>%
          mutate(distracter = ifelse(is.na(distracter), 0, distracter))
      } else {
        curr_test_props <- curr_test_props %>%
          mutate(distracter = 0)
      }
      
      curr_test_props <- curr_test_props %>%
        mutate(total = target + distracter, 
               prop = target/total)
      
      test_props <- bind_rows(test_props, curr_test_props)
      
      curr_frames <- curr_file %>%
        select(pid, age_m, order, trial, x, pupil, base_word:window) %>%
        filter(trial %in% usable_trials) %>%
        mutate(fix_num = case_when(
          fix == "distracter" ~ 0, 
          fix == "target" ~ 1, 
          TRUE ~ NA
        ))
      
      frames <- bind_rows(frames, curr_frames)
      
      }
    }
  }
  
    data <- frames %>%
      mutate(time = frame, 
             condition = factor(condition, levels = c("unnaturalv", "stable", "naturalv"))) %>%
      select(pid, age_m, order, trial, condition, base_word, pupil, time, fix_num, window, first_label, second_label) %>%
      group_by(pid) %>%
      mutate(block_order = list(unique(condition))) %>%
      ungroup() %>%
      mutate(condition_order = map2_int(block_order, condition, ~which(.x == .y)))
    write_csv(data, "prepped-data.csv")
    
    # export eyetrackingR version of data
    # eyetrackingr_data <- data %>%
    #   mutate(target = ifelse(fix_num == 1, 1, 0), 
    #          distracter = ifelse(fix_num == 0, 1, 0)) %>%
    #   select(-fix_num) %>%
    #   write_csv("eyetrackingr-data.csv")
    
    write_csv(test_props, "test-props.csv")
}

data <- read_csv("prepped-data.csv")
test_props <- read_csv("test-props.csv")


# corr between repeated measures for power --------------------------------

library(tidyverse)

# Load data
test_props <- read_csv("test-props.csv")

rep.measures.corrs <- test_props %>%
  group_by(pid, condition, base_word) %>%
  select(pid, condition, base_word, prop) %>%
  mutate(n = paste0("trial", row_number())) %>%
  pivot_wider(names_from = n, values_from = prop) %>%
  ungroup()

library(rmcorr)

# Repeated measures correlation between trial1 and trial2 across base_words within each participant
rmcorr_result <- rmcorr(participant = pid, measure1 = trial1, measure2 = trial2, dataset = rep.measures.corrs)
rmcorr_result

# analysis ----------------------------------------------------------------

bypid.bywindow.condition.summary <- data %>%
  filter(!is.na(window)) %>%
  mutate(condition = factor(condition, levels = c("stable", "naturalv", "unnaturalv")), 
         window = factor(window, levels = c("first", "second"), 
                         labels = c("post-naming window 1", 
                                    "post-naming window 2"))) %>%
  group_by(pid, condition, window, age_m) %>%
  summarize(accuracy = mean(fix_num)) %>%
  ungroup()

data %>%
  filter(!is.na(window)) %>%
  mutate(condition = factor(condition, levels = c("stable", "naturalv", "unnaturalv")), 
         window = factor(window, levels = c("first", "second"), 
                         labels = c("post-naming window 1", 
                                    "post-naming window 2"))) %>%
  group_by(pid, condition, window) %>%
  summarize(accuracy = mean(fix_num)) %>%
  ungroup() %>%
  group_by(condition, window) %>%
  summarize(mean = mean(accuracy), 
            sd = sd(accuracy), 
            ci = 1.96*sd/sqrt(n())) %>%
  ggplot(aes(x = condition, y = mean, color = condition, fill = condition)) + 
  facet_grid(.~window) +
  geom_hline(yintercept = 0.5, color = "gray", linetype = "dashed", size = 1.25) + 
  geom_bar(stat = "identity", alpha = 0.6, size = 2) + 
  geom_jitter(data = bypid.bywindow.condition.summary, 
              mapping = aes(x = condition, y = accuracy, color = condition), 
              width = 0.1, size = 3) + 
  geom_errorbar(aes(ymin = mean-ci, ymax = mean+ci), color = "black", 
                width = 0.2, size = 1.25) + 
  scale_color_manual(values = condition.colors) + 
  scale_fill_manual(values = condition.colors) + 
  scale_x_discrete(labels = c("consistent", "typical\nvariation", "atypical\nvariation")) +
  labs(x = "Condition", y = "Accuracy\n(Proportion looking to target object)") + 
  theme_test(base_size = 20) + 
  theme(legend.position = "none")
ggsave("plots/bywindow-condition.png")

data %>%
  filter(!is.na(window)) %>%
  mutate(condition = factor(condition, levels = c("stable", "naturalv", "unnaturalv")), 
         window = factor(window, levels = c("first", "second"), 
                         labels = c("post-naming window 1", 
                                    "post-naming window 2"))) %>%
  group_by(pid, condition, window) %>%
  summarize(accuracy = mean(fix_num)) %>%
  ungroup() %>%
  group_by(condition, window) %>%
  summarize(mean = mean(accuracy), 
            sd = sd(accuracy), 
            ci = 1.96*sd/sqrt(n())) %>%
  filter(window == "post-naming window 1") %>%
  ggplot(aes(x = condition, y = mean, color = condition, fill = condition)) + 
  facet_grid(.~window) +
  geom_hline(yintercept = 0.5, color = "gray", linetype = "dashed", size = 1.25) + 
  geom_bar(stat = "identity", alpha = 0.6, size = 2) + 
  geom_jitter(data = filter(bypid.bywindow.condition.summary, 
                            window == "post-naming window 1"), 
              mapping = aes(x = condition, y = accuracy, color = condition), 
              width = 0.1, size = 3) + 
  geom_errorbar(aes(ymin = mean-ci, ymax = mean+ci), color = "black", 
                width = 0.2, size = 1.25) + 
  scale_color_manual(values = condition.colors) + 
  scale_fill_manual(values = condition.colors) + 
  scale_x_discrete(labels = c("consistent", "typical\nvariation", "atypical\nvariation")) +
  labs(x = "Condition", y = "Accuracy\n(Proportion looking to target object)") + 
  theme_test(base_size = 20) + 
  theme(legend.position = "none", 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank())
ggsave("plots/bywindow-condition1.png", width = 4, height = 8)

data %>%
  filter(!is.na(window)) %>%
  mutate(condition = factor(condition, levels = c("stable", "naturalv", "unnaturalv")), 
         window = factor(window, levels = c("first", "second"), 
                         labels = c("post-naming window 1", 
                                    "post-naming window 2"))) %>%
  group_by(pid, condition, window) %>%
  summarize(accuracy = mean(fix_num)) %>%
  ungroup() %>%
  group_by(condition, window) %>%
  summarize(mean = mean(accuracy), 
            sd = sd(accuracy), 
            ci = 1.96*sd/sqrt(n())) %>%
  filter(window == "post-naming window 2") %>%
  ggplot(aes(x = condition, y = mean, color = condition, fill = condition)) + 
  facet_grid(.~window) +
  geom_hline(yintercept = 0.5, color = "gray", linetype = "dashed", size = 1.25) + 
  geom_bar(stat = "identity", alpha = 0.6, size = 2) + 
  geom_jitter(data = filter(bypid.bywindow.condition.summary, 
                            window == "post-naming window 2"), 
              mapping = aes(x = condition, y = accuracy, color = condition), 
              width = 0.1, size = 3) + 
  geom_errorbar(aes(ymin = mean-ci, ymax = mean+ci), color = "black", 
                width = 0.2, size = 1.25) + 
  scale_color_manual(values = condition.colors) + 
  scale_fill_manual(values = condition.colors) + 
  scale_x_discrete(labels = c("consistent", "typical\nvariation", "atypical\nvariation")) +
  labs(x = "Condition", y = "Accuracy\n(Proportion looking to target object)") + 
  theme_test(base_size = 20) + 
  theme(legend.position = "none", 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        axis.title.y = element_blank())
ggsave("plots/bywindow-condition2.png", width = 4, height = 8)

# ggplot(test_props, aes(x = condition, y = prop, color = condition, fill = condition)) +
#   facet_wrap(. ~ pid) +
#   stat_summary(geom = "bar", stat = "mean") +
#   stat_summary(geom = "errorbar", stat = "mean_se", width = 0.2, color = "black") +
#   scale_color_manual(values = condition.colors) +
#   scale_fill_manual(values = condition.colors) +
#   geom_hline(yintercept = 0.5) +
#   labs(x = "Condition", y = "Accuracy") +
#   theme_test(base_size = 10) +
#   theme(legend.position = "none")
# ggsave("plots/byparticipant-condition.png")
# 
# ggplot(test_props, aes(x = condition, y = prop, color = condition, fill = condition)) +
#   facet_wrap(. ~ base_word) +
#   stat_summary(geom = "bar", stat = "mean") +
#   stat_summary(geom = "errorbar", stat = "mean_se", width = 0.2, color = "black") +
#   scale_color_manual(values = condition.colors) +
#   scale_fill_manual(values = condition.colors) +
#   geom_hline(yintercept = 0.5)
# ggsave("plots/byword-condition.png")

bypid.condition.summary <- data %>%
  filter(!is.na(window)) %>%
  filter(window == "second") %>%
  mutate(condition = factor(condition, levels = c("stable", "naturalv", "unnaturalv"))) %>%
  group_by(pid, condition, age_m) %>%
  summarize(accuracy = mean(fix_num)) %>%
  ungroup()

condition.summary <- bypid.condition.summary %>%
  group_by(condition) %>%
  summarize(mean = mean(accuracy), 
            sd = sd(accuracy), 
            ci = 1.96*sd/sqrt(n()))

condition.summary.age <- bypid.condition.summary %>%
  # mutate(age_bin = case_when(
  #   age_m < 24 ~ "18-24", 
  #   age_m >= 24 & age_m < 30 ~ "24-30", 
  #   age_m >= 30 ~ "30-36"
  # )) %>%
  mutate(age_bin = ifelse(age_m <= median(age_m), "younger half", "older half"), 
         age_bin = factor(age_bin, levels = c("younger half", "older half")), 
         condition = factor(condition, levels = c("stable", "naturalv", "unnaturalv"))) %>%
  group_by(condition, age_bin) %>%
  summarize(mean = mean(accuracy), 
            sd = sd(accuracy), 
            ci = 1.96*sd/sqrt(n()))

ggplot(bypid.condition.summary, aes(x = condition, y = accuracy, 
                                    color = condition, fill = condition)) + 
  geom_hline(yintercept = 0.5, color = "gray", linetype = "dashed", size = 1.25) +
  ggrain::geom_rain(boxplot.args = list(fill = NA, color = "black", outlier.shape = NA),
                    violin.args = list(alpha = 0.5, color = "black"),
                    point.args = list(color = "black", shape = 21)) +
  # geom_segment(aes(x = as.numeric(condition[1]) - 0.2, xend = as.numeric(condition[1]) + 0.5,
  #                  y = filter(condition.summary, condition == "stable")$mean), size = 2, 
  #              color = condition.colors[1]) +
  # geom_segment(aes(x = as.numeric(condition[2]) - 0.2, xend = as.numeric(condition[2]) + 0.5,
  #                  y = filter(condition.summary, condition == "naturalv")$mean), size = 2,
  #              color = condition.colors[3]) +
  # geom_segment(aes(x = as.numeric(condition[3]) - 0.2, xend = as.numeric(condition[3]) + 0.5,
  #                  y = filter(condition.summary, condition == "unnaturalv")$mean), size = 2,
  #              color = condition.colors[2]) +
  labs(x = "Condition", y = "Accuracy\n(Proportion looking to target object)") + 
  scale_color_manual(values = condition.colors) + 
  scale_fill_manual(values = condition.colors) + 
  geom_text(aes(x = "stable", y = 1.1, label = "***")) +
  scale_x_discrete(labels = c("consistent", "typical\nvariation", "atypical\nvariation")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
  theme_classic(base_size = 20) + 
  theme(legend.position = "none", 
        title = element_text(face = "bold"))

ggsave("pres/descriptive-condition-ALL.png", width = 8, height = 6)

ggplot(condition.summary, aes(x = condition, y = mean, color = condition, fill = condition)) + 
  geom_hline(yintercept = 0.5, color = "gray", linetype = "dashed", size = 1.25) + 
  #geom_bar(stat = "identity", alpha = 0.6, size = 2) + 
  geom_pointrange(aes(ymin = mean-ci, ymax = mean+ci), width = 0.2) + 
  geom_jitter(data = bypid.condition.summary, 
             mapping = aes(x = condition, y = accuracy, color = condition), 
             width = 0.1, size = 1, alpha = 0.5) + 
  scale_color_manual(values = condition.colors) + 
  scale_fill_manual(values = condition.colors) + 
  scale_x_discrete(labels = c("consistent", "typical\nvariation", "atypical\nvariation")) + 
  labs(x = "Condition", y = "Accuracy\n(Proportion looking to target object)") + 
  theme_test(base_size = 20) + 
  theme(legend.position = "none")
ggsave("plots/overall-condition.png")

ggplot(condition.summary.age, aes(x = condition, y = mean, color = condition, fill = condition)) + 
  facet_wrap(. ~ age_bin) + 
  geom_hline(yintercept = 0.5, color = "gray", linetype = "dashed", size = 1.25) + 
  geom_pointrange(aes(ymin = mean-ci, ymax = mean+ci)) + 
  geom_jitter(data = bypid.condition.summary %>%
                # mutate(age_bin = case_when(
                #   age_m < 24 ~ "18-24", 
                #   age_m >= 24 & age_m < 30 ~ "24-30", 
                #   age_m >= 30 ~ "30-36"
                # )), 
                mutate(age_bin = ifelse(age_m <= median(age_m), "younger half", "older half"), 
                       age_bin = factor(age_bin, levels = c("younger half", "older half"))), 
              mapping = aes(x = condition, y = accuracy), 
              width = 0.1, size = 1, alpha = 0.5) + 
  scale_color_manual(values = condition.colors) + 
  scale_fill_manual(values = condition.colors) + 
  scale_x_discrete(labels = c("consistent", "typical\nvariation", "atypical\nvariation")) + 
  labs(x = "Condition", y = "Accuracy\n(Proportion looking to target object)") + 
  theme_test(base_size = 20) + 
  theme(legend.position = "none")
ggsave("plots/overall-condition-byage.png")

model.data <- data %>%
  #mutate(time = frame) %>%
  filter(window %in% c("first", "second")) %>%
  mutate(age = datawizard::center(age_m), 
         pid = as.factor(pid),
         condition = factor(condition, levels = c("unnaturalv", "stable", "naturalv")), 
         condition_unnat_vs_nat = ifelse(condition == "unnaturalv", -0.5,
                                       ifelse(condition == "naturalv", 0.5, 0)), 
         condition_unnat_vs_stab = ifelse(condition == "unnaturalv", -0.5,
                                          ifelse(condition == "stable", 0.5, 0)))

model.data.prop <- model.data %>%
  group_by(pid, age, condition, trial, base_word, order) %>%
  summarize(prop = mean(fix_num), 
            frames = n()) %>%
  ungroup() %>%
  mutate(frames = datawizard::center(frames)) %>%
  group_by(pid) %>%
  mutate(block_order = list(unique(condition))) %>%
  ungroup() %>%
  mutate(condition_order = map2_int(block_order, condition, ~which(.x == .y))) 


library(lme4)
library(lmerTest)

m <- glmer(fix_num ~ condition_unnat_vs_stab + condition_unnat_vs_nat + age + (1|pid),
         family = "binomial",
         data = model.data)

# maximal model
m <- lmer(prop ~ condition + age + (1+condition|pid) + (1+condition+age|base_word),
           #family = binomial(link = "logit"),
           control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5)),
           data = model.data.prop)

# simplified model that does converge and doesn't have singular fit
m <- lmer(prop ~ condition * age + (1 + condition|pid),
          control = control,
          data = model.data.prop)

summary(m)

performance::check_model(m)

model2 <- brm(
  bf(
    prop ~ condition + age + (1 + condition + age | order) + (1 + condition | pid),  # The mean of the 0-1 values, or mu
    phi ~ condition + age + (1 + condition + age | order) + (1 + condition | pid),  # The precision of the 0-1 values, or phi
    zoi ~ condition + age + (1 + condition + age | order) + (1 + condition | pid),  # The zero-or-one-inflated part, or alpha
    coi ~ condition + age + (1 + condition + age | order) + (1 + condition | pid)  # The one-inflated part, conditional on the 0s, or gamma
  ),
  save_pars=save_pars(all=TRUE),
  family = zero_one_inflated_beta(),
  data = model.data.prop,
  prior = prior(normal(0.5, 0.25), class = b),
  chains = 4,
  core = 4,
  warmup = 3000,
  iter = 10000
)

# estimating corr between repeated measures for BL power analysis

# Assuming your model is called 'beta_model'
post_samples <- posterior_samples(beta_model, pars = "^sd_")

# Compute ICC for each posterior sample
icc_samples <- post_samples$sd_pid__Intercept^2 / 
  (post_samples$sd_pid__Intercept^2 + (pi^2 / 3))  # Approximation for residual variance

# Summarize ICC estimates
icc_mean <- mean(icc_samples)
icc_ci <- quantile(icc_samples, probs = c(0.025, 0.975))

# Print results
print(paste0("Estimated ICC: ", round(icc_mean, 3)))
print(paste0("95% CI: [", round(icc_ci[1], 3), ", ", round(icc_ci[2], 3), "]"))

# model2 <- brm(
#   bf(
#     prop ~ condition + age + (1 + condition + age | base_word) + (1 + condition | pid),  # The mean of the 0-1 values, or mu
#     phi ~ condition + age,  # The precision of the 0-1 values, or phi
#     zoi ~ condition + age,  # The zero-or-one-inflated part, or alpha
#     coi ~ condition + age  # The one-inflated part, conditional on the 0s, or gamma
#   ),
#   save_pars=save_pars(all=TRUE),
#   family = zero_one_inflated_beta(),
#   data = model.data.prop,
#   prior = prior(normal(0.5, 0.25), class = b),
#   chains = 4,
#   core = 4,
#   warmup = 3000,
#   iter = 10000
# )
# 
# model3 <-  brm(
#   bf(
#     prop ~ condition + age + (1 + condition + age | base_word) + (1 + condition | pid),  # The mean of the 0-1 values, or mu
#     zoi ~ condition + age,  # The zero-or-one-inflated part, or alpha
#     coi ~ condition + age  # The one-inflated part, conditional on the 0s, or gamma
#   ),
#   save_pars=save_pars(all=TRUE),
#   family = zero_one_inflated_beta(),
#   data = model.data.prop,
#   prior = prior(normal(0.5, 0.25), class = b),
#   chains = 4,
#   core = 4,
#   warmup = 3000,
#   iter = 10000
# )

# 
# bayesfactor_models(model1, model2)
# bayesfactor_models(model1, model3)

beta_model <- model2
effects <- conditional_effects(beta_model)
plot(effects)

describe_posterior(
  beta_model,
  effects = "fixed",
  component = "all",
  centrality = "all",
  standardize = "refit"
)

plot(effects)$condition

diagnostic_posterior(beta_model)

bayesplot::color_scheme_set("mix-blue-red")
bayesplot::mcmc_trace(beta_model, pars = c("b_conditionnaturalv"), 
                      facet_args = list(ncol = 1, strip.position = "left"))

bayesplot::mcmc_trace(beta_model, pars = c("b_conditionstable"), 
                      facet_args = list(ncol = 1, strip.position = "left"))

pp_check(beta_model, ndraws=1000)
pp_check(beta_model, type = "stat", stat = "mean")

describe_posterior(
  beta_model,
  effects = "fixed",
  component = "all",
  centrality = "all"
)

# cohen's d
describe_posterior(
  beta_model,
  effects = "fixed",
  component = "all",
  standardize = "refit"
)

# power analysis for next study
library(brms)
library(bayestestR)

# 1. Extract draws
posterior_draws <- as.data.frame(as_draws_df(beta_model))

# 2. Simulate new data for a new sample size (e.g., 60 participants)
simulate_data <- function(n_participants, effect = 0.3) {
  pid <- factor(rep(1:n_participants, each = 6))            # length: n_participants * 6
  condition <- rep(c("unnaturalv", "stable", "naturalv"), times = n_participants)  # length: n_participants * 3
  condition <- factor(condition, levels = c("unnaturalv", "stable", "naturalv"))
  
  # Simulate prop outcome based on effect size (on logit scale)
  eta <- ifelse(condition == "stable", effect,
                ifelse(condition == "naturalv", effect, 0))
  prob <- plogis(eta)
  
  prop <- rbeta(n_participants * 3, shape1 = prob * 20, shape2 = (1 - prob) * 20)  # simulate proportions
  
  tibble(pid, condition, prop)
}

# 3. Run multiple simulations
library(tibble)

simulate_data <- function(n_participants, effect = 0.3) {
  # Repeat each participant 3 times (for 3 conditions)
  pid <- rep(1:n_participants, times = 3)
  condition <- rep(c("unnaturalv", "stable", "naturalv"), times = n_participants)
  condition <- factor(condition, levels = c("unnaturalv", "stable", "naturalv"))
  
  # Logit-based effect (assuming 0 for unnaturalv)
  eta <- ifelse(condition == "stable", effect,
                ifelse(condition == "naturalv", effect, 0))
  prob <- plogis(eta)
  
  # Simulate beta-distributed proportions (mu = prob, precision ~ 20)
  prop <- rbeta(length(prob), shape1 = prob * 20, shape2 = (1 - prob) * 20)
  
  tibble(pid = factor(pid), condition, prop)
}

library(brms)
library(bayestestR)

run_sim <- function(n, effect = 0.3) {
  dat <- simulate_data(n, effect)
  
  fit <- suppressMessages(brm(
    prop ~ condition + (1 | pid),
    data = dat,
    family = Beta(),
    prior = prior(normal(0, 1), class = "b"),
    chains = 2, iter = 2000, refresh = 0, silent = TRUE
  ))
  
  # Extract pd (probability of direction) for the stable effect
  result <- describe_posterior(fit, test = "p_direction")
  pd_stable <- result[result$Parameter == "conditionstable", "pd"]
  return(pd_stable > 0.95)
}


n_reps <- 20  # start small for testing; increase to 500+ for final estimate
power_estimates <- replicate(n_reps, run_sim(n = 60))  # for N=60

mean(power_estimates)  # estimated power


# Assuming your brms model is called beta_model
# Simulate new data from the posterior

library(brms)

# Simulate a dataset with the same structure
sim_data <- posterior_predict(beta_model, draws = 1, summary = FALSE)
sim_data <- as.data.frame(t(sim_data))  # rows = observations, cols = iterations
sim_data$prop <- sim_data[[1]]  # take one simulation for simplicity

# Use your original data structure and replace the outcome
sim_df <- model.data.prop  # your original data
sim_df$prop <- sim_data$prop

library(lme4)

# Fit a simplified linear mixed model (LMM) as an approximation
lmm_model <- lmer(prop ~ condition + age + (1 | pid), data = sim_df)

library(simr)

# Extend the model to simulate power for a larger sample (e.g., 80 participants)
extended_model <- extend(lmm_model, along = "pid", n = 80)

# Run power simulation
powerSim(extended_model, fixed("condition", "anova"), nsim = 100)

pc <- powerCurve(extended_model, along = "pid", breaks = seq(30, 100, by = 10), fixed("condition", "anova"))
plot(pc)

# Simulate power for the fixed effect of condition
# This will test whether the *condition* effect (overall) is significant
power_result <- powerSim(extended_model, fixed("condition", "anova"), nsim = 100)
print(power_result)



# Fit reduced model with only condition (drop age and random slopes)
beta_model_cond <- brm(
  formula = prop ~ condition + (1 | pid),
  family = zero_one_inflated_beta(),
  data = model.data.prop,
  chains = 4, iter = 4000, warmup = 1000, cores = 4
)

# Now compute R² for this simpler model
bayes_R2(beta_model_cond)

install.packages("pwr")  # run if not already installed
library(pwr)

# Run the power analysis
pwr_result <- pwr.anova.test(
  k = 3,         # number of groups/conditions
  f = 0.29,      # Cohen's f from your delta R²
  sig.level = 0.05,  # alpha
  power = 0.90       # desired power
)

# View required sample size
print(pwr_result)



# probability of direction
describe_posterior(
  beta_model,
  effects = "fixed",
  component = "all",
  test=c("p_direction"), 
  centrality = "all")

describe_posterior(beta_model, effects = "fixed", component = "mu", standardize = "refit")


pred <- avg_predictions(beta_model, variables = "condition") %>%  posterior_draws() %>%
  mutate(condition = factor(condition, levels = c("unnaturalv", 
                                                  "naturalv", 
                                                  "stable")))

ggplot(pred, aes(x = draw, y = condition, fill = condition)) +
  geom_vline(xintercept = 0.5, color = "gray", linetype = "dashed") + 
  stat_halfeye(.width=c(.95), alpha = 0.8)  +
  scale_fill_manual(values = condition.colors) + 
  scale_y_discrete(labels = c("atypical\nvariation", "typical\nvariation", "stable")) + 
  guides(fill = "none") +
  labs(x ="Accuracy\n(Proportion looking to target object)", y = "Condition",
       caption = "95% credible intervals shown in black") +
  theme_test(base_size=20) 
ggsave("plots/overall-condition-bayes-BU.png")


effects_age <- conditional_effects(beta_model, effects = "age")$age

effects_condition <- conditional_effects(beta_model, effects = "condition")$condition

# Plot the effects

effects_condition %>%
  mutate(condition = factor(condition, levels = c("stable", "naturalv", 
                                                  "unnaturalv"))) %>%
  ggplot(aes(x = condition, y = `estimate__`, color = condition)) +
  geom_hline(yintercept = 0.5, color = "gray", linetype = "dashed", size = 1.25) +
  geom_point(size = 4) + 
  geom_errorbar(aes(ymin = `lower__`, ymax = `upper__`), width = 0.2, size = 1.25) + 
  scale_color_manual(values = condition.colors) + 
  scale_x_discrete(labels = c("consistent", "typical\nvariation", "atypical\nvariation")) +
  scale_y_continuous(limits = c(0.25, 0.75), breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
  guides(color = "none") +
  labs(y ="Model-predicted accuracy\n(Proportion looking to target object)", x = "Condition",
       caption = "Error bars represent 95% credible intervals") +
  theme_test(base_size=20) + 
  theme(legend.position = "none", 
        title = element_text(face = "bold"))
  
ggsave("pres/bayes-condition-BU.png", height = 6, width = 8)

effects_condition <- effects_condition %>%
  mutate(condition = factor(condition,
                            levels = c("stable", "naturalv", "unnaturalv")))

ggplot(bypid.condition.summary, aes(x = condition, y = accuracy, 
                                    color = condition, fill = condition)) + 
  geom_hline(yintercept = 0.5, color = "gray", linetype = "dashed", size = 1.25) +
  ggrain::geom_rain(boxplot.args = list(fill = NA, color = "black", outlier.shape = NA),
                    violin.args = list(alpha = 0.5, color = "black"),
                    point.args = list(color = "black", shape = 21)) +
  geom_errorbar(data = effects_condition, 
                mapping = aes(x = as.numeric(condition) - 0.2, y = `estimate__`, ymin = `lower__`, ymax = `upper__`), 
                width = 0.2, size = 1.5) + 
  geom_point(data = effects_condition, 
                mapping = aes(x = as.numeric(condition) - 0.2, y = `estimate__`), 
                size = 3) + 
  labs(x = "Condition", y = "Accuracy\n(Proportion looking to target object)") + 
  scale_color_manual(values = condition.colors) + 
  scale_fill_manual(values = condition.colors) + 
  scale_x_discrete(labels = c("consistent", "typical\nvariation", "atypical\nvariation")) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
  theme_classic(base_size = 20) + 
  theme(legend.position = "none", 
        title = element_text(face = "bold"))

ggsave("pres/BU2-NEW.png", height = 6, width = 8)


mean_age <- select(data, pid, age_m) %>%
  distinct() %>% 
  summarize(mean = mean(age_m)) %>%
  pull(mean)

effects_age %>%
  ggplot(aes(x = age, y = `estimate__`)) +
  geom_hline(yintercept = 0.5, color = "gray", linetype = "dashed", size = 1.25) +
  geom_line(size = 1.25) + 
  geom_ribbon(aes(ymin = `lower__`, ymax = `upper__`), alpha = 0.25) + 
  scale_x_continuous(labels = function(x) round(x + mean_age + 1)) + 
  scale_y_continuous(limits = c(0.25, 0.75), breaks = c(0.25, 0.5, 0.75)) + 
  labs(y ="Model-predicted accuracy\n(Proportion looking to target object)", x = "Age (months)",
       caption = "Shaded region represents 95% credible interval") +
  theme_test(base_size=20) + 
  theme(legend.position = "none", 
        title = element_text(face = "bold"))
ggsave("pres/bayes-age-NEW.png", height = 6, width = 8)


# new_data <- expand.grid(
#   condition = unique(model.data.prop$condition), # Replace 'model.data.prop' with your actual data frame name
#   age = seq(min(model.data.prop$age), max(model.data.prop$age), length.out = 100) # Adjust seq() as needed
# )
# 
# new_data$predicted <- posterior_epred(beta_model, newdata = new_data) %>% 
#   rowMeans() 

# marginaleffects::avg_slopes(
#   beta_model,
#   variables  = "age", 
#   by="condition") 

ggemmeans(beta_model, terms=c("age", "condition")) %>% 
  plot() + 
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black") + 
  geom_line(size = 1.5) + 
  scale_color_manual(values = condition.colors, labels = c("unnatural variation", "stable", "natural variation")) + 
  scale_fill_manual(values = condition.colors, labels = c("unnatural variation", "stable", "natural variation")) + 
  scale_y_continuous(limits = c(0,1), breaks = c(0, 0.25, 0.5, 0.75, 1)) + 
  labs(x = "Age (months)",  y = "Model-Predicted Accuracy\n(Proportion looking to target object)") +  
  theme_classic(base_size = 15) + 
  scale_x_continuous(labels = function(x) round(x + mean_age))

ggsave("age-effect-bayes.png", width = 8, height = 6)
