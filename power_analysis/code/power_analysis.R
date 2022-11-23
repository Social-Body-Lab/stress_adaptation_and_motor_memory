"
Authors: Arran J. Davis
Emails: arran.davis@anthro.ox.ac.uk | davis.arran@gmail.com
Affiliation: Social Body Lab, Institute of Human Sciences, University of Oxford
Date: 21 November 2022
"

library(extrafont)
library(faux)
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(sjPlot)
library(tidyr)
library(broom.mixed)
library(purrr)
library(censReg)
library(plm)
library(reshape2)
library(data.table)

#clean environment
rm(list = ls())

#set working directory
setwd(getSrcDirectory()[1])
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#load plot themes
source("plot_theme.R")

################################################################################################################################################

### CREATE DATA ###

#create the neighbourhood variable 
between = list(neighbourhood = c(Scampia = "Scampia", 
                                 Roma = "Roma camp",
                                 Pozzuoli = "Pozzuoli "))

#create the test variable
within = list(test = c("Short-term memory in reference to motor stimuli",
                       "Working memory in reference to motor stimuli"))

within_short = list(test = c("motor-short-term",
                             "motor-working"))

#create scores for each neighbourhood for short-term memory in reference to motor stimuli
scampia_motor_stm = 4.65
roma_motor_stm = 4.65
pozzuoli_motor_stm = 5.15

#create scores for each neighbourhood for working memory in reference to motor stimuli
scampia_motor_wm = 3.65
roma_motor_wm = 3.65
pozzuoli_motor_wm = 4.15

#create the mean scores for each neighbourhood on each memory test 
hood_means = list(Scampia = c(motor_stm = scampia_motor_stm, 
                              motor_wm = scampia_motor_wm),
                  Roma = c(motor_stm = roma_motor_stm,
                           motor_wm = roma_motor_wm),
                  Pozzuoli = c(motor_stm = pozzuoli_motor_stm,
                               motor_wm = pozzuoli_motor_wm))

#create the standard deviation scores for each neighbourhood on each memory test 
hood_sds = list(Scampia = c(motor_stm = 1, 
                            motor_wm = 1),
                Roma = c(motor_stm = 1,
                         motor_wm = 1),
                Pozzuoli = c(motor_stm = 1,
                             motor_wm = 0.75))

#set the correlation of test score results for each neighbourhood
hood_cors = list(Scampia = .6, Roma = .6, Pozzuoli = .7)

#create the dataframe
dat = sim_design(within_short, between, n = 50, 
                 mu = hood_means, sd = hood_sds, r = hood_cors,
                 empirical = FALSE, plot = FALSE)

#make the data long format
dat_long = melt(setDT(dat), id.vars = c("id","neighbourhood"), variable.name = "test")

#ensure no values are below 2 or above 8
dat_long$value = ifelse(dat_long$value < 2, 2,
                        ifelse(dat_long$value > 8, 8, dat_long$value))

#set contrasts of neighbourhood so that "Pozzuoli" is the reference
contrasts(dat_long$neighbourhood) = contr.treatment(3, base = 3)

#round the test scores to a whole number
dat_long$value = round(dat_long$value, 0)

### ### ###

#plot the data and save it
test_labels = c("Short-term memory in\nreference to motor stimuli",
                "Working memory in\nreference to motor stimuli")

hood_colors = c("#999999", "#009E73", "#F0E442")

scores_by_hood = ggplot(dat_long, aes(x = test, y = value)) + 
                  geom_violin(aes(fill = neighbourhood), trim = TRUE, position = position_dodge(0.9)) +
                  geom_boxplot(aes(fill = neighbourhood), width = 0.15, position = position_dodge(0.9)) +
                  scale_x_discrete(labels = test_labels) +
                  scale_y_continuous(limits = c(0,8), breaks = c(seq(0,8,1))) +
                  scale_fill_manual(values = hood_colors) +
                  ylab("Test score") +
                  xlab("Test type") +
                  labs(fill = "Neighbourhood") +
                  aes(ymin = 0) +
                  avenir_theme

ggsave("../plots/plotted_simulated_data_study_3.jpg", scores_by_hood, width = 10, height = 5)

################################################################################################################################################

### ANALYSES ###

#create a dataset for each test type
motor_stm_df = droplevels(subset(dat_long, dat_long$test == "motor-short-term"))
motor_wm_df = droplevels(subset(dat_long, dat_long$test == "motor-working"))

### ### ###

#run models for each test type (first set the contrast for neighbourhood)
contrasts(motor_stm_df$neighbourhood) = contr.treatment(3, base = 3)
motor_stm_mod = lm(value ~ neighbourhood, data = motor_stm_df)
summary(motor_stm_mod)

contrasts(motor_wm_df$neighbourhood) = contr.treatment(3, base = 3)
motor_wm_mod = lm(value ~ neighbourhood, data = motor_wm_df)
summary(motor_wm_mod)

### ### ###

library(effsize)

#create groups for each neighbourhood and test
roma_motor_stm = subset(motor_stm_df, motor_stm_df$neighbourhood == "Roma")
roma_motor_stm_list = as.numeric(roma_motor_stm$value)

pozzuoli_motor_stm = subset(motor_stm_df, motor_stm_df$neighbourhood == "Pozzuoli")
pozzuoli_motor_stm_list = as.numeric(pozzuoli_motor_stm$value)

#get Cohen's d for the comparisons between groups
cohen.d(roma_motor_stm_list, pozzuoli_motor_stm_list) 

################################################################################################################################################

### DATA SIMULATION FOR POWER ANALYSIS ###

#create a function for power simulations (based on above data creation procedure); `participant_n` is the number of participants per neighbourhood
power_simulation = function(participant_n = 40,
                            pozzuoli_motor_stm = 5.15,
                            pozzuoli_motor_wm = 4.15,
                            high_stress_motor_effect = -0.3,
                            pozzuoli_stm_sd = 1,
                            pozzuoli_wm_sd = 0.75,
                            high_stress_stm_sd = 1,
                            high_stress_wm_sd = 1,
                            test_score_cors_pozzuoli = 0.7,
                            test_score_cors_high_stress = 0.6,
                            ...) {
  
  #create the neighbourhood variable 
  between = list(neighbourhood = c(Scampia = "Scampia", 
                                   Roma = "Roma camp",
                                   Pozzuoli = "Pozzuoli "))
  
  #create the test variable
  within = list(test = c("Short-term memory in reference to motor stimuli",
                         "Working memory in reference to motor stimuli"))
  
  within_short = list(test = c("motor-short-term",
                               "motor-working"))
  
  #create scores for each neighbourhood for social short-term memory
  scampia_motor_stm = pozzuoli_motor_stm + high_stress_motor_effect
  roma_motor_stm = pozzuoli_motor_stm + high_stress_motor_effect
  pozzuoli_motor_stm = pozzuoli_motor_stm
  
  #create scores for each neighbourhood for verbal working memory
  scampia_motor_wm = pozzuoli_motor_wm + high_stress_motor_effect
  roma_motor_wm = pozzuoli_motor_wm + high_stress_motor_effect
  pozzuoli_motor_wm = pozzuoli_motor_wm
  
  #create the mean scores for each neighbourhood on each memory test 
  hood_means = list(Scampia = c(motor_stm = scampia_motor_stm, 
                                motor_wm = scampia_motor_wm),
                    Roma = c(motor_stm = roma_motor_stm,
                             motor_wm = roma_motor_wm),
                    Pozzuoli = c(motor_stm = pozzuoli_motor_stm,
                                 motor_wm = pozzuoli_motor_wm))
  
  #create the standard deviation scores for each neighbourhood on each memory test 
  hood_sds = list(Scampia = c(motor_stm = high_stress_stm_sd, 
                              motor_wm = high_stress_wm_sd),
                  Roma = c(motor_stm = high_stress_stm_sd,
                           motor_wm = high_stress_wm_sd),
                  Pozzuoli = c(motor_stm = pozzuoli_stm_sd,
                               motor_wm = pozzuoli_wm_sd))
  
  #set the correlation of test score results for each neighbourhood
  hood_cors = list(Scampia = test_score_cors_high_stress, Roma = test_score_cors_high_stress, Pozzuoli = test_score_cors_pozzuoli)
  
  #create the dataframe
  dat = sim_design(within_short, between, n = participant_n, 
                   mu = hood_means, sd = hood_sds, r = hood_cors,
                   empirical = FALSE, plot = FALSE)
  
  #make the data long format
  dat_long = melt(setDT(dat), id.vars = c("id","neighbourhood"), variable.name = "test")
  
  #ensure no values are below 2 or above 8
  dat_long$value = ifelse(dat_long$value < 2, 2,
                          ifelse(dat_long$value > 8, 8, dat_long$value))
  
  #set contrasts of neighbourhood so that "Pozzuoli" is the reference
  contrasts(dat_long$neighbourhood) = contr.treatment(3, base = 3)
  
  #round the test scores to a whole number
  dat_long$value = round(dat_long$value, 0)  
  
  ### ### ###

  #create a dataset for each test type
  motor_stm_df = droplevels(subset(dat_long, dat_long$test == "motor-short-term"))
  motor_wm_df = droplevels(subset(dat_long, dat_long$test == "motor-working"))

  ### ### ###
  
  #set the contrast for neighbourhood for each dataset
  contrasts(motor_stm_df$neighbourhood) = contr.treatment(3, base = 3)
  contrasts(motor_wm_df$neighbourhood) = contr.treatment(3, base = 3)
  
  #run model for each test type
  motor_stm_mod = lm(value ~ neighbourhood, data = motor_stm_df)
  motor_wm_mod = lm(value ~ neighbourhood, data = motor_wm_df)
  
  ### ### ###
  
  #create groups for one of the high-stress neighbourhoods (they are drawn from the same population) and calculate effect sizes for each outcome
  roma_motor_stm = subset(motor_stm_df, motor_stm_df$neighbourhood == "Roma")
  roma_motor_stm_list = as.numeric(roma_motor_stm$value)
  pozzuoli_motor_stm = subset(motor_stm_df, motor_stm_df$neighbourhood == "Pozzuoli")
  pozzuoli_motor_stm_list = as.numeric(pozzuoli_motor_stm$value)
  
  roma_motor_wm = subset(motor_wm_df, motor_wm_df$neighbourhood == "Roma")
  roma_motor_wm_list = as.numeric(roma_motor_wm$value)
  pozzuoli_motor_wm = subset(motor_wm_df, motor_wm_df$neighbourhood == "Pozzuoli")
  pozzuoli_motor_wm_list = as.numeric(pozzuoli_motor_wm$value)
  
  #get Cohen's d for the comparisons between groups
  effect_size_results_motor_stm = cohen.d(roma_motor_stm_list, pozzuoli_motor_stm_list)
  effect_size_motor_stm = as.numeric(effect_size_results_motor_stm$estimate)
  effect_size_motor_stm_lower = as.numeric(effect_size_results_motor_stm$conf.int[1])
  effect_size_motor_stm_upper = as.numeric(effect_size_results_motor_stm$conf.int[2])
  
  effect_size_results_motor_wm = cohen.d(roma_motor_wm_list, pozzuoli_motor_wm_list)
  effect_size_motor_wm = as.numeric(effect_size_results_motor_wm$estimate)
  effect_size_motor_wm_lower = as.numeric(effect_size_results_motor_wm$conf.int[1])
  effect_size_motor_wm_upper = as.numeric(effect_size_results_motor_wm$conf.int[2])
  
  ### ### ###

  #return a dataframe of the model results 
  motor_stm_results = broom.mixed::tidy(motor_stm_mod)
  motor_stm_results$outcome = "Short-term memory in reference to motor stimuli"
  motor_stm_results$high_stress_effect = c(NA, effect_size_motor_stm, NA)
  motor_stm_results$high_stress_effect_lower = c(NA, effect_size_motor_stm_lower, NA)
  motor_stm_results$high_stress_effect_upper = c(NA, effect_size_motor_stm_upper, NA)
  
  motor_wm_results = broom.mixed::tidy(motor_wm_mod)
  motor_wm_results$outcome = "Working memory in reference to motor stimuli"
  motor_wm_results$high_stress_effect = c(NA, effect_size_motor_wm, NA)
  motor_wm_results$high_stress_effect_lower = c(NA, effect_size_motor_wm_lower, NA)
  motor_wm_results$high_stress_effect_upper = c(NA, effect_size_motor_wm_upper, NA)
  
  do.call("rbind", list(motor_stm_results, motor_wm_results))

}

#run repeated simulations with dataset variants for each effect size
simulations = crossing(replications = 1:1000,
                       participant_n = c(50, 75, 100, 125, 150),
                       pozzuoli_motor_stm = 5.15,
                       pozzuoli_motor_wm = 4.15,
                       high_stress_motor_effect = c(-0.1, -0.2, -0.3, -0.4, -0.5),
                       pozzuoli_stm_sd = 1,
                       pozzuoli_wm_sd = 0.75,
                       high_stress_stm_sd = 1,
                       high_stress_wm_sd = 1,
                       test_score_cors_pozzuoli = 0.7,
                       test_score_cors_high_stress = 0.6) %>% mutate(analysis = pmap(., power_simulation)) %>% unnest(analysis)

#save the dataframes
write.csv(simulations, "../data/study3_motor_memory_data_simulations.csv")

################################################################################################################################################

### POWER ANALYSIS FOR BETWEEN GROUP COMPARISONS ###

#subset the data to each outcome type
stm_sims = subset(simulations, simulations$outcome == "Short-term memory in reference to motor stimuli")
wm_sims = subset(simulations, simulations$outcome == "Working memory in reference to motor stimuli")

#create dataset to plot power analysis for neighbourhood main effects (both high-stress neighbourhoods were estimated to have the same effect)
simulation_results_neighbourhood_stm = filter(stm_sims, term == "neighbourhood1") %>%
                                            group_by(high_stress_motor_effect, participant_n) %>% 
                                            summarise(power = mean(p.value < .05), .groups = "drop")


simulation_results_neighbourhood_wm = filter(wm_sims, term == "neighbourhood1") %>%
                                           group_by(high_stress_motor_effect, participant_n) %>% 
                                           summarise(power = mean(p.value < .05), .groups = "drop")

### ### ###

#plot power analysis for neighbourhood effect on short-term memory in reference to motor stimuli
stm_pa = ggplot(aes(as.character(high_stress_motor_effect), participant_n, fill = power), 
                data = simulation_results_neighbourhood_stm) +
          geom_tile() +
          geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
          scale_y_continuous(breaks = c(50, 75, 100, 125, 150)) +
          scale_fill_viridis_c(name = "Power",
                               limits = c(0, 1), 
                               breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                               labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
          xlab("High-stress neighbourhood effect on short-term memory in reference to motor stimuli") + 
          ylab("Participant sample size (per neighbourhood)") +
          avenir_theme

ggsave("../plots/power_analysis_motor_short_term_memory.jpg", stm_pa, width = 10, height = 5)

### ### ###

#plot power analysis for neighbourhood effect on working memory in reference to motor stimuli
verb_pa = ggplot(aes(as.character(high_stress_motor_effect), participant_n, fill = power), 
                 data = simulation_results_neighbourhood_wm) +
           geom_tile() +
           geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
           scale_y_continuous(breaks = c(50, 75, 100, 125, 150)) +
           scale_fill_viridis_c(name = "Power",
                                limits = c(0, 1), 
                                breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
           xlab("High-stress neighbourhood effect on working memory in reference to motor stimuli") + 
           ylab("Participant sample size (per neighbourhood)") +
           avenir_theme

ggsave("../plots/power_analysis_motor_working_memory.jpg", verb_pa, width = 10, height = 5)

################################################################################################################################################

### EFFECT SIZE ESTIMATES ###

#create datasets
simulation_results_effect_sizes_stm = stm_sims[complete.cases(stm_sims), ]
simulation_results_effect_sizes_wm = wm_sims[complete.cases(wm_sims), ]

#select a participant sample size (per neighbourhood)
neighbourhood_n = 125
effect = -0.4

#get the mean effect size and CI for the assumed effect and a sample size of 125 for short-term memory (to estimate comparing effect sizes)
assumed_stm_effect = subset(simulation_results_effect_sizes_stm,
                            simulation_results_effect_sizes_stm$participant_n == neighbourhood_n & 
                            simulation_results_effect_sizes_stm$high_stress_motor_effect == effect)

assumed_wm_effect = subset(simulation_results_effect_sizes_wm,
                           simulation_results_effect_sizes_wm$participant_n == neighbourhood_n & 
                           simulation_results_effect_sizes_wm$high_stress_motor_effect == effect)


combined = data.frame(stm_lower = assumed_stm_effect$high_stress_effect_lower,
                      stm_effect = assumed_stm_effect$high_stress_effect,
                      stm_upper = assumed_stm_effect$high_stress_effect_upper,
                      wm_lower = assumed_wm_effect$high_stress_effect_lower,
                      wm_effect = assumed_wm_effect$high_stress_effect,
                      wm_upper = assumed_wm_effect$high_stress_effect_upper)

#add a variable that is whether or not the confidence intervals overlap
combined$overlapping_CI = ifelse(combined$stm_lower <= combined$stm_upper, TRUE, FALSE)
print(paste0("PERCENTAGE OF SIGNIFICANT EFFECT SIZE DIFFERENCES: ",
      table(combined$overlapping_CI)[1] / sum(table(combined$overlapping_CI)) * 100, "%"))

print(paste0("MEAN EFEFCT SIZES: ", mean(assumed_stm_effect$high_stress_effect)))
print(paste0("MEAN EFEFCT SIZES: ", mean(assumed_wm_effect$high_stress_effect)))

################################################################################################################################################

### COMPARISON WITH STUDY 1 TEST SCORES (RELATIVE CHANGE) - EXAMPLE DATA ###

#set difference scores, standard deviations, and correlations
high_stress_difference_standard = -0.5
high_stress_difference_motor = 5
sd = 1
test_score_cors_pozzuoli = 0.7
test_score_cors_high_stress = 0.6

#get high and low stress neighbourhood means for both test types
high_stress_standard_mean = high_stress_difference_standard / 3
low_stress_standard_mean = high_stress_standard_mean - high_stress_difference_standard

high_stress_motor_mean = high_stress_difference_motor / 3
low_stress_motor_mean = high_stress_motor_mean - high_stress_difference_motor  

#set expected means and SDs (data will be scaled within test type) for each group and test
pozzuoli_standard_mean = low_stress_standard_mean
pozzuoli_standard_sd = sd
roma_standard_mean = high_stress_standard_mean
roma_standard_sd = sd
scampia_standard_mean = high_stress_standard_mean
scampia_standard_sd = sd

pozzuoli_motor_mean = low_stress_motor_mean
pozzuoli_motor_sd = sd
roma_motor_mean = high_stress_motor_mean
roma_motor_sd = sd
scampia_motor_mean = high_stress_motor_mean
scampia_motor_sd = sd

#create a dataframe of just 20 participants, some with two and one measures (for illustrative purposes)
data = data.frame(test_score = c(rnorm(20, pozzuoli_standard_mean, pozzuoli_standard_sd),
                                 rnorm(19, roma_standard_mean, roma_standard_sd),
                                 rnorm(18, scampia_standard_mean, scampia_standard_sd),
                                 rnorm(20, pozzuoli_standard_mean, pozzuoli_standard_sd),
                                 rnorm(19, roma_standard_mean, roma_standard_sd),
                                 rnorm(18, scampia_standard_mean, scampia_standard_sd),
                                 rnorm(17, pozzuoli_motor_mean, pozzuoli_motor_sd),
                                 rnorm(18, roma_motor_mean, roma_motor_sd),
                                 rnorm(16, scampia_motor_mean, scampia_motor_sd)),
                  test_type = c(rep("visuospatial", 20),
                                rep("visuospatial", 19),
                                rep("visuospatial", 18),
                                rep("verbal", 20),
                                rep("verbal", 19),
                                rep("verbal", 18),
                                rep("motor", 17),
                                rep("motor", 18),
                                rep("motor", 16)),
                  hood =  c(rep("Pozzuoli", 20),
                            rep("Roma", 19),
                            rep("Scampia", 18),
                            rep("Pozzuoli", 17),
                            rep("Roma", 18),
                            rep("Scampia", 16)),
                  participant_id = c(seq(1, 20, 1),
                                     seq(21, 24, 1), seq(26, 40, 1),
                                     41, 42, 43, seq(45, 56, 1), 58, 59, 60,
                                     seq(1, 20, 1),
                                     seq(21, 24, 1), seq(26, 40, 1),
                                     41, 42, 43, seq(45, 56, 1), 58, 59, 60,
                                     seq(1, 7, 1), seq(9, 13, 1), seq(15, 19, 1),
                                     21, 22, seq(24, 37, 1), 39, 40,
                                     42, 44, seq(45, 53, 1), 55, 56, 57, 59, 60)) 

#make relevant variables factors
data$hood = factor(data$hood, levels = c("Scampia", "Roma", "Pozzuoli"))
data$test_type = factor(data$test_type, levels = c("visuospatial", "verbal", "motor"))

#run the model (first set contrasts)
contrasts(data$hood) = contr.treatment(3, base = 3)
contrasts(data$test_type) = contr.treatment(3, base = 1)
comparison_model = lmer(test_score ~ hood*test_type + (1 | participant_id), data = data,
                        control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
#check results
summary(comparison_model)

### ### ###

#plot the data
scores_by_hood_and_test = ggplot(data, aes(x = test_type, y = test_score)) + 
                            geom_violin(aes(fill = hood), trim = TRUE, position = position_dodge(0.9)) +
                            geom_boxplot(aes(fill = hood), width = 0.15, position = position_dodge(0.9)) +
                            scale_fill_manual(values = hood_colors) +
                            scale_x_discrete(labels = c("Visuospatial\nshort-term memory",
                                                        "Verbal\nshort-term memory",
                                                        "Short-term memory in\nreference to motor stimuli")) +
                            ylab("Test score") +
                            xlab("Test type") +
                            labs(fill = "Neighbourhood") +
                            aes(ymin = 0) +
                            avenir_theme

ggsave("../plots/plotted_simulated_data_for_test_type_comparisons.jpg", scores_by_hood_and_test, width = 10, height = 5)

#run the model (first set contrasts)
contrasts(data$hood) = contr.treatment(3, base = 3)
contrasts(data$test_type) = contr.treatment(2, base = 1)

comparison_model = lmer(test_score ~ hood*test_type + (1 | participant_id), data = data)
summary(comparison_model)

### ### ###

library(emmeans)
  
#compare results across neighbourhood
em_dat = emmeans(comparison_model, data = data, pairwise ~ test_type | hood)

#check results
emmeans(comparison_model, data = data, pairwise ~ test_type | hood)

################################################################################################################################################

### COMPARISON WITH STUDY 1 TEST SCORES (RELATIVE CHANGE) - EXAMPLE DATA ###

#create a function for power simulations (based on above data creation procedure); assuming less than 20 participants have data for both measures
power_simulation_test_comparison = function(participant_n = 30,
                                            difference_between_hoods = 1,
                                            high_stress_performs_better_on = "motor",
                                            sd = 1,
                                            test_score_cors_pozzuoli = 0.7,
                                            test_score_cors_high_stress = 0.6,
                                            ...) {
  
  #set mean differences
  if (high_stress_performs_better_on == "motor"){
    high_stress_difference_standard = difference_between_hoods * -1
    high_stress_difference_motor = difference_between_hoods
  } else {
    high_stress_difference_standard = difference_between_hoods
    high_stress_difference_motor = difference_between_hoods * -1
  }
  
  #get high and low stress neighbourhood means for both test types
  high_stress_standard_mean = high_stress_difference_standard / 3
  low_stress_standard_mean = high_stress_standard_mean - high_stress_difference_standard

  high_stress_motor_mean = high_stress_difference_motor / 3
  low_stress_motor_mean = high_stress_motor_mean - high_stress_difference_motor  
  
  #set expected means and SDs (data will be scaled within test type) for each group and test
  pozzuoli_standard_mean = low_stress_standard_mean
  pozzuoli_standard_sd = sd
  roma_standard_mean = high_stress_standard_mean
  roma_standard_sd = sd
  scampia_standard_mean = high_stress_standard_mean
  scampia_standard_sd = sd
  
  pozzuoli_motor_mean = low_stress_motor_mean
  pozzuoli_motor_sd = sd
  roma_motor_mean = high_stress_motor_mean
  roma_motor_sd = sd
  scampia_motor_mean = high_stress_motor_mean
  scampia_motor_sd = sd
  
  #create a dataframe
  data = data.frame(test_score = c(rnorm(participant_n, pozzuoli_standard_mean, pozzuoli_standard_sd),
                                   rnorm(participant_n, roma_standard_mean, pozzuoli_standard_sd),
                                   rnorm(participant_n, scampia_standard_mean, scampia_standard_sd),
                                   rnorm(participant_n, pozzuoli_standard_mean, pozzuoli_standard_sd),
                                   rnorm(participant_n, roma_standard_mean, pozzuoli_standard_sd),
                                   rnorm(participant_n, scampia_standard_mean, scampia_standard_sd),
                                   rnorm(participant_n, pozzuoli_motor_mean, pozzuoli_motor_sd),
                                   rnorm(participant_n, roma_motor_mean, roma_motor_sd),
                                   rnorm(participant_n, scampia_motor_mean, scampia_motor_sd)),
                    test_type = c(rep("visuospatial", participant_n),
                                  rep("visuospatial", participant_n),
                                  rep("visuospatial", participant_n),
                                  rep("verbal", participant_n),
                                  rep("verbal", participant_n),
                                  rep("verbal", participant_n),
                                  rep("motor", participant_n),
                                  rep("motor", participant_n),
                                  rep("motor", participant_n)),
                    hood =  c(rep("Pozzuoli", participant_n),
                              rep("Roma", participant_n),
                              rep("Scampia", participant_n),
                              rep("Pozzuoli", participant_n),
                              rep("Roma", participant_n),
                              rep("Scampia", participant_n),
                              rep("Pozzuoli", participant_n),
                              rep("Roma", participant_n),
                              rep("Scampia", participant_n)),
                    participant_id = c(seq(1, participant_n, 1),
                                       seq(participant_n + 1, participant_n * 2, 1),
                                       seq((participant_n * 2) + 1, participant_n * 3, 1)))

  #make relevant variables factors
  data$hood = factor(data$hood, levels = c("Scampia", "Roma", "Pozzuoli"))
  data$test_type = factor(data$test_type, levels = c("visuospatial", "verbal", "motor"))
  
  #run the model (first set contrasts)
  contrasts(data$hood) = contr.treatment(3, base = 3)
  contrasts(data$test_type) = contr.treatment(3, base = 1)
  comparison_model = lmer(test_score ~ hood*test_type + (1 | participant_id), data = data,
                          control=lmerControl(check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))

  #compare results across neighbourhood
  em_dat = emmeans(comparison_model, data = data, pairwise ~ test_type | hood)
  
  #extract results (visuospatial - motor test comparison for Scampia and Pozzuoli)
  high_stress_test_comparison = as.data.frame(em_dat$contrasts[2])
  low_stress_test_comparison = as.data.frame(em_dat$contrasts[8])
  
  do.call("rbind", list(high_stress_test_comparison, low_stress_test_comparison))
  
}

#run repeated simulations with dataset variants for each effect size
simulations = crossing(simulations = 1:1000,
                       participant_n = c(25, 50, 75, 100, 125),
                       difference_between_hoods = c(1.25, 1, 0.75, 0.5, 0.25),
                       high_stress_performs_better_on = "motor",
                       sd = 1,
                       test_score_cors_pozzuoli = 0.7,
                       test_score_cors_high_stress = 0.6) %>% mutate(analysis = pmap(., power_simulation_test_comparison)) %>% unnest(analysis)

#save the dataframe
write.csv(simulations, "../data/study3_motor_and_standard_memory_comparisons_data_simulations.csv")

################################################################################################################################################

### POWER ANALYSIS FOR COMPARISONS WITH STUDY 1 TEST SCORES (RELATIVE CHANGE) ###

#create dataset to plot power analysis (both high-stress neighbourhoods were estimated to have the same effect)
simulation_results_roma = filter(simulations, hood == "Scampia") %>%
                                    group_by(difference_between_hoods, participant_n) %>% 
                                    summarise(power = mean(p.value < .05), .groups = "drop")

### ### ###

#plot power analysis for test score differences for the Scampia neighbourhood
stm_comp_pa = ggplot(aes(as.character(difference_between_hoods), participant_n, fill = power), 
                     data = simulation_results_roma) +
                geom_tile() +
                geom_text(aes(label = sprintf("%.2f", power)), color = "black", size = 5) +
                scale_y_continuous(breaks = c(25, 50, 75, 100, 125, 150)) +
                scale_fill_viridis_c(name = "Power",
                                     limits = c(0, 1), 
                                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                                     labels = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                xlab("Difference in standardised 'traditional' memory test scores and motor memory test scores") + 
                ylab("Number of participants from the neighbourhood") +
                avenir_theme

ggsave("../plots/power_analysis_traditional_and_motor_memory_comparison.jpg", stm_comp_pa, width = 10, height = 5)
