# Competitive ability predicts the spread of novel carnivory in California ground squirrels 


# 1) Load libraries -------------------------------------------------------

# install devtools if not already
if (!require("devtools")) install.packages("devtools")
devtools::install_github("michaelchimento/STbayes")
library(brms)
library(car)
library(dplyr)
library(flextable)
library(ggplot2)

library(STbayes)
library(tidyr)
library(dplyr)
library(posterior)
library(officer)
library(emmeans)
library(performance)
library(stringr)
library(igraph)
library(grid)
library(gridGraphics)
library(patchwork)
library(ggplotify)
library(qgraph)



# 2) Load data and create summary tables ------------------------------------------------------------

# create a summary of the genome data
genome_data <- read.csv("Input Data/genome_summary.csv", row.names=1)

# how many samples
length(genome_data$mark)
# 100 samples

# by how many different individuals
length(unique(genome_data$mark))
# 56

total_samples_with_dna <- genome_data %>%
  mutate(
    any_animal_dna = if_any(
      Microtus.californicus:Peromyscus.truei,
      ~ . > 0
    )
  ) %>%
  summarise(
    total = sum(any_animal_dna),
    
    adults = sum(any_animal_dna & age == "A"),
    juveniles = sum(any_animal_dna & age == "P"),
    
    females = sum(any_animal_dna & sex == "F"),
    males = sum(any_animal_dna & sex == "M")
  )

total_samples_with_dna
# total adults juveniles females males
# 1    35     18        17      21    14

# 35 samples with DNA


# how many of which species
species_summary <- genome_data %>%
  summarise(across(Microtus.californicus:Peromyscus.truei,
                   ~ mean(. > 0))) %>%
  tidyr::pivot_longer(everything(),
                      names_to = "species",
                      values_to = "percent_samples")

# # A tibble: 5 × 2
# species                   percent_samples
# <chr>                               <dbl>
#   1 Microtus.californicus                0.22
# 2 Thomomys.bottae                      0.11
# 3 Reithrodontomys.megalotis            0.02
# 4 Mus.musculus                         0.01
# 5 Peromyscus.truei                     0.01

prop_carnivory <- genome_data %>%
  mutate(
    any_animal_dna = if_any(
      Microtus.californicus:Peromyscus.truei,
      ~ . > 0
    )
  ) %>%
  group_by(age, sex) %>%
  summarise(
    n_total = n(),
    n_positive = sum(any_animal_dna),
    prop_positive = n_positive / n_total,
    .groups = "drop"
  )

prop_carnivory


# this shows evidence for carniovry from molecular data across age-sex classes

# # A tibble: 4 × 5
# age   sex   n_total n_positive prop_positive
# <chr> <chr>   <int>      <int>         <dbl>
#   1 A     F          41         13         0.317
# 2 A     M          23          5         0.217
# 3 P     F          19          8         0.421
# 4 P     M          17          9         0.529


str(genome_data)


# Summarise by age and sex
summary_table <- genome_data %>%
  mutate(has_mammal_DNA = rowSums(across(Microtus.californicus:Peromyscus.truei)) > 0) %>%
  group_by(age, sex) %>%
  summarise(
    unique_individuals = n_distinct(mark),
    n_fecal_samples = n(),
    n_mammal_DNA = sum(has_mammal_DNA),
    prevalence = round(n_mammal_DNA / n_fecal_samples, 3),
    .groups = "drop"
  )

# Total row
total_row <- genome_data %>%
  mutate(has_mammal_DNA = rowSums(across(Microtus.californicus:Peromyscus.truei)) > 0) %>%
  summarise(
    age = "Total",
    sex = "",
    unique_individuals = n_distinct(mark),
    n_fecal_samples = n(),
    n_mammal_DNA = sum(has_mammal_DNA),
    prevalence = round(n_mammal_DNA / n_fecal_samples, 3)
  )

# Bind together
final_table <- bind_rows(summary_table, total_row)

# Recode age and sex labels
final_table <- final_table %>%
  mutate(
    age = case_match(age, "A" ~ "Adult", "P" ~ "Juvenile", "Total" ~ "Total"),
    sex = case_match(sex, "M" ~ "Male", "F" ~ "Female", "" ~ "")
  )

print(final_table)


# # A tibble: 5 × 6
# age      sex      unique_individuals n_fecal_samples n_mammal_DNA prevalence
# <chr>    <chr>                 <int>           <int>        <int>      <dbl>
#   1 Adult    "Female"                 18              41           13      0.317
# 2 Adult    "Male"                   15              23            5      0.217
# 3 Juvenile "Female"                 11              19            8      0.421
# 4 Juvenile "Male"                   12              17            9      0.529
# 5 Total    ""                       56             100           35      0.35


# and load diffusion data and individual covariates
event_data <- read.csv("Input Data/event.data.csv")
edge_list <- read.csv("Input Data/edge_list.csv", row.names = 1)
ILV_data <- read.csv("Input Data/ILVs.csv", row.names = 1)
ILV_tv <- read.csv("Input Data/ILV_tv.csv", row.names = 1)


# 3) Modeling frequency of carnivory --------------------------------------


# is the rate of carnivory predicted by individual covariates? we combine both molecular and obsrvational data for this analysis. We need to account for effort (number of fecal samples per individual, plus the number of formal observation days present in the study area)

ILV_data$effort<- ILV_data$n_trap_days+ILV_data$n_days_obs+ILV_data$n_feces


# 3.1) descriptive stats --------------------------------------------------

# some descriptive stats:

summary_table <- ILV_data %>%
  group_by(age, sex) %>%
  summarise(
    n_total = n_distinct(id),
    
    DNA_positive = n_distinct(id[n_mammals_DNA > 0]),
    
    obs_positive = n_distinct(id[n_vole_obs > 0]),
    
    either_positive = n_distinct(
      id[n_mammals_DNA > 0 | n_vole_obs > 0]
    ),
    
    .groups = "drop"
  ) %>%
  mutate(
    age = dplyr::recode(age, A = "Adult", P = "Juvenile"),
    sex = dplyr::recode(sex, F = "Female", M = "Male")
  )

# total across all age-sex classes
total_row <- ILV_data %>%
  summarise(
    age = "Total",
    sex = "",
    
    n_total = n_distinct(id),
    
    DNA_positive = n_distinct(id[n_mammals_DNA > 0]),
    
    obs_positive = n_distinct(id[n_vole_obs > 0]),
    
    either_positive = n_distinct(
      id[n_mammals_DNA > 0 | n_vole_obs > 0]
    )
  )

summary_table <- bind_rows(summary_table, total_row)

summary_table # used for Table S2
# # A tibble: 5 × 6
# age      sex      n_total DNA_positive obs_positive either_positive
# <chr>    <chr>      <int>        <int>        <int>           <int>
#   1 Adult    "Female"      21            8           10              13
# 2 Adult    "Male"        15            4            7               9
# 3 Juvenile "Female"      14            5            6               8
# 4 Juvenile "Male"         9            3            3               3
# 5 Total    ""            59           20           26              33

# for vole carnivory

summary_table <- ILV_data %>%
  mutate(
    age_sex_class = case_when(
      age_sex == "A_F" ~ "Adult female",
      age_sex == "A_M" ~ "Adult male",
      age_sex == "P_F" ~ "Juvenile female",
      age_sex == "P_M" ~ "Juvenile male",
      TRUE ~ NA_character_
    ),
    vole_obs_pos = n_vole_obs > 0,
    microtus_pos = n_microtus_positive > 0,
    either_pos = vole_obs_pos | microtus_pos
  ) %>%
  group_by(age_sex_class) %>%
  summarise(
    n_either = sum(either_pos, na.rm = TRUE),
    n_molecular = sum(microtus_pos, na.rm = TRUE),
    n_observational = sum(vole_obs_pos, na.rm = TRUE),
    .groups = "drop"
  )

summary_table

# summary_table
# # A tibble: 4 × 4
# age_sex_class   n_either n_molecular n_observational
# <chr>              <int>       <int>           <int>
#   1 Adult female          11           3              10
# 2 Adult male             8           2               7
# 3 Juvenile female        7           3               6
# 4 Juvenile male          3           3               3

# 3.2) Run model with interaction ----------------------------------------------------------

# assess VIFs

model_vif <- glm(
  formula = total_carnivory ~ agon_PCA + hum_reactivity + n_vole_burrows + age + sex + offset(log(effort)),
  data = ILV_data,
  family = poisson()
)

vif(model_vif)

# vif(model_vif)
# agon_PCA hum_reactivity n_vole_burrows            age            sex 
# 1.448755       1.054466       1.088135       1.372134       1.177833


# we first fit the full model including interactions
model_carnivory_interaction <- brm(
  formula = total_carnivory ~ agon_PCA * age + hum_reactivity*age + n_vole_burrows  + sex*age + offset(log(effort)),
  data = ILV_data,
  family = poisson(),
  chains = 4,
  cores = 4,
  iter = 4000
)

#save(model_carnivory_interaction, file="output/model_carnivory_interaction.RDA")
load("output/model_carnivory_interaction.RDA")


summary(model_carnivory)

# Family: poisson 
# Links: mu = log 
# Formula: total_carnivory ~ agon_PCA * age + hum_reactivity * age + n_vole_burrows + sex * age + offset(log(effort)) 
# Data: ILV_data (Number of observations: 59) 
# Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept              -2.06      0.22    -2.50    -1.64 1.00     6421     6205
# agon_PCA                0.33      0.16     0.02     0.64 1.00     5989     4995
# ageP                    0.13      0.39    -0.65     0.89 1.00     5529     5535
# hum_reactivity         -1.24      0.75    -2.81     0.14 1.00     5877     5285
# n_vole_burrows         -0.02      0.04    -0.11     0.07 1.00     7217     5477
# sexM                   -0.04      0.30    -0.64     0.55 1.00     6493     5718
# agon_PCA:ageP           0.19      0.28    -0.36     0.73 1.00     6744     5463
# ageP:hum_reactivity     1.98      1.09    -0.11     4.13 1.00     5221     4807
# ageP:sexM              -0.32      0.51    -1.34     0.67 1.00     5892     5501
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

# we can see that we do not find evidence for interactive effects, so we repeat the simpler model with just the main effects (after exporting a table with these results)


tab <- as.data.frame(summary(model_carnivory)$fixed)
tab$Parameter <- rownames(tab)
tab <- tab %>%
  rename(
    Estimate = Estimate,
    CI_low = 'l-95% CI',
    CI_high = 'u-95% CI',
    Rhat = Rhat
  ) %>%
  
  # Relabel parameters
  mutate(Parameter = case_match(Parameter,
                                "Intercept" ~ "Baseline carnivory rate",
                                "agon_PCA" ~ "Agonistic tendency (PCA)",
                                "hum_reactivity" ~ "Reactivity to humans",
                                "n_vole_burrows" ~ "Vole burrow prevalence",
                                "ageP" ~ "Age (juvenile vs adult)",
                                "sexM" ~ "Sex (male vs female)",
                                "agon_PCA:ageP" ~ "Agonistic tendency × Age",
                                "ageP:hum_reactivity" ~ "Reactivity to humans × Age",
                                "ageP:sexM" ~ "Sex × Age"
  )) %>%
  
  # Convert to rate ratios
  mutate(
    Estimate = exp(Estimate),
    CI_low = exp(CI_low),
    CI_high = exp(CI_high)
  ) %>%
  
  # Round
  mutate(across(c(Estimate, CI_low, CI_high, Rhat), ~round(.x, 3))) %>%
  
  select(Parameter, Estimate, CI_low, CI_high, Rhat)

colnames(tab) <- c("Predictor", "Rate ratio", "Lower 95% CI", "Upper 95% CI", "R̂")

ft <- flextable(tab) %>%
  autofit()

doc <- read_docx() %>%
  body_add_flextable(ft)
print(doc, target = "Table/model_table_carnivory_w_interactions.docx")



# 3.3.1) Run model without interaction ------------------------------------------------------

# Fit the main effects model
model_carnivory <- brm(
  formula = total_carnivory ~ agon_PCA + hum_reactivity + n_vole_burrows + age + sex + offset(log(effort)),
  data = ILV_data,
  family = poisson(),
  chains = 4,
  cores = 4,
  iter = 4000
)
#save(model_carnivory, file="output/model_carnivory.RDA")
load("output/model_carnivory.RDA")


summary(model_carnivory)

# Family: poisson 
# Links: mu = log 
# Formula: total_carnivory ~ agon_PCA + hum_reactivity + n_vole_burrows + age + sex + offset(log(effort)) 
# Data: ILV_data (Number of observations: 59) 
# Draws: 4 chains, each with iter = 4000; warmup = 2000; thin = 1;
# total post-warmup draws = 8000
# 
# Regression Coefficients:
#   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
# Intercept         -2.21      0.19    -2.60    -1.84 1.00     8549     6285
# agon_PCA           0.39      0.14     0.13     0.66 1.00     6396     5744
# hum_reactivity    -0.33      0.51    -1.36     0.63 1.00     8702     5853
# n_vole_burrows    -0.00      0.04    -0.09     0.08 1.00     8856     5030
# ageP               0.53      0.26     0.02     1.04 1.00     6733     5749
# sexM              -0.12      0.24    -0.61     0.36 1.00     8062     6254
# 
# Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
# and Tail_ESS are effective sample size measures, and Rhat is the potential
# scale reduction factor on split chains (at convergence, Rhat = 1).

tab <- as.data.frame(summary(model_carnivory)$fixed)
tab$Parameter <- rownames(tab)
tab <- tab %>%
  rename(
    Estimate = Estimate,
    CI_low = 'l-95% CI',
    CI_high = 'u-95% CI',
    Rhat = Rhat
  ) %>%
  
  mutate(Parameter = case_match(Parameter,
                                "Intercept" ~ "Baseline carnivory rate",
                                "agon_PCA" ~ "Agonistic tendency (PCA)",
                                "hum_reactivity" ~ "Reactivity to humans",
                                "n_vole_burrows" ~ "Vole burrow prevalence",
                                "ageP" ~ "Age (juvenile vs adult)",
                                "sexM" ~ "Sex (male vs female)"
  )) %>%
  
  mutate(
    Estimate = exp(Estimate),
    CI_low = exp(CI_low),
    CI_high = exp(CI_high)
  ) %>%
  
  mutate(across(c(Estimate, CI_low, CI_high, Rhat), ~round(.x, 3))) %>%
  
  select(Parameter, Estimate, CI_low, CI_high, Rhat)

colnames(tab) <- c("Predictor", "Rate ratio", "Lower 95% CI", "Upper 95% CI", "R̂")

ft <- flextable(tab) %>%
  autofit()

doc <- read_docx() %>%
  body_add_flextable(ft)
print(doc, target = "Table/model_table_carnivory.docx")


# 3.3.2) Extract contrasts --------------------------------------------------


# extract the emmeans for the main effect model


emmeans(
  model_carnivory,
  ~ age,
  type = "response",
  offset = log(1)   # per unit effort
)


# age  rate lower.HPD upper.HPD
# A   0.095     0.064     0.131
# P   0.161     0.109     0.220
# 
# Results are averaged over the levels of: sex 
# Point estimate displayed: median 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95 

# predict carnivory rates for 25th and 75th percentil of competitive ability


q <- quantile(ILV_data$agon_PCA, c(0.25, 0.75))

emmeans(
  model_carnivory,
  ~ agon_PCA,
  at = list(
    agon_PCA = q,
    effort = 1
  ),
  type = "response"
)

# agon_PCA   rate lower.HPD upper.HPD
# -0.952 0.0852    0.0515     0.123
# 0.616 0.1572    0.1210     0.195
# 
# Results are averaged over the levels of: age, sex 
# Point estimate displayed: median 
# Results are back-transformed from the log scale 
# HPD interval probability: 0.95


# 4) NBDA - diffusion of carnivory on voles -------------------------------


# 4.1) Prepare data -------------------------------------------------------

# subset the ILV data frame
ILV_data_sub <- select(ILV_data, c("id", "age_sex", "hum_reactivity", "agon_PCA", "n_vole_burrows"))

# make sex a factor
ILV_data_sub$age_sex <- as.factor(ILV_data_sub$age_sex)

str(ILV_data_sub)
# make sure categorical variables are factors

# we don't need the spaital network, so we remove it from the edge-list
edge_list <- edge_list[,c(1:4)]

# 4.2) Create a data list ---------------------------------------------------

data_list <- import_user_STb(event_data = event_data, 
                             networks = edge_list ,
                             network_type = "undirected",
                             ILV_c = ILV_data_sub,
                             ILV_tv = ILV_tv,
                             ILVi = c("age_sex", "hum_reactivity", "agon_PCA", "n_vole_burrows", "presence"),
                             ILVs = c("age_sex", "n_vole_burrows", "presence"))

# ILVi - variables influencing asocial learning rate
# ILVs - variables influencing social learning rate

# generate model list
model_full_constant <- generate_STb_model(data_list, gq = T, est_acqTime = T, data_type = "discrete_time", intrinsic_rate = "constant")


# 4.3) Run the models -------------------------------------------------------

# social model with constant baseline
full_fit_constant <- fit_STb(data_list,
                    model_full_constant,
                    parallel_chains = 4,
                    chains = 4,
                    cores = 4,
                    iter = 2000,
                    refresh=1000
)


# social model with weibull baseline:

model_full_weibull <- generate_STb_model(data_list, gq = T, est_acqTime = T, data_type = "discrete_time", intrinsic_rate = "weibull")

full_fit_weibull <- fit_STb(data_list,
                             model_full_weibull,
                             parallel_chains = 4,
                             chains = 4,
                             cores = 4,
                             iter = 2000,
                             refresh=1000
)

# asocial model with constant baseline: 

model_full_constant_asocial <- generate_STb_model(data_list, gq = T, est_acqTime = T, data_type = "discrete_time", intrinsic_rate = "constant", model_type = "asocial")

full_fit_asocial_constant <- fit_STb(data_list,
                                    model_full_constant_asocial,
                                    parallel_chains = 4,
                                    chains = 4,
                                    cores = 4,
                                    iter = 2000,
                                    refresh=1000
)


# asocial model with weibull baseline: 
model_full_weibull_asocial <- generate_STb_model(data_list, gq = T, est_acqTime = T, data_type = "discrete_time", intrinsic_rate = "weibull", model_type = "asocial")

full_fit_asocial_weibull <- fit_STb(data_list,
                                    model_full_weibull_asocial,
                                    parallel_chains = 4,
                                    chains = 4,
                                    cores = 4,
                                    iter = 2000,
                                    refresh=1000
)



# 4.4) Model comparison -----------------------------------------------------


loo_output = STb_compare(full_fit_weibull, full_fit_constant, full_fit_asocial_weibull, full_fit_asocial_constant, method="loo-psis")

loo_output$comparison
#                           elpd_diff se_diff
# full_fit_weibull            0.0       0.0  
# full_fit_constant         -13.2       6.1  
# full_fit_asocial_constant -49.1       7.6  
# full_fit_asocial_weibull  -75.9      14.7

# weibull social model performs best, we therefore continue with that one

# look at the output
output_weibull <- STb_summary(full_fit_weibull, digits = 3)
output_weibull

# Parameter Median   MAD CI_Lower CI_Upper ess_bulk ess_tail  Rhat
# 1         log_lambda_0_mean -2.395 0.415   -3.221   -1.580 4644.477 3262.061 1.001
# 2          log_s_prime_mean -4.000 2.084   -8.329   -0.179 6709.280 2438.405 1.002
# 10                 lambda_0  0.091 0.038    0.027    0.182 4644.541 3262.061 1.001
# 12                        s  0.198 0.271    0.000    6.472 6525.793 2690.354 1.001
# 4  beta_ILVi_hum_reactivity -0.399 0.583   -1.521    0.808 6032.841 2395.588 1.001
# 5        beta_ILVi_agon_PCA  0.655 0.234    0.224    1.141 5445.967 3353.444 1.000
# 6  beta_ILVi_n_vole_burrows  0.068 0.075   -0.087    0.210 7113.676 2691.165 1.001
# 7        beta_ILVi_presence  0.048 0.033   -0.018    0.117 6545.883 3149.473 1.002
# 8  beta_ILVs_n_vole_burrows -0.015 1.011   -1.924    1.959 8447.321 2843.971 1.001
# 9        beta_ILVs_presence -0.015 1.030   -1.927    1.987 9874.911 2502.570 1.000
# 14     beta_ILVi_age_sex[1] -0.268 0.467   -1.133    0.666 5481.189 3327.219 1.000
# 17     beta_ILVs_age_sex[1]  0.006 0.969   -1.935    1.931 7569.997 2869.558 1.000
# 15     beta_ILVi_age_sex[2]  0.600 0.478   -0.353    1.533 5389.573 3010.188 1.000
# 18     beta_ILVs_age_sex[2]  0.021 1.012   -1.836    1.904 7484.423 2587.966 1.003
# 16     beta_ILVi_age_sex[3]  0.065 0.561   -1.100    1.159 5578.222 3149.619 1.000
# 19     beta_ILVs_age_sex[3]  0.013 1.038   -1.942    2.034 8800.070 2663.550 1.001
# 13            percent_ST[1]  0.000 0.000    0.000    0.000       NA       NA    NA
# 3                 log_gamma -0.757 0.182   -1.114   -0.406 6671.074 3367.527 1.001
# 11                    gamma  0.469 0.086    0.305    0.634 6671.106 3367.527 1.001

# 4.5) create a table -------------------------------------------------------


final_table <- output_weibull %>%
  
  # 1. Remove only log-scale parameters
  filter(!str_detect(Parameter, "^log_")) %>%
  
  # 2. Select relevant columns
  select(
    parameter = Parameter,
    median = Median,
    ci_lower = CI_Lower,
    ci_upper = CI_Upper,
    rhat = Rhat
  ) %>%
  
  # 3. Clean percent_ST naming (remove [1])
  mutate(
    parameter = str_replace(parameter, "percent_ST\\[1\\]", "percent_ST")
  ) %>%
  
  # 4. Add grouping column
  mutate(
    group = case_when(
      str_detect(parameter, "^beta_ILVi") ~ "ILV [asocial]",
      str_detect(parameter, "^beta_ILVs") ~ "ILV [social]",
      parameter %in% c("lambda_0", "s", "gamma", "percent_ST") ~ "Baseline",
      TRUE ~ "Other"
    )
  ) %>%
  
  # 5. Rename parameters (core ones + ILVs optional later)
  mutate(
    parameter = case_when(
      parameter == "lambda_0" ~ "Baseline asocial learning rate",
      parameter == "s" ~ "Social transmission rate",
      parameter == "gamma" ~ "Shape parameter",
      parameter == "percent_ST" ~ "Proportion social transmission",
      TRUE ~ parameter
    )
  ) %>%
  
  # 6. Order groups
  mutate(
    group = factor(group, levels = c("Baseline", "ILV [asocial]", "ILV [social]", "Other"))
  ) %>%
  
  arrange(group, parameter)

final_table <- final_table %>%
  
  # 1. Create ordering variable (ONLY for sorting)
  mutate(
    param_order = case_when(
      parameter == "Baseline asocial learning rate" ~ 1,
      parameter == "Social transmission rate" ~ 2,
      parameter == "Proportion social transmission" ~ 3,
      parameter == "Shape parameter" ~ 4,
      TRUE ~ 5
    )
  ) %>%
  
  # 2. Arrange using that
  arrange(group, param_order, parameter) %>%
  
  # 3. Drop helper column
  select(-param_order) %>%
  
  # 4. Move group to front
  select(group, everything())


final_table <- final_table %>%
  mutate(
    parameter = dplyr::recode(parameter,
                       "beta_ILVi_hum_reactivity" = "Reactivity to humans",
                       "beta_ILVi_agon_PCA" = "Agonistic tendencies",
                       "beta_ILVi_n_vole_burrows" = "Vole burrow prevalence",
                       "beta_ILVi_presence" = "Presence in study area",
                       
                       "beta_ILVs_n_vole_burrows" = "Vole burrow prevalence",
                       "beta_ILVs_presence" = "Presence in study area",
                       "beta_ILVi_age_sex[1]" = "Age/sex [adult male]",
                       "beta_ILVi_age_sex[2]" = "Age/sex [juvenile female]",
                       "beta_ILVi_age_sex[3]" = "Age/sex [juvenile male]",
                       "beta_ILVs_age_sex[1]" = "Age/sex [adult male]",
                       "beta_ILVs_age_sex[2]" = "Age/sex [juvenile female]",
                       "beta_ILVs_age_sex[3]" = "Age/sex [juvenile male]",
    )
  )

colnames(final_table) <- c("Group", "Parameter", "Median", "Lower CI", "Upper CI", "Rhat")


# convert to flextable
ft <- flextable(final_table)

# basic formatting
ft <- ft %>%
  autofit() %>%
  theme_booktabs() %>%
  bold(part = "header") %>%
  align(align = "center", part = "all") %>%
  align(j = "Parameter", align = "left", part = "all") %>%
  align(j = "Group", align = "left", part = "all")



# create Word document
doc <- read_docx() %>%
  body_add_par("NBDA Results", style = "heading 1") %>%
  body_add_flextable(ft)

# save
print(doc, target = "Table/NBDA_results.docx")


# 5) Figures --------------------------------------------------------------


# 5.1) Social network --------------------------------------------------------------------



# remove self associations
edge_list_plot <- edge_list[edge_list$assoc >= 0.03 & edge_list$focal != edge_list$other, ]

# list of those with vole carnivory evidence
learners <- event_data[event_data$time<39,"id"]


g <- graph_from_data_frame(
  d = edge_list_plot,
  directed = FALSE
)

e <- as_edgelist(g,names=FALSE)


# Fruchtermann Reigold algorithm
lay <- qgraph.layout.fruchtermanreingold(
  e,
  vcount = vcount(g)
)


lay_scaled <- norm_coords(lay, xmin = -1, xmax = 1, ymin = -1, ymax = 1)

# node color
node_colors <- ifelse(V(g)$name %in% learners, "#e1a136", "grey80")

agon_vals <- ILV_data$agon_PCA
names(agon_vals) <- ILV_data$id

node_size <- agon_vals[V(g)$name]

node_size_scaled <- scales::rescale(node_size, to = c(4, 8))

edge_width <- scales::rescale(E(g)$assoc, to = c(1, 5))


# ---- Convert igraph/base plot into a ggplot object ----
igraph_plot <- as.ggplot(function() {
  par(mar = c(0, 0, 0, 0))   # removes base R outer margins
  plot(
    g,
    layout = lay_scaled,
    vertex.size = node_size_scaled,
    vertex.label = NA,
    edge.width = edge_width,
    vertex.color = node_colors,
    edge.color = "grey40",
    margin = 0               # removes igraph internal margin
  )
})



# 5.2) Cumulative diffusion curve -------------------------------------------------------------------

# we need to store num inds per trial to refer to later
event_data_plot <- event_data %>%
  group_by(trial) %>%
  mutate(n_trial = n())

plot_data_obs <- event_data_plot %>%
  filter(time > 0, time <= t_end) %>% # exclude demonstrators (time == 0) and censored (time > t_end)
  group_by(trial) %>%
  arrange(time, .by_group = TRUE) %>%
  mutate(
    cum_prop = row_number() / n_trial, # this denominator needs to be the number of individuals per trial
    type = "observed"
  ) %>%
  select(trial, time, cum_prop, type) %>%
  ungroup()

# add in 0,0 starting point
plot_data_obs <- bind_rows(
  plot_data_obs,
  plot_data_obs %>%
    distinct(trial) %>%
    mutate(time = 0, cum_prop = 0, type = "observed")
) %>%
  arrange(trial, time)


draws_df <- as_draws_df(full_fit_weibull$draws(variables = "acquisition_time", inc_warmup = FALSE))

# pivot longer
ppc_long <- draws_df %>%
  select(starts_with("acquisition_time[")) %>%
  pivot_longer(
    cols = everything(),
    names_to = c("trial", "ind"),
    names_pattern = "acquisition_time\\[(\\d+),(\\d+)\\]",
    values_to = "time"
  ) %>%
  mutate(
    trial = as.integer(trial),
    ind = as.integer(ind),
    draw = rep(1:(nrow(draws_df)), 
               each = length(unique(.$trial)) * length(unique(.$ind)))
  )

# thin sample for plotting
sample_idx <- sample(c(1:max(ppc_long$draw)), 100)
ppc_long <- ppc_long %>% filter(draw %in% sample_idx)

# same as before, we need a way to reference the number of individuals in each trial
ppc_long <- ppc_long %>%
  group_by(draw, trial) %>%
  mutate(n_trial = n())

# we also need to remove individuals predicted as censored, which
# have value of -1 in predicted data
ppc_long <- ppc_long %>%
  filter(time > -1)

# build cumulative curves per draw
plot_data_ppc <- ppc_long %>%
  group_by(draw, trial, time) %>%
  summarise(n = n(), n_trial = first(n_trial), .groups = "drop") %>%
  group_by(draw, trial) %>%
  arrange(time) %>%
  mutate(cum_prop = cumsum(n) / n_trial)

# add in 0,0 starting point
plot_data_ppc <- bind_rows(
  plot_data_ppc,
  plot_data_ppc %>%
    distinct(trial, draw) %>%
    mutate(time = 0, cum_prop = 0, type = "ppc")
) %>%
  arrange(trial, time)

tiff("Figures/Cumulative learning curve_model check.tiff",
     width = 5,
     height = 5,
     units = "in",
     res = 600,
     compression = "lzw")


# plot it
ggplot() +
  geom_line(data = plot_data_ppc, 
            aes(x = time, y = cum_prop, 
                group = interaction(draw, trial)), alpha = .1) +
  geom_line(data = plot_data_obs, aes(x = time, y = cum_prop), linewidth = 1) +
  labs(x = "Time [days]", y = "Cumulative proportion\nof informed individuals", color = "Trial") +
  theme_minimal()+
  ylim(c(0,1))


dev.off()

# only cumulative prop of informed individuals


diff_curve <- ggplot() +
  # geom_line(data = plot_data_ppc, 
  #           aes(x = time, y = cum_prop, 
  #               group = interaction(draw, trial)), alpha = .1) +
  geom_line(data = plot_data_obs, aes(x = time, y = cum_prop), linewidth = 1) +
  labs(x = "Time [days]", y = "Cumulative proportion\nof informed individuals", color = "Trial") +
  theme_minimal()+
  ylim(c(0,1))


tiff("Figures/soc_net_Cumulative learning curve.tiff",
 width = 8,
 height = 4,
 units = "in",
 res = 600,
 compression = "lzw")



# ---- Combine with our ggplot ----
combined_plot <- igraph_plot + diff_curve +
  plot_annotation(tag_levels = "a")


# Show figure
combined_plot



dev.off()



# 5.3) Aggression/age overall carnivory -----------------------------------


library(ggpubr)
library(ggplot2)

# Data from emmeans
age_dat <- data.frame(
  age = c("A", "P"),
  rate = c(0.095, 0.161),
  lower = c(0.064, 0.109),
  upper = c(0.131, 0.220)
)

# Age panel
p_age <- ggplot(age_dat, aes(x = age, y = rate)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1) +
  theme_bw() +
  xlab("Age category") +
  ylab("") +
  scale_x_discrete(labels = c("A" = "Adult", "P" = "Juvenile"))


p.carn <- plot(conditional_effects(model_carnivory))

tiff("Figures/Figure 3.tiff",
     width = 8,
     height = 4,
     units = "in",
     res = 600,
     compression = "lzw")

ggarrange(
  # emmeans-based age panel
  p.carn$agon_PCA +
    theme_bw() +
    xlab("Competitive ability [PCA]") +
    ylab("Carnivory rate"),
  p_age,
  labels = c("c", "d")
)

dev.off()






