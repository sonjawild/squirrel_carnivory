# Competitive ability predicts the spread of novel carnivory in California ground squirrels 


# 1) Load libraries -------------------------------------------------------

# install devtools if not already
if (!require("devtools")) install.packages("devtools")
#devtools::install_github("michaelchimento/STbayes")
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
edge_list <- read.csv("Input Data/edge_list.csv", row.names = 1)
ILV_data <- read.csv("Input Data/ILVs.csv", row.names = 1)
event_data_row <- read.csv("Input Data/rowwise_event.csv", row.names = 1)


# 3) Modeling mammal carnivory (binary) -----------------------------------


# is the rate of carnivory predicted by individual covariates? we combine both molecular and obsrvational data for this analysis.


# 3.1) descriptive stats --------------------------------------------------

# some descriptive stats on age-sex classes that consumed mammals and/or voles, respectively:
# we use the ILV_data frame since that contains a summary of the by event data frame


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


# 3.2) Calculate VIFs -----------------------------------------------------

event_data_row %>% count(opportunity_type, carnivory_pos)

# opportunity_type carnivory_pos   n
# 1              obs             0 430
# 2              obs             1  47
# 3   trap_no_sample             0 175
# 4   trap_no_sample             1   6
# 5 trap_with_sample             0  54
# 6 trap_with_sample             1  24

# since the different opportunities for vole detection have diffrent detection probabilities, we include them in our model as a covariate. 


model_vif <- glm(
  formula = carnivory_pos ~ age + sex + hum_reactivity + agon_PCA + n_vole_burrows + 
    opportunity_type,
  data = event_data_row,
  family = binomial(link = "logit")
)

car::vif(model_vif)
# GVIF Df GVIF^(1/(2*Df))
# age              1.403818  1        1.184828
# sex              1.233214  1        1.110502
# hum_reactivity   1.057848  1        1.028517
# agon_PCA         1.555810  1        1.247321
# n_vole_burrows   1.099067  1        1.048364
# opportunity_type 1.011510  2        1.002865

# we run a main effects model only 

model_carnivory <- brm(
  carnivory_pos ~ age + sex + hum_reactivity + agon_PCA + n_vole_burrows + 
    opportunity_type + (1 | id),
  data = event_data_row,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  iter = 2000,
  seed = 123
)
#save(model_carnivory, file="output/model_carnivory.RDA")
load("output/model_carnivory.RDA")

# and a interactive effects model with age

model_carnivory_int <- brm(
  carnivory_pos ~ age*sex + hum_reactivity*age + agon_PCA*age + n_vole_burrows + opportunity_type  + (1 | id),
  data = event_data_row,
  family = bernoulli(link = "logit"),
  chains = 4,
  cores = 4,
  iter = 2000,
  seed = 123
)

#save(model_carnivory_int, file="output/model_carnivory_int.RDA")
load("output/model_carnivory_int.RDA")

# we compare the two 
loo_compare(loo(model_carnivory), loo(model_carnivory_int))
# elpd_diff se_diff
# model_carnivory      0.0       0.0   
# model_carnivory_int -0.6       1.1 

# adding the age interactions does not improve the model, so we are going with the simpler main effects model

# we create tables for both

# first interaction model
tab_int <- as.data.frame(summary(model_carnivory_int)$fixed)
tab_int$Parameter <- rownames(tab_int)

tab_int <- tab_int %>%
  rename(
    Estimate = Estimate,
    CI_low = 'l-95% CI',
    CI_high = 'u-95% CI',
    Rhat = Rhat
  ) %>%
  mutate(Parameter = case_match(Parameter,
                                "Intercept" ~ "Baseline carnivory (adult reference)",
                                "agon_PCA" ~ "Agonistic tendency (PCA), adults",
                                "hum_reactivity" ~ "Reactivity to humans, adults",
                                "n_vole_burrows" ~ "Vole burrow prevalence",
                                "ageP" ~ "Age (juvenile vs adult)",
                                "sexM" ~ "Sex (male vs female), adults",
                                "opportunity_typetrap_no_sample" ~ "Evidence (trap no sample vs observation)",
                                "opportunity_typetrap_with_sample" ~ "Evidence (trap with sample vs observation)",
                                "ageP:sexM" ~ "Sex × Age",
                                "ageP:hum_reactivity" ~ "Reactivity to humans × Age",
                                "ageP:agon_PCA" ~ "Agonistic tendency × Age"
  )) %>%
  mutate(
    Estimate = exp(Estimate),
    CI_low = exp(CI_low),
    CI_high = exp(CI_high)
  ) %>%
  mutate(across(c(Estimate, CI_low, CI_high, Rhat), ~round(.x, 3))) %>%
  select(Parameter, Estimate, CI_low, CI_high, Rhat)


colnames(tab_int) <- c("Predictor", "Odds ratio", "Lower 95% CI", "Upper 95% CI", "R̂")

rownames(tab_int) <- NULL

ft <- flextable(tab_int) %>%
  autofit()

doc <- read_docx() %>%
  body_add_flextable(ft)
print(doc, target = "output/model_table_carnivory_w_interactions.docx")

## main effect model
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
                                "Intercept" ~ "Baseline carnivory (adult reference)",
                                "agon_PCA" ~ "Agonistic tendency (PCA), adults",
                                "hum_reactivity" ~ "Reactivity to humans, adults",
                                "n_vole_burrows" ~ "Vole burrow prevalence",
                                "ageP" ~ "Age (juvenile vs adult)",
                                "sexM" ~ "Sex (male vs female), adults",
                                "opportunity_typetrap_no_sample" ~ "Evidence (trap no sample vs observation)",
                                "opportunity_typetrap_with_sample" ~ "Evidence (trap with sample vs observation)"
  )) %>%
  mutate(
    Estimate = exp(Estimate),
    CI_low = exp(CI_low),
    CI_high = exp(CI_high)
  ) %>%
  mutate(across(c(Estimate, CI_low, CI_high, Rhat), ~round(.x, 3))) %>%
  select(Parameter, Estimate, CI_low, CI_high, Rhat)


colnames(tab) <- c("Predictor", "Odds ratio", "Lower 95% CI", "Upper 95% CI", "R̂")

rownames(tab) <- NULL

ft <- flextable(tab) %>%
  autofit()

doc <- read_docx() %>%
  body_add_flextable(ft)
print(doc, target = "output/model_table_carnivory.docx")


# 3.3) Extract contrasts --------------------------------------------------

# get contrasts for opportunity types 

emm_type <- emmeans(
  model_carnivory,
  ~ opportunity_type,
  type = "response"
)
emm_type
# opportunity_type response lower.HPD upper.HPD
# obs                0.0662   0.03236    0.1005
# trap_no_sample     0.0189   0.00379    0.0399
# trap_with_sample   0.2615   0.14191    0.3970
# 
# Results are averaged over the levels of: age, sex 
# Point estimate displayed: median 
# Results are back-transformed from the logit scale 
# HPD interval probability: 0.95 

contrast(emm_type, "pairwise")
# contrast                          odds.ratio lower.HPD upper.HPD
# obs / trap_no_sample                  3.6526    1.2692     8.560
# obs / trap_with_sample                0.1988    0.0953     0.362
# trap_no_sample / trap_with_sample     0.0541    0.0125     0.124
# 
# Results are averaged over the levels of: age, sex 
# Point estimate displayed: median 
# Results are back-transformed from the log odds ratio scale 
# HPD interval probability: 0.95 

# predict carnivory rates for 25th and 75th percentil of competitive ability


q <- quantile(ILV_data$agon_PCA, c(0.25, 0.75))

emm <- emmeans(
  model_carnivory,
  ~ agon_PCA,
  at = list(agon_PCA = q, opportunity_type = "obs"),
  type = "response"
)
emm

# agon_PCA response lower.HPD upper.HPD
# -0.952   0.0305   0.00936    0.0592
# 0.616   0.0923   0.04832    0.1429
# 
# Results are averaged over the levels of: age, sex 
# Point estimate displayed: median 
# Results are back-transformed from the logit scale 
# HPD interval probability: 0.95 

contrast(emm, "revpairwise", type = "response")

# contrast                                                 odds.ratio lower.HPD upper.HPD
# agon_PCA0.615777272631032 / (agon_PCA-0.951729812201686)       3.22      1.11      6.63
# 
# Results are averaged over the levels of: age, sex 
# Point estimate displayed: median 
# Results are back-transformed from the log odds ratio scale 
# HPD interval probability: 0.95 

# Individuals at the 75th percentile of agonistic tendency have 3.22 x higher odds of detected carnivory than those at the 25th percentile



# 3.4) Model performance and checks ---------------------------------------

pp_check(model_carnivory)
# looks like a great fit

plot(model_carnivory)
# check for stationarity and mixing (also looks good)

r2_bayes(model_carnivory)
# # Bayesian R2 with Compatibility Interval
# 
# Conditional R2: 0.178 (95% CI [0.111, 0.246])
# Marginal R2: 0.100 (95% CI [0.037, 0.176])

# look at the output visually
plot(conditional_effects(model_carnivory))


# 4) NBDA - diffusion of carnivory on voles -------------------------------


# 4.1) Prepare data -------------------------------------------------------

# subset the ILV data frame
ILV_data_sub <- select(ILV_data, c("id", "age_sex", "hum_reactivity", "agon_PCA", "n_vole_burrows"))

# make sex a factor
ILV_data_sub$age_sex <- as.factor(ILV_data_sub$age_sex)

str(ILV_data_sub)
# make sure categorical variables are factors

# from our row-wise event data, we extract the first events of vole carnivory

event_data_row <- event_data_row %>%
  mutate(date = as.Date(date)) %>%
  arrange(id, date) %>%
  group_by(id) %>%
  mutate(
    first_vole_evidence = vole_pos == 1 &
      date == min(date[vole_pos == 1], na.rm = TRUE)
  ) %>%
  ungroup()

# Acquisition date per learner ---
acquisition_dates <- event_data_row %>%
  filter(first_vole_evidence) %>%
  distinct(id, date)

nrow(acquisition_dates)   # expect 29

# time scale = calendar days since the start of the study
study_start <- min(event_data_row$date, na.rm = TRUE)
study_end   <- max(event_data_row$date, na.rm = TRUE)

event_time_map <- tibble(date = seq(study_start, study_end, by = "day")) %>%
  mutate(time = as.integer(date - study_start) + 1)

nrow(event_time_map)      # 62
t_end_int <- max(event_time_map$time)
t_end_int                 # expect 62

# Acquisition dates per learner, looked up against the full time scale ---
acquisition_dates <- event_data_row %>%
  filter(first_vole_evidence) %>%
  distinct(id, date) %>%
  left_join(event_time_map, by = "date")

nrow(acquisition_dates)       # expect 29

# Full event_data: all IDs, censored individuals filled in ---
all_ids <- event_data_row %>% distinct(id)

event_data_new <- all_ids %>%
  left_join(acquisition_dates, by = "id") %>%
  mutate(
    time  = if_else(is.na(time), t_end_int + 1, time),
    date  = as.character(date),
    t_end = t_end_int,
    trial = 1
  ) %>%
  select(id, time, date, t_end, trial)

str(event_data_new)                          # 59 obs; time now spans up to t_end_int+1
sum(event_data_new$time == t_end_int + 1)    # expect 30
sum(event_data_new$time <= t_end_int)        # expect 29




# --- discretization grid: only the distinct time values that actually appear in event_data ---
event_times <- sort(unique(event_data_new$time))
length(event_times)   # expect ~17

# map each of these time labels back to its cutoff calendar date
time_map <- tibble(time = event_times) %>%
  mutate(
    t_end = if_else(
      time == t_end_int + 1,
      study_end,                        # censoring step -> cutoff = last observed day
      study_start + (time - 1)          # real elapsed-day time labels -> corresponding calendar date
    )
  )

nrow(time_map)                                 # expect 17
anyNA(time_map$t_end)                          # FALSE
n_distinct(time_map$time) == nrow(time_map)    # TRUE

# --- id x time scaffold, now only at real event-boundary times ---
ids <- event_data_row %>% distinct(id)

ILV_tv_new <- time_map %>%
  tidyr::crossing(ids) %>%
  arrange(time, id)

nrow(ILV_tv_new)   # expect length(event_times) * 59, ~1003

# --- cumulative counts per opportunity type ---
count_type <- function(id_val, t_end_val, type_val) {
  sum(event_data_row$id == id_val &
        event_data_row$date <= t_end_val &
        event_data_row$opportunity_type == type_val)
}

ILV_tv_new <- ILV_tv_new %>%
  rowwise() %>%
  mutate(
    obs_days_cum            = count_type(id, t_end, "obs"),
    trap_no_sample_days_cum = count_type(id, t_end, "trap_no_sample"),
    trap_sample_days_cum    = count_type(id, t_end, "trap_with_sample")
  ) %>%
  ungroup()

# --- weighted composite ---
w_obs              <- 1
w_trap_no_sample   <- exp(fixef(model_carnivory)["opportunity_typetrap_no_sample", "Estimate"])
w_trap_with_sample <- exp(fixef(model_carnivory)["opportunity_typetrap_with_sample", "Estimate"])

ILV_tv_new <- ILV_tv_new %>%
  mutate(effective_effort = w_obs * obs_days_cum +
           w_trap_no_sample * trap_no_sample_days_cum +
           w_trap_with_sample * trap_sample_days_cum)

ILV_tv_new$effective_effort <- scale(ILV_tv_new$effective_effort)

ILV_tv_weighted <- ILV_tv_new %>%
  transmute(trial = 1, id = id, time = time, effective_effort = effective_effort)

str(ILV_tv_weighted)   # expect 1003 obs of 4 variables
hist(ILV_tv_weighted$effective_effort)

setdiff(unique(event_data_new$time), unique(ILV_tv_weighted$time))   # should be integer(0)
setdiff(event_data_new$id, ILV_tv_weighted$id)                        # should be character(0)


# 4.2) Create a data list ---------------------------------------------------

data_list <- import_user_STb(event_data = event_data_new, 
                             networks = edge_list ,
                             network_type = "undirected",
                             ILV_c = ILV_data_sub,
                             ILV_tv = ILV_tv_weighted,
                             ILVi = c("age_sex", "hum_reactivity", "agon_PCA", "n_vole_burrows", "effective_effort"),
                             ILVs = c("age_sex", "n_vole_burrows", "effective_effort"))

# ILVi - variables influencing asocial learning rate
# ILVs - variables influencing social learning rate

# the warning about the asymmetry of the network can be safely ignored. The supplied network only contains the upper half of the matrix and we let ST Bayes fill out the lower half (symmetric network)




# 4.3) Run the models -------------------------------------------------------

# generate model list
model_full_constant <- generate_STb_model(data_list, gq = T, est_acqTime = T, data_type = "discrete_time", intrinsic_rate = "constant")


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
# model elpd_diff se_diff p_worse diag_diff      diag_elpd
# full_fit_weibull       0.0     0.0      NA                         
# full_fit_constant      -0.3     0.8    0.63   N < 100 2 k_psis > 0.7
# full_fit_asocial_constant     -29.8     6.3    1.00   N < 100 1 k_psis > 0.7
# full_fit_asocial_weibull     -30.3     6.4    1.00   N < 100 4 k_psis > 0.7
# 
# Diagnostic flags present.
# See ?`loo-glossary` (sections `diag_diff` and `diag_elpd`)
# or https://mc-stan.org/loo/reference/loo-glossary.html.


# Weibull and constant are pretty much the same - so we go with the simpler model (cosntant) for reporting

# when looking at the output of the weibull, we can also see that the confidence itnervals for the shape parameter include 0 (log_gamma), which further suggests that it is no different to a constant bseline

output_constant <- STb_summary(full_fit_constant, digits = 3)
output_constant
# Parameter Median   MAD CI_Lower CI_Upper ess_bulk ess_tail  Rhat
# 1           log_lambda_0_mean -5.139 0.470   -6.180   -4.341 1947.763 2434.302 1.000
# 2            log_s_prime_mean -3.717 0.912   -5.845   -2.041 1672.833 1795.908 1.001
# 9                    lambda_0  0.006 0.003    0.002    0.012 1947.768 2434.302 1.000
# 10                          s  4.156 4.334    0.005   24.790 1420.255 1722.648 1.001
# 3    beta_ILVi_hum_reactivity -0.769 0.709   -2.094    0.627 4478.786 3353.912 1.000
# 4          beta_ILVi_agon_PCA  0.898 0.283    0.366    1.463 3353.685 2812.636 1.001
# 5    beta_ILVi_n_vole_burrows -0.147 0.147   -0.473    0.104 2236.896 2066.314 1.001
# 6  beta_ILVi_effective_effort  0.565 0.237    0.053    1.029 3599.745 2144.638 1.000
# 7    beta_ILVs_n_vole_burrows  0.424 0.180    0.051    0.833 2474.720 2237.985 1.001
# 8  beta_ILVs_effective_effort  0.075 0.467   -0.923    1.010 3918.830 2759.244 1.001
# 12       beta_ILVi_age_sex[1]  0.040 0.609   -1.254    1.165 3052.269 2968.619 1.000
# 15       beta_ILVs_age_sex[1] -0.436 0.681   -1.864    0.941 4203.330 2501.755 1.001
# 13       beta_ILVi_age_sex[2]  0.333 0.598   -0.938    1.487 3574.763 2739.819 1.001
# 16       beta_ILVs_age_sex[2]  0.424 0.758   -1.087    1.961 3279.388 2990.692 1.000
# 14       beta_ILVi_age_sex[3] -0.039 0.710   -1.463    1.227 4462.475 2648.576 1.001
# 17       beta_ILVs_age_sex[3] -0.273 0.797   -1.877    1.228 4329.560 3076.936 1.000
# 11              percent_ST[1]  0.344 0.086    0.164    0.521 1412.652 1484.787 1.001



# 4.5) create a table -------------------------------------------------------


final_table <- output_constant %>%
  
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
                       "beta_ILVi_effective_effort" = "Sampling opportunity",
                       
                       "beta_ILVs_n_vole_burrows" = "Vole burrow prevalence",
                       "beta_ILVs_effective_effort" = "Sampling opportunity",                       "beta_ILVi_age_sex[1]" = "Age/sex [adult male]",
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
print(doc, target = "output/NBDA_results.docx")


# 5) Figures --------------------------------------------------------------


# 5.1) Social network --------------------------------------------------------------------



# remove self associations
edge_list_plot <- edge_list[edge_list$assoc >= 0.03 & edge_list$focal != edge_list$other, ]

# list of those with vole carnivory evidence
learners <- event_data_new$id[event_data_new$time<57]


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


library(ggplot2)
library(cowplot)

# Create a custom legend
legend_plot <- ggplot(data.frame(
  x = c(1.0, 2.2),
  y = c(1, 1),
  label = c("vole-eaters", "non-vole eaters"),
  colour = c("orange", "grey70")
)) +
  geom_point(aes(x = x, y = y, colour = colour), size = 4) +
  geom_text(aes(x = x + 0.12, y = y, label = label),
            hjust = 0, size = 4) +
  scale_colour_identity() +
  coord_cartesian(xlim = c(0.8, 3.5), ylim = c(0.8, 1.2),
                  clip = "off") +
  theme_void()

# Combine network and legend
soc.net <- plot_grid(
  igraph_plot,
  legend_plot,
  ncol = 1,
  rel_heights = c(1, 0.02)
)

# 5.2) Cumulative diffusion curve -------------------------------------------------------------------

# we need to store num inds per trial to refer to later
event_data_plot <- event_data_new %>%
  group_by(trial) %>%
  mutate(n_trial = n())

# Keep individual acquisition events (no summarise!)
plot_data_obs <- event_data_plot %>%
  filter(time <= t_end) %>%
  arrange(trial, time, id) %>%
  group_by(trial) %>%
  mutate(
    cum_prop = row_number() / first(n_trial),
    type = "observed"
  ) %>%
  ungroup() %>%
  select(trial, time, cum_prop, type)

# Add starting and closing points
plot_data_obs <- bind_rows(
  plot_data_obs,
  
  plot_data_obs %>%
    distinct(trial) %>%
    mutate(time = 0, cum_prop = 0, type = "observed"),
  
  event_data_plot %>%
    distinct(trial, t_end) %>%
    left_join(
      plot_data_obs %>%
        group_by(trial) %>%
        summarise(final_cum_prop = max(cum_prop), .groups = "drop"),
      by = "trial"
    ) %>%
    transmute(
      trial,
      time = t_end,
      cum_prop = final_cum_prop,
      type = "observed"
    )
) %>%
  arrange(trial, time, cum_prop)












plot_data_obs <- bind_rows(
  plot_data_obs,
  # 0,0 starting point (your existing code)
  plot_data_obs %>%
    distinct(trial) %>%
    mutate(time = 0, cum_prop = 0, type = "observed"),
  # NEW: flat-line closing point at t_end
  event_data_plot %>%
    distinct(trial, t_end) %>%
    left_join(
      plot_data_obs %>%
        group_by(trial) %>%
        summarise(final_cum_prop = max(cum_prop), .groups = "drop"),
      by = "trial"
    ) %>%
    transmute(trial, time = t_end, cum_prop = final_cum_prop, type = "observed")
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
  labs(x = "Time steps [days]", y = "Cumulative proportion\nof informed individuals", color = "Trial") +
  theme_minimal()+
  ylim(c(0,1))


dev.off()

# only cumulative prop of informed individuals


diff_curve <- ggplot() +
  # geom_line(data = plot_data_ppc, 
  #           aes(x = time, y = cum_prop, 
  #               group = interaction(draw, trial)), alpha = .1) +
  geom_line(data = plot_data_obs, aes(x = time, y = cum_prop), linewidth = 1) +
  labs(x = "Time steps [days]", y = "Cumulative proportion\nof informed individuals", color = "Trial") +
  theme_minimal()+
  ylim(c(0,1))




# 5.3) Aggression/age overall carnivory -----------------------------------


library(ggpubr)
library(ggplot2)


p.age <- plot(conditional_effects(model_carnivory))


tiff("Figures/Figure 3.tiff",
     width = 8,
     height = 4,
     units = "in",
     res = 600,
     compression = "lzw")

ggarrange(
  p.age$age+
    theme_bw()+
    xlab("Age category")+
    ylab("Carnivory probability") +
    scale_x_discrete(labels = c("A" = "Adult", "P" = "Juvenile"))+
    ylim(c(0,0.5))
   # theme(
      # Makes the text invisible but keeps its bounding box space
   #   axis.text.y = element_text(color = "transparent") 
   # )
  ,
  # emmeans-based age panel
  p.age$agon_PCA +
    theme_bw() +
    xlab("Competitive ability [PC score]") +
    ylab("")+
    ylim(c(0,0.5))+
    theme(
  #  Makes the text invisible but keeps its bounding box space
    axis.text.y = element_text(color = "transparent") 
    )
 ,
  labels = c("c", "d")
)

dev.off()



# 5.4) Individual trajectories --------------------------------------------



evidence_df <- event_data_row %>%
  mutate(
    date = as.Date(date),
    evidence = case_when(
      opportunity_type == "obs" & vole_pos == 1                     ~ "obs_positive",
      opportunity_type == "obs" & vole_pos == 0                     ~ "obs_negative",
      opportunity_type == "trap_no_sample" & vole_pos == 1          ~ "trap_positive",
      opportunity_type == "trap_no_sample" & vole_pos == 0          ~ "trap_negative",
      opportunity_type == "trap_with_sample" & vole_pos == 1        ~ "fecal_positive",
      opportunity_type == "trap_with_sample" & vole_pos == 0        ~ "fecal_negative",
      TRUE ~ NA_character_
    )
  )

id_order <- evidence_df %>%
  group_by(id) %>%
  summarise(first_evidence_date = min(date[first_vole_evidence], na.rm = TRUE)) %>%
  mutate(first_evidence_date = if_else(is.infinite(first_evidence_date), as.Date(NA), first_evidence_date)) %>%
  arrange(first_evidence_date) %>%   # earliest learners first
  pull(id)

# numeric label 1:59 following that same order (1 = earliest learner)
id_lookup <- tibble(id = id_order, id_num = seq_along(id_order))

evidence_df <- evidence_df %>%
  left_join(id_lookup, by = "id") %>%
  mutate(id_num = factor(id_num, levels = 1:length(id_order)))  # NOT reversed -> first learner (1) sits at bottom

ind.plot <- ggplot(evidence_df, aes(x = date, y = id_num,
                        shape = evidence,
                        colour = first_vole_evidence)) +
  geom_point(size = 1.8) +
  scale_shape_manual(
    values = c(
      obs_negative   = 1,
      obs_positive   = 16,
      trap_negative  = 2,
      trap_positive  = 17,
      fecal_negative = 0,
      fecal_positive = 15
    ),
    labels = c(
      obs_negative   = "Observation (-)",
      obs_positive   = "Observation (+)",
      trap_negative  = "Trap day, no sample (-)",
      trap_positive  = "Trap day, no sample (+)",
      fecal_negative = "Fecal sample (-)",
      fecal_positive = "Fecal sample (+)"
    ),
    name = NULL
  ) +
  scale_colour_manual(
    values = c(
      `FALSE` = "black",
      `TRUE`  = "#e1a136"
    ),
    guide = "none"
  ) +
  labs(
    x = "Date",
    y = "Individual"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 6)
  )


# Arrange the three plots


tiff("Figures/soc_net_Cumulative learning curve.tiff",
     width = 7,
     height = 7.5,
     units = "in",
     res = 600,
     compression = "lzw")



# ---- Combine with our ggplot ----
library(patchwork)
layout <- "
AB
CC
"

ind.plot <- ind.plot + theme(legend.position = "bottom")

combined_plot <- soc.net + diff_curve + ind.plot +
  plot_layout(design = layout,
              widths  = c(1.1, 1),   # A wider than B
              heights = c(1, 1.7)) + # row 2 (C) taller than row 1
  plot_annotation(tag_levels = "a")


# Show figure
combined_plot



dev.off()


# how many opportunities were needed to become vole positive
opportunity_to_first_evidence <- event_data_row %>%
  mutate(date = as.Date(date)) %>%
  arrange(id, date) %>%
  group_by(id) %>%
  mutate(opportunity_number = row_number()) %>%   # rank of this opportunity among all of this id's opportunities
  filter(first_vole_evidence) %>%
  summarise(opportunity_number = opportunity_number, .groups = "drop")

# how many individuals were positive on their very first opportunity
n_first_opportunity <- opportunity_to_first_evidence %>%
  filter(opportunity_number == 1) %>%
  nrow()

n_first_opportunity
nrow(opportunity_to_first_evidence)   # total learners, should be 29, for context

# average number of opportunities needed (mean, SD, range) across all learners
summary_stats <- opportunity_to_first_evidence %>%
  summarise(
    mean_opportunities = mean(opportunity_number),
    median_opportunities = median(opportunity_number),
    
        sd_opportunities    = sd(opportunity_number),
    min_opportunities   = min(opportunity_number),
    max_opportunities   = max(opportunity_number)
  )

summary_stats

# # A tibble: 1 × 5
# mean_opportunities median_opportunities sd_opportunities min_opportunities max_opportunities
# <dbl>                <int>            <dbl>             <int>             <int>
#   1               7.62                    7             4.73                 2                17
