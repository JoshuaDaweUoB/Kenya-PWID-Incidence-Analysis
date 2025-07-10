# load packages
pacman::p_load(dplyr, tidyr, readr, readxl, lubridate, sandwich, lmtest)

# set working directory
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Publications/Sex work and risk of HIV and HCV/Emails to authors/Kenya data/data")

# load hiv data
kenya_clean <- read_csv("kenya_clean.csv")

# convert character variables to factors
kenya_clean <- kenya_clean %>%
  mutate(
    incarc_life = as.factor(incarc_life),
    OSTACCESSED = as.factor(OSTACCESSED),
    Hbinary = as.factor(Hbinary),
    sw_ever_partner = as.factor(sw_ever_partner),
    sw_ever = as.factor(sw_ever),
    sw_ever_gender = as.factor(sw_ever_gender),
    sw_30d = as.factor(sw_30d),
    msm_30d = as.factor(msm_30d)
)

# stratify by sex
kenya_males <- kenya_clean %>% filter(GENDER == "Male")
kenya_females <- kenya_clean %>% filter(GENDER == "Female")

# filter to sequence = 1 and create summary
baseline_data <- kenya_clean %>% filter(id_seq == 1)

# stratify by sex
baseline_data_males <- kenya_males %>% filter(id_seq == 1)
baseline_data_females <- kenya_females %>% filter(id_seq == 1)

# Function to create summary table for baseline data
create_baseline_summary <- function(data, dataset_name = "") {
  baseline_data <- data %>% filter(id_seq == 1)
  
  cat("\n=== Baseline Summary:", dataset_name, "===\n")
  print(summary(baseline_data[c("incarc_life", "AGEBEST", "OSTACCESSED", "Hbinary",
                               "sw_ever_partner", "sw_ever", "sw_ever_gender", "sw_30d", 
                               "inj_freq_days_30d", "YRSINJ")]))
  
  return(baseline_data)
}

# Apply to all datasets
baseline_all <- create_baseline_summary(kenya_clean, "All participants")
baseline_males <- create_baseline_summary(kenya_males, "Males only")
baseline_females <- create_baseline_summary(kenya_females, "Females only")

# Count incident HIV cases 
num_incident_cases_hiv <- sum(kenya_clean$hiv_test_rslt == "Positive", na.rm = TRUE)

# Calculate total person-time in years
total_person_years_hiv <- sum(kenya_clean$lead_time, na.rm = TRUE)

# Incidence rate per 100 person-years
incidence_rate_hiv <- (num_incident_cases_hiv / total_person_years_hiv) * 100

# Poisson 95% CI for incidence rate
ir_lower_hiv <- (qchisq(0.025, 2 * num_incident_cases_hiv) / 2) / total_person_years_hiv * 100
ir_upper_hiv <- (qchisq(0.975, 2 * (num_incident_cases_hiv + 1)) / 2) / total_person_years_hiv * 100

num_incident_cases_hiv
total_person_years_hiv
incidence_rate_hiv
ir_lower_hiv
ir_upper_hiv

# Function to calculate incidence rates and rate ratios
calculate_incidence <- function(data, exposure_var, reference_level = NULL) {
  
  # Calculate incidence by exposure groups
  incidence_by_exposure <- data %>%
    group_by(!!sym(exposure_var)) %>%
    summarise(
      num_incident_cases = sum(hiv_test_rslt == "Positive", na.rm = TRUE),
      total_person_years = sum(lead_time, na.rm = TRUE),
      incidence_rate = (num_incident_cases / total_person_years) * 100,
      # Poisson 95% CI for incidence rate
      ir_lower = (qchisq(0.025, 2 * num_incident_cases) / 2) / total_person_years * 100,
      ir_upper = (qchisq(0.975, 2 * (num_incident_cases + 1)) / 2) / total_person_years * 100,
      .groups = 'drop'
    )
  
  print(incidence_by_exposure)
  
  # Get unique exposure levels
  exposure_levels <- unique(data[[exposure_var]])
  exposure_levels <- exposure_levels[!is.na(exposure_levels)]
  
  # Determine reference and comparison levels
  if (is.null(reference_level)) {
    if (is.numeric(exposure_levels)) {
      ref_level <- min(exposure_levels)
      comp_level <- max(exposure_levels)
    } else {
      ref_level <- exposure_levels[1]
      comp_level <- exposure_levels[2]
    }
  } else {
    ref_level <- reference_level
    comp_level <- setdiff(exposure_levels, reference_level)[1]
  }
  
  # Extract values for rate ratio calculation
  cases_0 <- incidence_by_exposure %>% filter(!!sym(exposure_var) == ref_level) %>% pull(num_incident_cases)
  py_0    <- incidence_by_exposure %>% filter(!!sym(exposure_var) == ref_level) %>% pull(total_person_years)
  cases_1 <- incidence_by_exposure %>% filter(!!sym(exposure_var) == comp_level) %>% pull(num_incident_cases)
  py_1    <- incidence_by_exposure %>% filter(!!sym(exposure_var) == comp_level) %>% pull(total_person_years)
  

  # Calculate rate ratio and CI
  rate_0 <- cases_0 / py_0
  rate_1 <- cases_1 / py_1
  rate_ratio <- rate_1 / rate_0
  
  se_log_rr <- sqrt(1/cases_1 + 1/cases_0)
  log_rr <- log(rate_ratio)
  ci_lower <- exp(log_rr - 1.96 * se_log_rr)
  ci_upper <- exp(log_rr + 1.96 * se_log_rr)
  
  cat("\n", exposure_var, "Rate Ratio (", comp_level, "vs", ref_level, "):", 
      round(rate_ratio, 3), "95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3), "\n\n")
  
  return(incidence_by_exposure)
}

# apply incidence function
calculate_incidence(kenya_clean, "sw_30d")
calculate_incidence(kenya_clean, "sw_ever")
calculate_incidence(kenya_males, "sw_30d")
calculate_incidence(kenya_males, "sw_ever")
calculate_incidence(kenya_females, "sw_30d")
calculate_incidence(kenya_females, "sw_ever")

# Function to run adjusted Poisson regression models
run_poisson_models <- function(data, exposure_var) {
  
  # Model 1: Adjusted for homelessness and incarceration
  formula1 <- as.formula(paste("hiv_test_rslt == 'Positive' ~", exposure_var, "+ Hbinary + incarc_life + offset(log(lead_time))"))
  model_adj1 <- glm(formula1, family = poisson, data = data)
  
  # Get robust standard errors
  robust_se_adj1 <- sqrt(diag(vcovHC(model_adj1, type = "HC0")))
  
  # Extract coefficients and calculate IRRs with robust CIs
  coef_adj1 <- coef(model_adj1)
  irr_adj1 <- exp(coef_adj1)
  ci_lower_adj1 <- exp(coef_adj1 - 1.96 * robust_se_adj1)
  ci_upper_adj1 <- exp(coef_adj1 + 1.96 * robust_se_adj1)
  
  # Model summary
  cat("\n=== Adjusted Poisson Regression:", exposure_var, "(homelessness + incarceration) ===\n")
  print(summary(model_adj1))
  cat("\nIncidence Rate Ratios (IRR) with robust 95% CI:\n")
  results_adj1 <- data.frame(
    Variable = names(coef_adj1),
    IRR = round(irr_adj1, 3),
    CI_Lower = round(ci_lower_adj1, 3),
    CI_Upper = round(ci_upper_adj1, 3)
  )
  print(results_adj1)
  
  # Model 2: Fully adjusted (homelessness, incarceration, and injection variables)
  formula2 <- as.formula(paste("hiv_test_rslt == 'Positive' ~", exposure_var, "+ Hbinary + incarc_life + inj_freq_days_30d + YRSINJ + offset(log(lead_time))"))
  model_adj2 <- glm(formula2, family = poisson, data = data)
  
  # Get robust standard errors
  robust_se_adj2 <- sqrt(diag(vcovHC(model_adj2, type = "HC0")))
  
  # Extract coefficients and calculate IRRs with robust CIs
  coef_adj2 <- coef(model_adj2)
  irr_adj2 <- exp(coef_adj2)
  ci_lower_adj2 <- exp(coef_adj2 - 1.96 * robust_se_adj2)
  ci_upper_adj2 <- exp(coef_adj2 + 1.96 * robust_se_adj2)
  
  # Model summary
  cat("\n=== Fully Adjusted Poisson Regression:", exposure_var, "(all covariates) ===\n")
  print(summary(model_adj2))
  cat("\nIncidence Rate Ratios (IRR) with robust 95% CI:\n")
  results_adj2 <- data.frame(
    Variable = names(coef_adj2),
    IRR = round(irr_adj2, 3),
    CI_Lower = round(ci_lower_adj2, 3),
    CI_Upper = round(ci_upper_adj2, 3)
  )
  print(results_adj2)
  
  # Return both models
  return(list(model_adj1 = model_adj1, model_adj2 = model_adj2, 
              results_adj1 = results_adj1, results_adj2 = results_adj2))
}

# apply poisson function
sw_30d_models <- run_poisson_models(kenya_clean, "sw_30d")
sw_30d_male_models <- run_poisson_models(kenya_males, "sw_30d")
sw_30d_female_models <- run_poisson_models(kenya_females, "sw_30d")
sw_ever_models <- run_poisson_models(kenya_clean, "sw_ever")
sw_ever_male_models <- run_poisson_models(kenya_males, "sw_ever")
sw_ever_female_models <- run_poisson_models(kenya_females, "sw_ever")

# manual calculation for male IRs

# Calculate incidence by exposure groups
  incidence_by_exposure <- kenya_males %>%
    group_by(sw_30d) %>%
    summarise(
      num_incident_cases = sum(hiv_test_rslt == "Positive", na.rm = TRUE) + 0.5,
      total_person_years = sum(lead_time, na.rm = TRUE),
      incidence_rate = (num_incident_cases / total_person_years) * 100,
      # Poisson 95% CI for incidence rate
      ir_lower = (qchisq(0.025, 2 * num_incident_cases) / 2) / total_person_years * 100,
      ir_upper = (qchisq(0.975, 2 * (num_incident_cases + 1)) / 2) / total_person_years * 100,
      .groups = 'drop'
    )

  print(incidence_by_exposure)

# Calculate rate ratio and 95% CI
ref_group <- incidence_by_exposure %>% filter(sw_30d == "No")
exp_group <- incidence_by_exposure %>% filter(sw_30d == "Yes")

if(nrow(ref_group) > 0 & nrow(exp_group) > 0) {
  rate_ratio <- exp_group$incidence_rate / ref_group$incidence_rate
  
  # Calculate SE of log rate ratio
  se_log_rr <- sqrt(1/exp_group$num_incident_cases + 1/ref_group$num_incident_cases)
  
  # Calculate 95% CI for rate ratio
  ci_lower <- exp(log(rate_ratio) - 1.96 * se_log_rr)
  ci_upper <- exp(log(rate_ratio) + 1.96 * se_log_rr)
  
  cat("\nRate Ratio (Yes vs No):", round(rate_ratio, 3), 
      "95% CI:", round(ci_lower, 3), "-", round(ci_upper, 3), "\n")
} else {
  cat("Cannot calculate rate ratio - missing exposure groups\n")
}