# load packages
pacman::p_load(dplyr, tidyr, readr, readxl, lubridate, sandwich, lmtest)

# set working directory
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Publications/Sex work and risk of HIV and HCV/Emails to authors/Kenya data/data")

# load hiv data
kenya_clean <- read_csv("kenya_clean.csv")
View(kenya_clean)

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

# filter to sequence = 1 and create summary
baseline_data <- kenya_clean %>% filter(id_seq == 1)

# summary table
summary(baseline_data[c("incarc_life", "AGEBEST", "OSTACCESSED", "Hbinary",
                       "sw_ever_partner", "sw_ever", "sw_ever_gender", "sw_30d", 
                       "inj_freq_days_30d", "YRSINJ")])

# Count incident HIV cases (assuming each "Positive" is an incident case)
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
  
  # Check if we have data for both groups
  if (length(cases_0) == 0 || length(cases_1) == 0 || cases_0 == 0 || cases_1 == 0) {
    cat("Cannot calculate rate ratio - insufficient data\n")
    return(incidence_by_exposure)
  }
  
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
calculate_incidence(kenya_clean, "msm_30d")
calculate_incidence(kenya_clean, "sw_30d")
calculate_incidence(kenya_clean, "sw_ever")
kenya_males <- kenya_clean %>% filter(GENDER == "Male")
calculate_incidence(kenya_males, "sw_30d")
calculate_incidence(kenya_males, "sw_ever")
kenya_females <- kenya_clean %>% filter(GENDER == "Female")
calculate_incidence(kenya_females, "sw_30d")
calculate_incidence(kenya_females, "sw_ever")
