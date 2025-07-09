# load packages
pacman::p_load(dplyr, tidyr, readr, readxl, lubridate)

# set working directory
setwd("C:/Users/vl22683/OneDrive - University of Bristol/Documents/Publications/Sex work and risk of HIV and HCV/Emails to authors/Kenya data/data")

# load raw data
load("indat_25June25.RData")  
kenya_raw <- get(load("indat_25June25.RData"))
View(kenya_raw)
# Fixed version - handle NAs properly
kenya_clean <- kenya_raw %>%
  rename(
    incarc_life = F01,
    sexually_active_30d = C03,
    gender_partner_30d = C04,
    sw_ever_partner = C12,
    sw_ever = C28,
    sw_ever_gender = C28A,
    sw_30d = C30,
    inj_freq_days_30d = B02,
    inj_freq_num_30d = B02A
  ) %>%
  group_by(id) %>%
  mutate(id_seq = row_number()) %>%
  ungroup() %>%
  mutate(
    # convert to character
    sw_ever = as.character(sw_ever),
    sw_30d = as.character(sw_30d),
    
    # update sw_ever and sw_recent
    sw_ever = ifelse(sw_ever_partner == "Yes", "Yes", sw_ever),
    sw_ever = ifelse(sw_30d == "Yes", "Yes", sw_ever),
    sw_ever = ifelse(is.na(sw_ever), "No", sw_ever),
    sw_30d = ifelse(is.na(sw_30d), "No", sw_30d),
    sw_30d = ifelse(sw_ever == "No", "No", sw_30d),
    
    # msm
    msm_30d = case_when(
        GENDER == "Female" ~ NA_real_,
        gender_partner_30d == "Male Participant, Male Main Partner" & !is.na(gender_partner_30d) ~ 1,
        sw_ever_gender == "Men" & !is.na(sw_ever_gender) ~ 1,
        TRUE ~ 0
)
  )

# convert character variables to factors
kenya_clean <- kenya_clean %>%
  mutate(
    incarc_life = as.factor(incarc_life),
    sexually_active_30d = as.factor(sexually_active_30d),
    gender_partner_30d = as.factor(gender_partner_30d),
    sw_ever_partner = as.factor(sw_ever_partner),
    sw_ever = as.factor(sw_ever),
    sw_ever_gender = as.factor(sw_ever_gender),
    sw_30d = as.factor(sw_30d),
    msm_30d = as.factor(msm_30d)
)
  
# lead of HIVSTATUS
kenya_clean <- kenya_clean %>%
  group_by(id) %>%
  mutate(hiv_test_rslt = lead(HIVSTATUS)) %>%
  ungroup() %>%
  filter(!is.na(hiv_test_rslt))

# save cleaned data
write_csv(kenya_clean, "kenya_clean.csv")