# library(tidyverse)  
# library(readxl)
# library(plm)
# library(ivreg)
# library(lmtest)
# library(sandwich) 
# library(data.table)
# library(mediation)
# library(broom)
# 
# DATA_DIR <- "Data"    
# OUTPUT_DIR <- "Output"    
# if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, showWarnings = FALSE)
# data_path <- function(...) file.path(DATA_DIR, ...)
##Set working directory to the project root
# setwd("path/to/project/root")


path_HICs_macro_malaria <- data_path("HICs_macro_data.csv")
HICs <- read.csv(path_HICs_macro_malaria, header = TRUE)

sickle_cell_path <- data_path("Sickle_cell","sickle.csv")
sickle_cell <- read.csv(sickle_cell_path, header = TRUE)

malaria_incidence_path <- data_path("World_malaria_report_2023_Annex","malaria_incidence.csv")
malaria_incidence <- read.csv(malaria_incidence_path, header = TRUE)

olr_annual_path <- data_path("ENSO","olr_annual.csv")
olr_annual <- read.csv(olr_annual_path, header = TRUE)



HICs_pdata <- pdata.frame(HICs, index = c("Country.Name", "year"))
HICs_panel_model_fixed <- plm(gdp_log ~  k_log + employment_log , data = HICs_pdata, model = "within")
summary(HICs_panel_model_fixed)

HICs_intercept_value <- as.data.frame(summary(fixef(HICs_panel_model_fixed)))
HICs_intercept_value <- HICs_intercept_value %>%
  rownames_to_column(var = "Country.Name") 

HICs <- HICs_intercept_value %>%
  select(Country.Name, Estimate) %>%
  left_join(HICs, by = "Country.Name")

colnames(HICs)[2] <- 'a_log'
HICs_res <- as.data.frame(resid(HICs_panel_model_fixed))
colnames(HICs_res)[1] <- 'residuals_log'
HICs_res <- HICs_res %>%
  rownames_to_column(var = "Country.Name-year") 
HICs_res <- separate(HICs_res, 'Country.Name-year', into = c("Country.Name", "year"), sep = "-")
HICs_res$year <- as.integer(HICs_res$year)
HICs <- HICs_res %>%
  full_join(HICs, by = c("Country.Name", "year"))

HICs$tfp_log = HICs$a_log + HICs$residuals_log


HICs_HbS_panel <- sickle_cell %>%
  inner_join(HICs, 'Country.Code')

HICs_HbS_panel <- malaria_incidence %>%
  select(Country.Name, year, case_incidence, cases_log, lag_case_log, death_incidence, case_incidence_log, lag_case_incidence_log, death_incidence_log) %>%
  inner_join(HICs_HbS_panel, c("Country.Name", "year"))

HICs_HbS_panel<- HICs_HbS_panel %>%
  filter(year >= '2000') %>%
  filter(year < '2020') %>%
  group_by(Country.Code) %>%
  ungroup()

HICs_shift_share_panel <- HICs_HbS_panel %>%
  left_join(olr_annual, by =  'year') %>%
  mutate(
    hbs_olr_o  = log(results_weighted) * log(olr_o)
  ) %>%
  arrange(Country.Code, year)



cat("\n==============================\n")
cat("Table A10. Unbalanced Panel Estimates of Malaria Impact on TFP in UMHICs \n")
cat("==============================\n\n")

HICs_shift_share_pdata <- pdata.frame(HICs_shift_share_panel, index = c("Country.Name", "year"))


HICs_first_stage_shiftshare <- plm(case_incidence_log ~ (hbs_olr_o) + agr_share_log + trade_openness_log + cpi_log + edu_log , data = HICs_shift_share_pdata, model = "within")

cat("\n==============================\n")
cat("Table A10. First stage \n")
cat("==============================\n\n")
print(summary(HICs_first_stage_shiftshare))

HICs_iv_shiftshare <- ivreg::ivreg(tfp_log ~ case_incidence_log  + agr_share_log + trade_openness_log + cpi_log  +  edu_log  + factor(Country.Name) |
                              (hbs_olr_o) + agr_share_log + trade_openness_log + cpi_log  +  edu_log  + factor(Country.Name), data = HICs_shift_share_pdata)

cat("\n==============================\n")
cat("Table A10. 2SLS \n")
cat("==============================\n\n")
print(summary(HICs_iv_shiftshare))




