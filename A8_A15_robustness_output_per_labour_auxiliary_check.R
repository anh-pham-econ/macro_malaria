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


path_macro_malaria <- data_path("LLMICs_all_data.csv")
shift_share_iv <- read.csv(path_macro_malaria, header = TRUE)

### A8_Alternative Productivity Measure Y/L

cat("\n==============================\n")
cat("Table A8. Alternative Productivity Measure Y/L \n")
cat("==============================\n\n")

tfp_data <- pdata.frame(shift_share_iv, index = c("Country.Name", "year"))
clust <- tfp_data$Country.Name

tfp_panel_model_fixed_robustness <- plm(log(Y_L) ~ case_incidence_log + agr_share_log + trade_openness_log + cpi_log + edu_log, data = tfp_data, model = "within")
summary(tfp_panel_model_fixed_robustness)

cat("\n==============================\n")
cat("Table A8. FE \n")
cat("==============================\n\n")
print(coeftest(tfp_panel_model_fixed_robustness,vcovHC(tfp_panel_model_fixed_robustness, type = 'HC1', cluster = 'group'))
)

tfp_panel_model_random_robustness <- plm(log(Y_L) ~ case_incidence_log + agr_share_log + trade_openness_log + cpi_log + edu_log , data = tfp_data, model = "random")
summary(tfp_panel_model_random_robustness)

cat("\n==============================\n")
cat("Table A8. RE \n")
cat("==============================\n\n")
print(coeftest(tfp_panel_model_random_robustness,vcovHC(tfp_panel_model_random_robustness, type = 'HC1', cluster = 'group')))

lagged_iv_model_pooled_robustness <- ivreg::ivreg(log(Y_L) ~ case_incidence_log + agr_share_log + trade_openness_log + cpi_log  + edu_log + factor(Country.Name) |
                                             lag_case_incidence_log  + agr_share_log + trade_openness_log + cpi_log  + edu_log  + factor(Country.Name) ,
                                           data = tfp_data)
summary(lagged_iv_model_pooled_robustness)

cat("\n==============================\n")
cat("Table A8. 2SLS Lag IV \n")
cat("==============================\n\n")
print(coeftest(lagged_iv_model_pooled_robustness,vcovHC(lagged_iv_model_pooled_robustness, type = 'HC1', cluster = clust)))


shift_share_iv_sickle_robustness <- ivreg::ivreg(log(Y_L) ~ case_incidence_log + agr_share_log + trade_openness_log + cpi_log  + edu_log+ factor(Country.Name) |
                                            (hbs_olr_o)  + agr_share_log + trade_openness_log + cpi_log  + edu_log + factor(Country.Name) ,
                                          data = tfp_data)
summary(shift_share_iv_sickle_robustness)

cat("\n==============================\n")
cat("Table A8. 2SLS Shift-share IV \n")
cat("==============================\n\n")
print(coeftest(shift_share_iv_sickle_robustness,vcovHC(shift_share_iv_sickle_robustness, type = 'HC1', cluster = clust))
)

### A15_auxiliary_check_capital_stock K

cat("\n==============================\n")
cat("\n====Auxiliary check===========\n")
cat("Table A15. Unbalanced Panel Estimates of Malaria Impact on Capital stocks \n")
cat("==============================\n\n")

shift_share_iv_sickle_robustness <- ivreg::ivreg(k_log ~ case_incidence_log + agr_share_log + trade_openness_log + cpi_log  + edu_log+ factor(Country.Name) |
                                                   (hbs_olr_o)  + agr_share_log + trade_openness_log + cpi_log  + edu_log + factor(Country.Name) ,
                                                 data = tfp_data)
print(summary(shift_share_iv_sickle_robustness))

