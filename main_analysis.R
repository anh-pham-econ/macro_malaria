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

#########################
# Main TFP regressions (FE + IV)
#########################

tfp_data <- pdata.frame(shift_share_iv, index = c("Country.Name", "year"))

# FE, OLS
tfp_panel_model_fixed <- plm(
  tfp_log ~ case_incidence_log + agr_share_log + trade_openness_log +
    cpi_log + edu_log,
  data  = tfp_data,
  model = "within"
)
summary(tfp_panel_model_fixed)
coeftest(tfp_panel_model_fixed,
         vcovHC(tfp_panel_model_fixed, type = "HC1", cluster = "group"))

# Get list of Countries Included in the Sample
intercept_value <- as.data.frame(summary(fixef(tfp_panel_model_fixed))) %>%
  rownames_to_column(var = "Country.Name")

country_list <- intercept_value %>%
  left_join(tfp_data %>% select(Country.Name, Country.Code,Region, IncomeGroup),
            by = "Country.Name") %>%
  distinct()
country_list <- country_list[,c(1,6,7,8)]

cat("\n==============================\n")
cat("Table A2. List of malaria-endemic LLMICs included in the estimation \n")
cat("==============================\n\n")
print(country_list)

#write.csv(country_list, file = "Data/sample_list.csv", row.names = FALSE)

cat("\n==============================\n")
cat("Table 2. FE \n")
cat("==============================\n\n")
print(coeftest(tfp_panel_model_fixed,
               vcovHC(tfp_panel_model_fixed, type = "HC1", cluster = "group"))
)

# RE
tfp_panel_model_random <- plm(
  tfp_log ~ case_incidence_log + agr_share_log + trade_openness_log +
    cpi_log + edu_log,
  data  = tfp_data,
  model = "random"
)
summary(tfp_panel_model_random)
coeftest(tfp_panel_model_random,
         vcovHC(tfp_panel_model_random, type = "HC1", cluster = "group"))

cat("\n==============================\n")
cat("Table 2. RE \n")
cat("==============================\n\n")
print(coeftest(tfp_panel_model_random,
               vcovHC(tfp_panel_model_random, type = "HC1", cluster = "group"))
)
#########################
# Lagged-IV
#########################
lagged_iv_model_pooled <- ivreg::ivreg(
  tfp_log ~ case_incidence_log + agr_share_log + trade_openness_log +
    cpi_log + edu_log + factor(Country.Name) |
    lag_case_incidence_log + agr_share_log + trade_openness_log +
    cpi_log + edu_log + factor(Country.Name),
  data = tfp_data
)
summary(lagged_iv_model_pooled, diagnostics = TRUE)
coeftest(lagged_iv_model_pooled,
         vcovHC(lagged_iv_model_pooled, type = "HC1"))


cat("\n==============================\n")
cat("Table 2. 2SLS Lag IV \n")
cat("==============================\n\n")
print(coeftest(lagged_iv_model_pooled,
               vcovHC(lagged_iv_model_pooled, type = "HC1"))
)
# KP rk F and partial R²
clust <- tfp_data$Country.Name

first_stage_lagged_iv <- lm(
  case_incidence_log ~ lag_case_incidence_log + agr_share_log +
    trade_openness_log + cpi_log + edu_log + factor(Country.Name),
  data = tfp_data
)
summary(first_stage_lagged_iv)

cat("\n==============================\n")
cat("Table A5. Lag IV First stage \n")
cat("==============================\n\n")
print(summary(first_stage_lagged_iv))

t_lag_case <- coeftest(
  first_stage_lagged_iv,
  vcovCL(first_stage_lagged_iv, type = "HC1", cluster = clust)
)

KP_rk_F_lag_IV <- (as.numeric(t_lag_case["lag_case_incidence_log", "t value"]))^2
cat("F-statistic \n")
print(KP_rk_F_lag_IV)

first_stage_no_lag <- lm(
  case_incidence_log ~ agr_share_log + trade_openness_log + cpi_log +
    edu_log + factor(Country.Name),
  data = tfp_data
)
summary(first_stage_no_lag)

SSR_with_lag <- sum(residuals(first_stage_lagged_iv)^2)
SSR_no_lag   <- sum(residuals(first_stage_no_lag)^2)
partial_R2_lag_IV   <- (SSR_no_lag - SSR_with_lag) / SSR_no_lag
partial_R2_lag_IV
cat("Partial R-squared \n")
print(partial_R2_lag_IV)

#########################
# Shift–share IV (HbS × OLR)
#########################

shift_share_iv_sickle <- ivreg::ivreg(
  tfp_log ~ case_incidence_log + agr_share_log + trade_openness_log +
    cpi_log + edu_log + factor(Country.Name) |
    hbs_olr_o + agr_share_log + trade_openness_log +
    cpi_log + edu_log + factor(Country.Name),
  data = tfp_data
)
summary(shift_share_iv_sickle, diagnostics = TRUE)
coeftest(shift_share_iv_sickle,
         vcovHC(shift_share_iv_sickle, type = "HC1"))

cat("\n==============================\n")
cat("Table 2. 2SLS Shift-share IV \n")
cat("==============================\n\n")
print(coeftest(shift_share_iv_sickle,
           vcovHC(shift_share_iv_sickle, type = "HC1"))
  )
# First stage (shift–share) and partial R²
first_stage_shift_share <- lm(
  case_incidence_log ~ hbs_olr_o + agr_share_log + trade_openness_log +
    cpi_log + edu_log + factor(Country.Name),
  data = tfp_data
)
summary(first_stage_shift_share)


cat("\n==============================\n")
cat("Table A5. Shift-share IV First stage \n")
cat("==============================\n\n")
print(summary(first_stage_shift_share))

t_hbs_olr <- coeftest(
  first_stage_shift_share,
  vcovCL(first_stage_shift_share, type = "HC1", cluster = clust)
)

KP_rk_F_shift_share <- (as.numeric(t_hbs_olr["hbs_olr_o", "t value"]))^2
KP_rk_F_shift_share
cat("F-statistic \n")
print(KP_rk_F_shift_share)

first_stage_no_shift_share <- lm(
  case_incidence_log ~ agr_share_log + trade_openness_log + cpi_log +
    edu_log + factor(Country.Name),
  data = tfp_data
)
summary(first_stage_no_shift_share)

SSR_withZ <- sum(residuals(first_stage_shift_share)^2)
SSR_noZ   <- sum(residuals(first_stage_no_shift_share)^2)
partial_R2_shift_share <- (SSR_noZ - SSR_withZ) / SSR_noZ
partial_R2_shift_share

cat("Partial R-squared \n")
print(partial_R2_shift_share)


#########################
# 12. Exposure-share diffuseness (HbS)
#########################

expo <- tfp_data %>%
  as.data.frame() %>%
  group_by(Country.Name) %>%
  summarise(hbs_prev = mean(results_weighted, na.rm = TRUE),
            .groups = "drop") %>%
  mutate(share = hbs_prev / sum(hbs_prev, na.rm = TRUE))

top_exposure_share <- max(expo$share, na.rm = TRUE)
HHI_exposure       <- sum(expo$share^2, na.rm = TRUE)
effective_exposers <- 1 / HHI_exposure

srt <- sort(expo$share, decreasing = TRUE)
top1_share <- srt[1]
top3_share <- sum(srt[seq_len(min(3, length(srt)))])
top5_share <- sum(srt[seq_len(min(5, length(srt)))])


cat("\n==============================\n")
cat("Table A4. Concentration Measures of Exposure Components \n")
cat("==============================\n\n")
print(c(
  top_exposure_share = top_exposure_share,
  HHI_exposure       = HHI_exposure,
  effective_exposers = effective_exposers,
  top1_share         = top1_share,
  top3_share         = top3_share,
  top5_share         = top5_share
))


#########################
# Shift diffuseness (OLR)
#########################

shift_tbl <- tfp_data %>%
  group_by(year) %>%
  summarise(mean_olr = mean(olr_o, na.rm = TRUE),
            .groups = "drop")%>%
  mutate(abs_olr = abs(mean_olr),
         share = abs_olr / sum(abs_olr, na.rm = TRUE))

top_shift_share <- max(shift_tbl$share, na.rm = TRUE)
HHI_shift <- sum(shift_tbl$share^2, na.rm = TRUE)
effective_shocks <- 1 / HHI_shift

srt <- sort(shift_tbl$share, decreasing = TRUE)
top1 <- srt[1]
top3 <- sum(srt[seq_len(min(3, length(srt)))])
top5 <- sum(srt[seq_len(min(5, length(srt)))])


cat("\n==============================\n")
cat("Table A4. Concentration Measures of Shift Components \n")
cat("==============================\n\n")

print(c(top_shift_share = top_shift_share,
  HHI_shift = HHI_shift,
  effective_shocks = effective_shocks,
  top1_share = top1,
  top3_share = top3,
  top5_share = top5))



#########################
# Combined IV
#########################

combined_iv_model <- ivreg::ivreg(
  tfp_log ~ case_incidence_log + agr_share_log + trade_openness_log +
    cpi_log + edu_log + factor(Country.Name) |
    hbs_olr_o + lag_case_incidence_log + log(olr_o) + agr_share_log +
    trade_openness_log + cpi_log + edu_log + factor(Country.Name),
  data = tfp_data
)


cat("\n==============================\n")
cat("Table A6. Panel 2SLS Estimates of Malaria Impact on TFP using Combined instruments \n")
cat("==============================\n\n")

print(summary(combined_iv_model, diagnostics = TRUE,
        vcov = vcovHC(combined_iv_model, type = "HC1")))


#########################
# Sensitivity to Factor Shares
#########################

factor_shares <- shift_share_iv
factor_shares$tfp_log_s1 <- factor_shares$gdp_log - 0.35*factor_shares$k_log - 0.65*factor_shares$employment_log
factor_shares$tfp_log_s2 <- factor_shares$gdp_log - 0.30*factor_shares$k_log - 0.65*factor_shares$employment_log
factor_shares$tfp_log_s3 <- factor_shares$gdp_log - 0.4*factor_shares$k_log - 0.65*factor_shares$employment_log
factor_shares_data <- pdata.frame(factor_shares, index = c("Country.Name", "year"))

#s1: \alpha = 0.35, \beta = 0.65
shift_share_iv_sickle_s1 <- ivreg::ivreg(tfp_log_s1 ~ case_incidence_log + agr_share_log + trade_openness_log + cpi_log  + edu_log +factor(Country.Name) |
                                           (hbs_olr_o)  + agr_share_log + trade_openness_log + cpi_log  + edu_log +factor(Country.Name),
                                         data = factor_shares_data)

cat("\n==============================\n")
cat("Table A17. Robustness of Malaria Impact on TFP to Alternative Factor Shares, Alpha = 0.35, Beta= 0.65 \n")
cat("==============================\n\n")

print(summary(shift_share_iv_sickle_s1, diagnostics = TRUE,
        vcov = vcovHC(shift_share_iv_sickle_s1, type = "HC1")))


#s2: \alpha = 0.3, \beta = 0.65
shift_share_iv_sickle_s2 <- ivreg::ivreg(tfp_log_s2 ~ case_incidence_log + agr_share_log + trade_openness_log + cpi_log  + edu_log +factor(Country.Name) |
                                           (hbs_olr_o)  + agr_share_log + trade_openness_log + cpi_log  + edu_log +factor(Country.Name),
                                         data = factor_shares_data)

cat("\n==============================\n")
cat("Table A17. Robustness of Malaria Impact on TFP to Alternative Factor Shares, Alpha = 0.30, Beta= 0.65 \n")
cat("==============================\n\n")
print(summary(shift_share_iv_sickle_s2, diagnostics = TRUE,
        vcov = vcovHC(shift_share_iv_sickle_s2, type = "HC1")))

#s3: \alpha = 0.4, \beta = 0.65
shift_share_iv_sickle_s3 <- ivreg::ivreg(tfp_log_s3 ~ case_incidence_log + agr_share_log + trade_openness_log + cpi_log  + edu_log +factor(Country.Name) |
                                           (hbs_olr_o)  + agr_share_log + trade_openness_log + cpi_log  + edu_log +factor(Country.Name),
                                         data = factor_shares_data)


cat("\n==============================\n")
cat("Table A17. Robustness of Malaria Impact on TFP to Alternative Factor Shares, Alpha = 0.40, Beta= 0.65 \n")
cat("==============================\n\n")

print(summary(shift_share_iv_sickle_s3, diagnostics = TRUE,
        vcov = vcovHC(shift_share_iv_sickle_s3, type = "HC1")))

vars <- factor_shares_data[, c("tfp_log", "tfp_log_s1", "tfp_log_s2", "tfp_log_s3")]
corr_mat <- round(cor(vars, use = "complete.obs"), 4)

cat("\n==============================\n")
cat("Table A16. Correlation Matrix of TFP Measures under Alternative Factor Shares \n")
cat("==============================\n\n")
print(corr_mat)
# 
# write.csv(corr_mat,
#           file = file.path(OUTPUT_DIR, "Table_A16_tfp_correlation.csv"),
#           row.names = TRUE)
