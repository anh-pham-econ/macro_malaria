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

sample_path <- data_path("sample_list.csv")
country_list <- read.csv(sample_path, header = TRUE)

path_macro_malaria <- data_path("LLMICs_all_data.csv")
shift_share_iv <- read.csv(path_macro_malaria, header = TRUE)

tfp_data <- pdata.frame(shift_share_iv, index = c("Country.Name", "year"))


country_full_sample <- unique(country_list$Country.Name)

compute_wald_components <- function(data, countries) {
  split_countries <- sample(countries, size = floor(length(countries) / 2))
  sample1 <- data %>% filter(Country.Name %in% split_countries)
  sample2 <- data %>% filter(!Country.Name %in% split_countries)
  
  # Collapse to country averages
  sample1_avg <- sample1 %>% group_by(Country.Name) %>%
    summarise(malaria = mean(case_incidence_log, na.rm = TRUE),
              shiftshare = mean((hbs_olr_o), na.rm = TRUE))
  
  sample2_avg <- sample2 %>% group_by(Country.Name) %>%
    summarise(ln_TFP = mean(tfp_log, na.rm = TRUE),
              shiftshare = mean((hbs_olr_o), na.rm = TRUE))
  
  first_stage <- lm(malaria ~ shiftshare, data = sample1_avg)
  second_stage <- lm(ln_TFP ~ shiftshare, data = sample2_avg)
  
  a_hat <- coef(first_stage)["shiftshare"]
  b_hat <- coef(second_stage)["shiftshare"]
  
  if (is.na(a_hat) || is.na(b_hat)) return(rep(NA, 3))
  
  return(c(wald = b_hat / a_hat, first_stage = a_hat, second_stage = b_hat))
}


cat("\n==============================\n")
cat("Table A11. Two-Sample Mendelian Randomisation Estimates \n")
cat("==============================\n\n")

set.seed(111)
wald_matrix <- replicate(1000, compute_wald_components(tfp_data, country_full_sample))
wald_df <- as.data.frame(t(na.omit(wald_matrix)))  # transpose and clean

mean_wald <- mean(wald_df$wald)
se_wald <- sd(wald_df$wald)
ci_wald <- quantile(wald_df$wald, probs = c(0.025, 0.975))

mean_first <- mean(wald_df$first_stage)
se_first <- sd(wald_df$first_stage)
ci_first <- quantile(wald_df$first_stage, probs = c(0.025, 0.975))

mean_second <- mean(wald_df$second_stage)
se_second <- sd(wald_df$second_stage)
ci_second <- quantile(wald_df$second_stage, probs = c(0.025, 0.975))

cat("Mean first-stage (shiftshare → malaria):", mean_first, "\n")
cat("Standard error (first-stage):", se_first, "\n")
cat("95% CI (first-stage):", ci_first[1], "to", ci_first[2], "\n\n")

cat("Mean second-stage (shiftshare → TFP):", mean_second, "\n")
cat("Standard error (second-stage):", se_second, "\n")
cat("95% CI (second-stage):", ci_second[1], "to", ci_second[2], "\n\n")

cat("Mean Wald estimate:", mean_wald, "\n")
cat("Standard error (Wald):", se_wald, "\n")
cat("95% CI (Wald):", ci_wald[1], "to", ci_wald[2], "\n\n")
