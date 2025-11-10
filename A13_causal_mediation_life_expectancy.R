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

life_expectancy_file_path <- data_path("life-expectancy","life-expectancy.csv")
life_expectancy_file_path_data <- read.csv(life_expectancy_file_path, header = TRUE)


names(life_expectancy_file_path_data)[2] <- "Country.Code"
names(life_expectancy_file_path_data)[3] <- "year"
names(life_expectancy_file_path_data)[4] <- "life_expectancy"
life_expectancy_file_path_data[, 1] <- NULL

life_expectancy_file_path_data$year <- as.integer(life_expectancy_file_path_data$year)
life_expectancy_file_path_data$life_expectancy_log = log(life_expectancy_file_path_data$life_expectancy)


life_expectancy <- life_expectancy_file_path_data %>%
  inner_join(shift_share_iv, c("Country.Code", "year"))

tfp_data_life_expectancy <- pdata.frame(life_expectancy, index = c("Country.Code", "year"))

# Mediator model: how malaria affects life expectancy
model.m <- lm(life_expectancy ~ case_incidence_log + agr_share_log+ trade_openness_log + cpi_log+ edu_log, data = tfp_data_life_expectancy)
# Outcome model: how life expectancy and malaria jointly affect TFP
model.y <- lm(tfp_log ~ case_incidence_log + life_expectancy + agr_share_log + trade_openness_log + cpi_log+ edu_log, data = tfp_data_life_expectancy)
# Run mediation with 1000 simulations (nonparametric bootstrap)
set.seed(111)
med.out <- mediate(model.m, model.y,
                   treat = "case_incidence_log",
                   mediator = "life_expectancy",
                   robustSE = TRUE,
                   sims = 1000)

cat("\n==============================\n")
cat("Table A13. Mediator: Life expectancy \n")
cat("==============================\n\n")

print(summary(med.out))


