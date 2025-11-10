############################################################
# master_macro_malaria.R
# Master script to run the malaria–TFP replication in order
#
# Order of scripts:
#   1. data_processing.R
#   2. main_analysis.R
#   3. A7_robustness_WDI_data_analysis.R
#   4. A8_A15_robustness_output_per_labour_auxiliary_check.R
#   5. A9_robustness_dengue.R
#   6. A10_robustness_HICs.R
#   7. A11_robustness_two_sample_MR.R
#   8. A12_robustness_institutional_quality.R
#   9. A13_causal_mediation_hci.R
#  10. A13_causal_mediation_life_expectancy.R
#  11. A14_auxiliary_check_age_group.R
#
# Run this file from the project root.
############################################################
library(tidyverse)  
library(readxl)
library(plm)
library(ivreg)
library(lmtest)
library(sandwich) 
library(data.table)
library(mediation)
library(broom)
rm(list = ls())

DATA_DIR <- "Data"    
OUTPUT_DIR <- "Output"    
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, showWarnings = FALSE)
data_path <- function(...) file.path(DATA_DIR, ...)

## Set working directory to the project root
# setwd("path/to/project/root")

# Start logging all warnings and messages
log_con <- file(file.path(OUTPUT_DIR, "replication_run_log.txt"), open = "wt")

sink(log_con, split = TRUE)
# warnings/messages (e.g. from readr, log(), etc.)
sink(log_con, type = "message")

cat("==== Starting replication ====\n")

cat("Running data_processing.R ...\n")
source("data_processing.R")

cat("Running main_analysis.R ...\n")
source("main_analysis.R")

cat("Running A7_robustness_WDI_data_analysis.R ...\n")
source("A7_robustness_WDI_data_analysis.R")

cat("Running A8_A15_robustness_output_per_labour_auxiliary_check.R ...\n")
source("A8_A15_robustness_output_per_labour_auxiliary_check.R")

cat("Running A9_robustness_dengue.R ...\n")
source("A9_robustness_dengue.R")

cat("Running A10_robustness_HICs.R ...\n")
source("A10_robustness_HICs.R")

cat("Running A11_robustness_two_sample_MR.R ...\n")
source("A11_robustness_two_sample_MR.R")

cat("Running A12_robustness_institutional_quality.R ...\n")
source("A12_robustness_institutional_quality.R")

cat("Running A13_causal_mediation_hci.R ...\n")
source("A13_causal_mediation_hci.R")

cat("Running A13_causal_mediation_life_expectancy.R ...\n")
source("A13_causal_mediation_life_expectancy.R")

cat("Running A14_auxiliary_check_age_group.R ...\n")
source("A14_auxiliary_check_age_group.R")

sink(type = "message")
sink()
close(log_con) 

cat("\nAll scripts completed successfully.\n")

