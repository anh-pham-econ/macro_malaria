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


sickle_cell_path <- data_path("Sickle_cell","sickle.csv")
sickle_cell <- read.csv(sickle_cell_path, header = TRUE)

path_macro <- data_path("all_countries_macro_data.csv")
all_countries_dengue <- read.csv(path_macro, header = TRUE)


dengue_path_1 <- data_path("VBDs/Dengue/panel","filtered_data_EURO_1755101430504.csv")
dengue_path_2 <- data_path("VBDs/Dengue/panel","filtered_data_EMRO_1755101384156.csv")
dengue_path_3 <- data_path("VBDs/Dengue/panel","filtered_data_PAHO_1755101283185.csv")
dengue_path_4 <- data_path("VBDs/Dengue/panel","filtered_data_SEARO_1755101247238.csv")
dengue_path_5 <- data_path("VBDs/Dengue/panel","filtered_data_WPRO_1755101324612.csv")

panel_dengue_paths <- c(dengue_path_1, dengue_path_2, dengue_path_3, dengue_path_4, dengue_path_5)

process_panel_dengue_file <- function(path) {
  df <- read.csv(path, header = TRUE)
  names(df) <- gsub("ISO_A0", "Country.code", names(df))
  df %>%
    group_by(Country.code, Year) %>%
    summarise(cases = sum(dengue_total, na.rm = TRUE), .groups = "drop")
}

total_dengue_panel <- panel_dengue_paths %>%
  lapply(process_panel_dengue_file) %>%
  bind_rows() %>%
  group_by(Country.code, Year) %>%
  summarise(dengue_cases = sum(cases, na.rm = TRUE), .groups = "drop")

total_dengue_panel$dengue_cases_log = log(total_dengue_panel$dengue_cases)

total_dengue_panel[total_dengue_panel == -Inf] <- NA
total_dengue_panel <- na.omit(total_dengue_panel)

dengue_tfpdata <- pdata.frame(all_countries_dengue, index = c("Country.Name", "year"))
dengue_panel_model_fixed <- plm(gdp_log ~  k_log + employment_log , data = dengue_tfpdata, model = "within")


dengue_intercept_value <- as.data.frame(summary(fixef(dengue_panel_model_fixed)))
dengue_intercept_value <- dengue_intercept_value %>%
  rownames_to_column(var = "Country.Name") 

all_countries_dengue <- dengue_intercept_value %>%
  select(Country.Name, Estimate) %>%
  left_join(all_countries_dengue, by = "Country.Name")

colnames(all_countries_dengue)[2] <- 'a_log'
dengue_res <- as.data.frame(resid(dengue_panel_model_fixed))
colnames(dengue_res)[1] <- 'residuals_log'
dengue_res <- dengue_res %>%
  rownames_to_column(var = "Country.Name-year") 
dengue_res <- separate(dengue_res, 'Country.Name-year', into = c("Country.Name", "year"), sep = "-")
dengue_res$year <- as.integer(dengue_res$year)
all_countries_dengue <- dengue_res %>%
  full_join(all_countries_dengue, by = c("Country.Name", "year"))

all_countries_dengue$tfp_log = all_countries_dengue$a_log + all_countries_dengue$residuals_log
all_countries_dengue <- sickle_cell %>%
  inner_join(all_countries_dengue, 'Country.Code')

all_countries_dengue<- all_countries_dengue %>%
  filter(year >= '2000') %>%
  filter(year < '2020') %>%
  group_by(Country.Code) %>%
  ungroup()

dengue_shift_share_panel <- total_dengue_panel %>%
  inner_join(
    all_countries_dengue,
    by = c("Country.code" = "Country.Code", "Year" = "year")
  )


path_olr <- data_path("ENSO", "olr_annual.csv")
olr_annual <- read.csv(path_olr, header = TRUE)


dengue_shift_share_panel <- dengue_shift_share_panel %>%
  left_join(olr_annual, by =  c("Year" = "year")) %>%
  mutate(
    hbs_olr_o  = log(results_weighted) * log(olr_o)
  ) %>%
  arrange(Country.code, Year)


cat("\n==============================\n")
cat("Table A9. Unbalanced Panel Estimates of Dengue Impact on TFP \n")
cat("==============================\n\n")


dengue_shift_share_pdata <- pdata.frame(dengue_shift_share_panel, index = c("Country.Name", "Year"))

dengue_first_stage_sickle <- plm(dengue_cases_log ~ results_weighted + agr_share_log + trade_openness_log + cpi_log + edu_log, data = dengue_shift_share_pdata, model = "pooling")

cat("\n==============================\n")
cat("Table A9. First stage Estimates of HbS frequency on Dengue \n")
cat("==============================\n\n")
print(summary(dengue_first_stage_sickle))

dengue_first_stage_shiftshare <- plm(dengue_cases_log ~ (hbs_olr_o) + agr_share_log + trade_openness_log + cpi_log + edu_log , data = dengue_shift_share_pdata, model = "within")

cat("\n==============================\n")
cat("Table A9. First stage Estimates of Shift-share IV on Dengue \n")
cat("==============================\n\n")
print(summary(dengue_first_stage_shiftshare))

dengue_shiftshare_iv <- ivreg::ivreg(
  tfp_log ~ dengue_cases_log + agr_share_log + trade_openness_log +
    cpi_log + edu_log + factor(Country.Name) |
    hbs_olr_o + agr_share_log + trade_openness_log +
    cpi_log + edu_log + factor(Country.Name),
  data = dengue_shift_share_pdata
)
summary(dengue_shiftshare_iv, diagnostics = TRUE)

cat("\n==============================\n")
cat("Table A9. 2SLS Estimates of Dengue on TFP \n")
cat("==============================\n\n")
print(coeftest(dengue_shiftshare_iv,
         vcovHC(dengue_shiftshare_iv, type = "HC1")))

