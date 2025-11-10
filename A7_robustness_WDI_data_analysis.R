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


# Load WB macro data ------------------

capital_alt <- read.csv(
  data_path("Capital_constant_2015", "gfcf_to_capital.csv"),
  header = TRUE
)

gdp_alt <- read.csv(
  data_path("GDP_PPP_constant_2017", "API_NY.GDP.MKTP.PP.KD_DS2_en_csv_v2_11551.csv"),
  header = TRUE, skip = 4
)

cpi_alt <- read.csv(
  data_path("GDP_deflator", "API_FP.CPI.TOTL.ZG_DS2_en_csv_v2_77.csv"),
  header = TRUE, skip = 4
)

employment_rate_alt <- read.csv(
  data_path("employment_to_pop", "API_SL.EMP.TOTL.SP.ZS_DS2_en_csv_v2_6303676.csv"),
  header = TRUE, skip = 4
)

population_alt <- read.csv(
  data_path("population", "API_SP.POP.1564.TO_DS2_en_csv_v2_169.csv"),
  header = TRUE, skip = 4
)

trade_openness_alt <- read.csv(
  data_path("trade_openess", "API_NE.TRD.GNFS.ZS_DS2_en_csv_v2_466.csv"),
  header = TRUE, skip = 4
)

agr_share_alt <- read.csv(
  data_path("Agriculture_constant_2015_USD", "agr_share.csv"),
  header = TRUE, skip = 4
)

edu_alt <- read.csv(
  data_path("Edu", "edu_series_completion_upper_secondary.csv"),
  header = TRUE
)

edu_alt <- read.csv(
  data_path("Edu", "edu_series_completion_upper_secondary.csv"),
  header = TRUE
)

income_group_alt <- read.csv(
  data_path("Income_group", "Metadata_Country_API_NE.GDI.FTOT.CD_DS2_en_csv_v2_6303108.csv"),
  header = TRUE
)

columns_to_remove <- c("Indicator.Code", "Indicator.Name", "X")

gdp_alt[, columns_to_remove]            <- NULL
cpi_alt[, columns_to_remove]            <- NULL
employment_rate_alt[, columns_to_remove] <- NULL
population_alt[, columns_to_remove]     <- NULL
trade_openness_alt[, columns_to_remove] <- NULL
agr_share_alt[, columns_to_remove]      <- NULL
edu_alt[, c("Series", "Series.Code")] <- NULL
income_group_alt[, c("SpecialNotes", "TableName", "X")] <- NULL

# Clean column names: drop "X" year prefixes, rename PWT country code
names(gdp_alt) <- gsub("X", "", names(gdp_alt))
names(cpi_alt) <- gsub("X", "", names(cpi_alt))
names(employment_rate_alt) <- gsub("X", "", names(employment_rate_alt))
names(population_alt) <- gsub("X", "", names(population_alt))
names(trade_openness_alt) <- gsub("X", "", names(trade_openness_alt))
names(agr_share_alt) <- gsub("X", "", names(agr_share_alt))
names(edu_alt) <- gsub("X", "", names(edu_alt))

gdp_alt <- gdp_alt %>%
  gather("year", "gdp", -c(Country.Code, Country.Name))

cpi_alt <- cpi_alt %>%
  gather("year", "cpi", -c(Country.Code, Country.Name))

employment_rate_alt <- employment_rate_alt %>%
  gather("year", "employment", -c(Country.Code, Country.Name))

population_alt <- population_alt %>%
  gather("year", "population", -c(Country.Code, Country.Name))

trade_openness_alt <- trade_openness_alt %>%
  gather("year", "trade_openness", -c(Country.Code, Country.Name))

agr_share_alt <- agr_share_alt %>%
  gather("year", "agr_share", -c(Country.Code, Country.Name))

edu_alt <- edu_alt %>%
  gather("year", "edu", -c(Country.Code, Country.Name))

capital_alt$year <- as.character(capital_alt$year)

edu_alt <- edu_alt %>%
  group_by(Country.Code) %>%
  fill(edu, .direction = "downup") %>%
  ungroup()

# 3.2 Merge macro data (WB-only) ----------------------------------------

merged_alt <- gdp_alt %>%
  full_join(cpi_alt,            by = c("Country.Code", "Country.Name", "year")) %>%
  full_join(employment_rate_alt, by = c("Country.Code", "Country.Name", "year")) %>%
  full_join(population_alt,     by = c("Country.Code", "Country.Name", "year")) %>%
  full_join(trade_openness_alt, by = c("Country.Code", "Country.Name", "year")) %>%
  full_join(agr_share_alt,      by = c("Country.Code", "Country.Name", "year")) %>%
  full_join(edu_alt,            by = c("Country.Code", "Country.Name", "year")) %>%
  full_join(income_group_alt,   by = "Country.Code") %>%
  full_join(
    capital_alt %>% select(Country.Name, Country.Code, capital, year),
    by = c("Country.Code", "Country.Name", "year")
  )

merged_alt$Region      <- factor(merged_alt$Region)
merged_alt$year        <- as.integer(merged_alt$year)

all_countries_alt <- merged_alt %>%
  filter(Region != "",
         year >= 1960,
         year <  2020) %>%
  group_by(Country.Code) %>%
  ungroup()

# 3.3 FE-based TFP using GDP/pop & capital/pop ---------------------------

all_countries_alt <- all_countries_alt %>%
  mutate(
    gdp_log = log(gdp / population),
    k_log   = log(capital / population),
    employment_log = log(employment),
    pop_log        = log(population),
    cpi_log        = log(cpi),
    trade_openness_log = log(trade_openness),
    agr_share_log  = log(agr_share),
    edu_log        = log(edu)
  )

# LLMIC sample only
all_countries_alt <- all_countries_alt %>%
  filter(!IncomeGroup %in% c("High income", "Upper middle income"))

pdata_alt <- pdata.frame(all_countries_alt, index = c("Country.Name", "year"))

panel_model_fixed_alt <- plm(
  gdp_log ~ k_log + employment_log,
  data  = pdata_alt,
  model = "within"
)

summary(panel_model_fixed_alt)
coeftest(panel_model_fixed_alt,
         vcovHC(panel_model_fixed_alt, type = "HC1", cluster = "group"))

# Country FE and residuals → TFP_alt
intercepts_alt <- as.data.frame(summary(fixef(panel_model_fixed_alt))) %>%
  rownames_to_column(var = "Country.Name")

all_countries_alt <- intercepts_alt %>%
  select(Country.Name, Estimate) %>%
  left_join(all_countries_alt, by = "Country.Name")

colnames(all_countries_alt)[2] <- "a_log"

res_alt <- as.data.frame(resid(panel_model_fixed_alt))
colnames(res_alt)[1] <- "residuals_log"

res_alt <- res_alt %>%
  rownames_to_column("Country.Name-year") %>%
  separate("Country.Name-year", into = c("Country.Name", "year"), sep = "-") %>%
  mutate(year = as.integer(year))

all_countries_alt <- res_alt %>%
  full_join(all_countries_alt, by = c("Country.Name", "year")) %>%
  mutate(
    tfp_log = a_log + residuals_log
  )

# Sickle cell (country-level HbS prevalence)
sickle_cell_path <- data_path("Sickle_cell","sickle.csv")
sickle_cell <- read.csv(sickle_cell_path, header = TRUE)

all_countries_alt <- sickle_cell %>%
  inner_join(all_countries_alt, by = "Country.Code")


# WHO malaria incidence

malaria_incidence_path <- data_path("World_malaria_report_2023_Annex","malaria_incidence.csv")
malaria_incidence <- read.csv(malaria_incidence_path, header = TRUE)
# Merge malaria into macro panel
all_countries_alt <- malaria_incidence %>%
  select(
    Country.Name, year,
    case_incidence, cases_log, lag_case_log,
    death_incidence, case_incidence_log, lag_case_incidence_log,
    death_incidence_log
  ) %>%
  inner_join(all_countries_alt, by = c("Country.Name", "year"))

# ENSO

path_olr <- data_path("ENSO", "olr_annual.csv")
olr_annual <- read.csv(path_olr, header = TRUE)

shift_share_iv_alt <- all_countries_alt %>%
  left_join(olr_annual, by = "year") %>%
  mutate(
    hbs_olr_o     = log(results_weighted) * log(olr_o)
  ) %>%
  arrange(Country.Code, year) %>%
  filter(year >= 2000, year < 2020)

tfp_data_alt <- pdata.frame(shift_share_iv_alt, index = c("Country.Name", "year"))


cat("\n==============================\n")
cat("Table A7. Unbalanced Panel Estimates of Malaria Impact on TFP using WDI data \n")
cat("==============================\n\n")

# FE, OLS
tfp_panel_model_fixed_alt <- plm(
  tfp_log ~ case_incidence_log + agr_share_log + trade_openness_log +
    cpi_log + edu_log,
  data  = tfp_data_alt,
  model = "within"
)
summary(tfp_panel_model_fixed_alt)

cat("\n==============================\n")
cat("Table A7. FE \n")
cat("==============================\n\n")
print(coeftest(tfp_panel_model_fixed_alt,
               vcovHC(tfp_panel_model_fixed_alt, type = "HC1", cluster = "group"))
)

# RE
tfp_panel_model_random_alt <- plm(
  tfp_log ~ case_incidence_log + agr_share_log + trade_openness_log +
    cpi_log + edu_log,
  data  = tfp_data_alt,
  model = "random"
)
summary(tfp_panel_model_random_alt)

cat("\n==============================\n")
cat("Table A7. RE \n")
cat("==============================\n\n")
print(coeftest(tfp_panel_model_random_alt,
         vcovHC(tfp_panel_model_random_alt, type = "HC1", cluster = "group")))

#########################
# Lagged-IV
#########################
lagged_iv_model_pooled_alt <- ivreg::ivreg(
  tfp_log ~ case_incidence_log + agr_share_log + trade_openness_log +
    cpi_log + edu_log + factor(Country.Name) |
    lag_case_incidence_log + agr_share_log + trade_openness_log +
    cpi_log + edu_log + factor(Country.Name),
  data = tfp_data_alt
)
summary(lagged_iv_model_pooled_alt, diagnostics = TRUE)

cat("\n==============================\n")
cat("Table A7. 2SLS Lag IV \n")
cat("==============================\n\n")
print(coeftest(lagged_iv_model_pooled_alt,
         vcovHC(lagged_iv_model_pooled_alt, type = "HC1")))

#########################
# Shift–share IV (HbS × OLR)
#########################

shift_share_iv_sickle_alt <- ivreg::ivreg(
  tfp_log ~ case_incidence_log + agr_share_log + trade_openness_log +
    cpi_log + edu_log + factor(Country.Name) |
    hbs_olr_o + agr_share_log + trade_openness_log +
    cpi_log + edu_log + factor(Country.Name),
  data = tfp_data_alt
)
summary(shift_share_iv_sickle_alt, diagnostics = TRUE)

cat("\n==============================\n")
cat("Table A7. 2SLS Shift-share IV \n")
cat("==============================\n\n")
print(coeftest(shift_share_iv_sickle_alt,
         vcovHC(shift_share_iv_sickle_alt, type = "HC1")))

