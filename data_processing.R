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

#########################
# Load core macro data 
#########################

pwt_file_path <- data_path(
  "PWT",
  "pwt1001.xlsx"
)
pwt <- read_excel(pwt_file_path, sheet = "Data")
pwt$year <- as.character(pwt$year)

cpi_file_path <- data_path(
  "GDP_deflator",
  "API_FP.CPI.TOTL.ZG_DS2_en_csv_v2_77.csv"
)
cpi <- read.csv(cpi_file_path, header = TRUE, skip = 4)

income_group_file_path <- data_path(
  "Income_group",
  "Metadata_Country_API_NE.GDI.FTOT.CD_DS2_en_csv_v2_6303108.csv"
)
income_group <- read.csv(income_group_file_path, header = TRUE)

trade_file_path <- data_path(
  "trade_openess",
  "API_NE.TRD.GNFS.ZS_DS2_en_csv_v2_466.csv"
)
trade_openness <- read.csv(trade_file_path, header = TRUE, skip = 4)

agr_file_path <- data_path(
  "Agriculture_constant_2015_USD",
  "agr_share.csv"
)
agr_share <- read.csv(agr_file_path, header = TRUE, skip = 4)

edu_file_path <- data_path(
  "Edu",
  "edu_series_completion_upper_secondary.csv"
)
edu <- read.csv(edu_file_path, header = TRUE)

#########################
# Tidy and reshape macro data
#########################

columns_to_remove <- c("Indicator.Code", "Indicator.Name", "X")

cpi[, columns_to_remove] <- NULL
trade_openness[, columns_to_remove] <- NULL
agr_share[, columns_to_remove] <- NULL
income_group[, c("SpecialNotes", "TableName", "X")] <- NULL
edu[, c("Series", "Series.Code")] <- NULL

# Clean column names: drop "X" year prefixes, rename PWT country code
names(cpi) <- gsub("X", "", names(cpi))
names(trade_openness) <- gsub("X", "", names(trade_openness))
names(agr_share) <- gsub("X", "", names(agr_share))
names(edu) <- gsub("X", "", names(edu))
names(pwt) <- gsub("countrycode", "Country.Code", names(pwt))

# Convert wide to long

cpi <- cpi %>%
  gather(key = "year", value = "cpi", -c(Country.Code, Country.Name))

trade_openness <- trade_openness %>%
  gather(key = "year", value = "trade_openness", -c(Country.Code, Country.Name))

agr_share <- agr_share %>%
  gather(key = "year", value = "agr_share", -c(Country.Code, Country.Name))

edu <- edu %>%
  gather(key = "year", value = "edu", -c(Country.Code, Country.Name))

# Fill missing education within country
edu <- edu %>%
  group_by(Country.Code) %>%
  fill(edu, .direction = "downup") %>%
  ungroup()

#########################
# 4. Merge macro data and construct panel
#########################

merged_df <- cpi %>%
  full_join(trade_openness, by = c("Country.Code", "Country.Name", "year")) %>%
  full_join(agr_share, by = c("Country.Code", "Country.Name", "year")) %>%
  full_join(income_group, by = "Country.Code") %>%
  full_join(edu, by = c("Country.Code", "Country.Name", "year")) %>%
  full_join(pwt, by = c("Country.Code", "year"))

merged_df$Region <- factor(merged_df$Region)
merged_df$IncomeGroup <- factor(merged_df$IncomeGroup)
merged_df$year <- as.integer(merged_df$year)

all_countries <- merged_df %>%
  filter(Region != "",
         year >= 1960,
         year < 2020) %>%
  group_by(Country.Code) %>%
  ungroup()

all_countries <- all_countries %>%
  mutate(
    gdp_log            = log(rgdpna),
    k_log              = log(rkna),
    employment_log     = log(emp),
    cpi_log            = log(cpi),
    trade_openness_log = log(trade_openness),
    agr_share_log      = log(agr_share),
    edu_log            = log(edu),
    Y_L = rgdpo / emp
  )

HICs <- all_countries %>%
  filter(IncomeGroup %in% c("High income", "Upper middle income"))

LLMICs <- all_countries %>%
  filter(!IncomeGroup %in% c("High income", "Upper middle income"))

# write.csv(HICs, file = "Data/HICs_macro_data.csv", row.names = FALSE)
# write.csv(LLMICs, file = "Data/LLMICs_macro_data.csv", row.names = FALSE)
# write.csv(all_countries, file = "Data/all_countries_macro_data.csv", row.names = FALSE)

#########################
# Cobb–Douglas TFP from PWT (country FE)
#########################

pdata <- pdata.frame(LLMICs, index = c("Country.Name", "year"))

panel_model_fixed <- plm(gdp_log ~ k_log + employment_log,
                         data = pdata, model = "within")

summary(panel_model_fixed)
production_function_reg_results <- coeftest(panel_model_fixed,
         vcovHC(panel_model_fixed, type = "HC1", cluster = "group"))

cat("\n==============================\n")
cat("Table 1. Robust Coefficient Test (HC1 Clustered)\n")
cat("==============================\n\n")
print(production_function_reg_results)

production_function_reg_table <- data.frame(
  Variable    = rownames(production_function_reg_results),
  Estimate    = production_function_reg_results[, 1],
  Std.Error   = production_function_reg_results[, 2],
  t.value     = production_function_reg_results[, 3],
  p.value     = production_function_reg_results[, 4]
)

#write.csv(production_function_reg_table,
#          file = file.path(OUTPUT_DIR, "Table_1_production_function_fixed_effects.csv"),
#          row.names = FALSE)

# Extract country fixed effects
intercept_value <- as.data.frame(summary(fixef(panel_model_fixed))) %>%
  rownames_to_column(var = "Country.Name")

LLMICs <- intercept_value %>%
  select(Country.Name, Estimate) %>%
  left_join(LLMICs, by = "Country.Name")

colnames(LLMICs)[2] <- "a_log"

# Residuals (time-varying TFP component)
res <- as.data.frame(resid(panel_model_fixed))
colnames(res)[1] <- "residuals_log"

res <- res %>%
  rownames_to_column(var = "Country.Name-year") %>%
  separate("Country.Name-year", into = c("Country.Name", "year"), sep = "-") %>%
  mutate(year = as.integer(year))

LLMICs <- res %>%
  full_join(LLMICs, by = c("Country.Name", "year"))

# TFP level and index (2017 = 1)
LLMICs <- LLMICs %>%
  mutate(
    tfp_log = a_log + residuals_log,
    tfp     = exp(tfp_log)
  ) %>%
  group_by(Country.Name) %>%
  mutate(
    tfp_base  = tfp[year == 2017],
    tfp_index = tfp / tfp_base
  ) %>%
  ungroup()

#########################
# Sickle cell distribution and malaria incidence
#########################

# Sickle cell (country-level HbS prevalence)
sickle_cell_path <- data_path("Sickle_cell","sickle.csv")
sickle_cell <- read.csv(sickle_cell_path, header = TRUE)

LLMICs <- sickle_cell %>%
  inner_join(LLMICs, by = "Country.Code")

# WHO malaria incidence
malaria_incidence_path <- data_path(
  "World_malaria_report_2023_Annex",
  "WMR2023_Annex_4F_all_countries.csv"
)
malaria_incidence <- read.csv(malaria_incidence_path, header = TRUE)

colnames(malaria_incidence)[3] <- "pop_denominator"
malaria_incidence$Country.Name[malaria_incidence$Country.Name == "Gambia"] <- "Gambia, The"
malaria_incidence$Country.Name[malaria_incidence$Country.Name == "Yemen"]  <- "Yemen, Rep."
malaria_incidence$year <- as.integer(malaria_incidence$year)

malaria_incidence <- malaria_incidence %>%
  mutate(
    case_incidence  = Cases  / pop_denominator * 100000,
    death_incidence = Deaths / pop_denominator * 100000,
    case_incidence_log  = log(case_incidence),
    death_incidence_log = log(death_incidence),
    cases_log  = log(Cases),
    deaths_log = log(Deaths)
  ) %>%
  group_by(Country.Name) %>%
  filter(any(!is.infinite(case_incidence_log))) %>%
  ungroup()

# Lags
setDT(malaria_incidence)[, lag_case_incidence_log := shift(case_incidence_log),
                         by = Country.Name]
setDT(malaria_incidence)[, lag_case_log := shift(cases_log),
                         by = Country.Name]

malaria_incidence[malaria_incidence == -Inf] <- NA
malaria_incidence <- na.omit(malaria_incidence)

#write.csv(malaria_incidence, file = "Data/World_malaria_report_2023_Annex/malaria_incidence.csv", row.names = FALSE)
# Merge malaria into macro panel
LLMICs <- malaria_incidence %>%
  select(
    Country.Name, year,
    case_incidence, cases_log, lag_case_log,
    death_incidence, case_incidence_log, lag_case_incidence_log,
    death_incidence_log
  ) %>%
  inner_join(LLMICs, by = c("Country.Name", "year"))

#########################
# 7. ENSO / OLR: construct annual OLR and shift–share instrument
#########################

path_olr <- data_path("ENSO", "olr.txt")
all_lines <- read_lines(path_olr)

extract_block <- function(lines, keyword) {
  start <- grep(keyword, lines)[1] + 3  # skip 3 header lines
  nxt <- grep("^\\s*OUTGOING LONG WAVE RADIATION", lines)
  end <- nxt[nxt > start][1]
  if (is.na(end)) end <- length(lines) + 1
  lines[start:(end - 1)]
}


blk_orig <- extract_block(all_lines, "ORIGINAL")
blk_anom <- extract_block(all_lines, "ANOMALY")
blk_std  <- extract_block(all_lines, "STANDARDIZED")

read_block <- function(block_lines) {
  df <- read_table2(
    paste(block_lines, collapse = "\n"),
    col_names = c(
      "year", "JAN", "FEB", "MAR", "APR", "MAY", "JUN",
      "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"
    ),
    col_types = cols(.default = col_double(), year = col_integer())
  )
  df[df == -999.9] <- NA
  df
}

olr_o <- read_block(blk_orig)
olr_a <- read_block(blk_anom)
olr_s <- read_block(blk_std)

annualize <- function(df, value_name) {
  df %>%
    pivot_longer(-year, names_to = "month", values_to = value_name) %>%
    group_by(year) %>%
    summarise(!!value_name := mean(.data[[value_name]], na.rm = TRUE),
              .groups = "drop")
}

olr_o_y <- annualize(olr_o, "olr_o")
olr_a_y <- annualize(olr_a, "olr_a")
olr_s_y <- annualize(olr_s, "olr_std")

olr_annual <- olr_o_y %>%
  left_join(olr_a_y, by = "year") %>%
  left_join(olr_s_y, by = "year") %>%
  arrange(year)
#write.csv(olr_annual, file = "Data/ENSO/olr_annual.csv", row.names = FALSE)

# Shift–share instruments
shift_share_iv <- LLMICs %>%
  left_join(olr_annual, by = c("year")) %>%
  mutate(
    hbs_olr_o     = log(results_weighted) * log(olr_o),
    ) %>%
  arrange(Country.Code, year) %>%
  filter(year >= 2000, year < 2020)


#write.csv(shift_share_iv, file = "Data/LLMICs_all_data.csv", row.names = FALSE)

### A3 Summary Statistics 

vars <- c("rgdpna",
          "rkna",
          "emp",
          "cpi",
          "trade_openness",
          "agr_share",
          "edu",
          "case_incidence",
          "results_weighted",
          "olr_o",
          "tfp")

id_var  <- "Country.Code"   

panel_summ <- lapply(vars, function(v) {
  x  <- shift_share_iv[[v]]
  id <- shift_share_iv[[id_var]]
  ok    <- !is.na(x) & !is.na(id)
  x_ok  <- x[ok]
  id_ok <- id[ok]
  
  # overall stats
  obs    <- length(x_ok)
  mean_x <- mean(x_ok)
  min_x  <- min(x_ok)
  max_x  <- max(x_ok)
  sd_tot <- sd(x_ok)
  
  # between: SD of country means
  country_means <- tapply(x_ok, id_ok, mean)
  sd_between    <- sd(country_means)
  
  # within: SD of demeaned values (x_it - \bar{x}_i)
  country_mean_vec <- country_means[match(id_ok, names(country_means))]
  x_demeaned       <- x_ok - country_mean_vec
  sd_within        <- sd(x_demeaned)
  
  c(Obs        = obs,
    Mean       = mean_x,
    Min        = min_x,
    Max        = max_x,
    SD         = sd_tot,
    Within_SD  = sd_within,
    Between_SD = sd_between)
})

panel_summ <- do.call(rbind, panel_summ)
rownames(panel_summ) <- vars
panel_summ_rounded <- round(panel_summ, 4)


cat("\n==============================\n")
cat("Table A3. Summary Statistics\n")
cat("==============================\n\n")
print(panel_summ_rounded)

# write.csv(panel_summ_rounded,
#           file = file.path(OUTPUT_DIR, "Table_A3_panel_summary.csv"),
#           row.names = TRUE)

