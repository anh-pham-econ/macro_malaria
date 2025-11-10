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

hci_file_path <- data_path("HDI","0e94ab23-cc2d-41ca-8cdb-262314eb111a_Data.csv")
hci <- read.csv(hci_file_path, header = TRUE)
hci <- hci[1:217, ]
names(hci) <- gsub("X", "", names(hci))
names(hci) <- gsub(".*(19[0-9]{2}|20[0-9]{2}).*", "\\1", names(hci))
hci[, c('Series.Name','Series.Code')] <- NULL
hci <- gather(hci, key = "year", value = "hci", -c(Country.Code,Country.Name))

hci <- hci %>%
  mutate(
    hci = na_if(hci, ".."),
    hci = as.numeric(hci)
  )

hci <- hci %>%
  group_by(Country.Code) %>%
  fill(hci, .direction = "downup")

hci$year <- as.integer(hci$year)

hci_analysis <- hci %>%
  inner_join(shift_share_iv, c("Country.Code","Country.Name", "year"))

tfp_data_hci <- pdata.frame(hci_analysis, index = c("Country.Code", "year"))

tfp_data_hci <- tfp_data_hci %>%
  filter(!is.na(hci))


cat("\n==============================\n")
cat("Table A13. Causal Mediation Analysis \n")
cat("==============================\n\n")

# Mediator model: how malaria affects HCI
model.m <- lm(hci ~ case_incidence_log + agr_share_log+ trade_openness_log + cpi_log+ edu_log, data = tfp_data_hci)
# Outcome model: how HCI and malaria jointly affect TFP
model.y <- lm(tfp_log ~ case_incidence_log + hci + agr_share_log + trade_openness_log + cpi_log+ edu_log, data = tfp_data_hci)
# Run mediation with 1000 simulations (nonparametric bootstrap)
set.seed(111)
med.out <- mediate(model.m, model.y,
                   treat = "case_incidence_log",
                   mediator = "hci",
                   robustSE = TRUE,
                   sims = 1000)


cat("\n==============================\n")
cat("Table A13. Mediator: HCI \n")
cat("==============================\n\n")

print(summary(med.out))





