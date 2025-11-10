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

wgi_file_path <- data_path("wgidataset_with_sourcedata_excel","wgidataset_with_sourcedata.xlsx")
wgi <- read_excel(wgi_file_path,na = c("...", "..", "") )
wgi$year <- as.integer(wgi$year)

wgi<- wgi %>%
  select(Country.Code, year, indicator, estimate)%>%
  filter(year >= '2000') %>%
  filter(year < '2020') %>%
  group_by(Country.Code) %>%
  ungroup()

wgi <- wgi %>% 
  pivot_wider(names_from  = indicator,    # or Indicator Name
              values_from = estimate)
wgi <- wgi %>%
  group_by(Country.Code) %>%
  fill(cc, rl, rq, .direction = "downup")

merged__wgi <- wgi %>%
  select(Country.Code, year, cc, rl, rq) %>%
  left_join(shift_share_iv, by = c("Country.Code","year"))


cat("\n==============================\n")
cat("Table A12. Controlling for Institutional Quality \n")
cat("==============================\n\n")

tfp_data_wgi <- pdata.frame(merged__wgi, index = c("Country.Code", "year"))
clust <- tfp_data_wgi$Country.Name

iv_shiftshare_1 <- ivreg::ivreg(tfp_log ~ case_incidence_log + agr_share_log + trade_openness_log + cpi_log  + edu_log + cc + factor(Country.Name)|
                         (hbs_olr_o)  + agr_share_log + trade_openness_log + cpi_log  + edu_log + cc+ factor(Country.Name),
                       data = tfp_data_wgi)
summary(iv_shiftshare_1)

cat("\n==============================\n")
cat("Table A12. Control of corruption (cc) \n")
cat("==============================\n\n")
print(coeftest(iv_shiftshare_1,vcovHC(iv_shiftshare_1, type = 'HC1', cluster = clust)))

iv_shiftshare_2 <- ivreg::ivreg(tfp_log ~ case_incidence_log + agr_share_log + trade_openness_log + cpi_log  + edu_log + rl + factor(Country.Name)|
                                  (hbs_olr_o)  + agr_share_log + trade_openness_log + cpi_log  + edu_log + rl+ factor(Country.Name),
                                data = tfp_data_wgi)
summary(iv_shiftshare_2)

cat("\n==============================\n")
cat("Table A12. Rule of law (rl) \n")
cat("==============================\n\n")
print(coeftest(iv_shiftshare_2,vcovHC(iv_shiftshare_2, type = 'HC1', cluster = clust)))

iv_shiftshare_3 <- ivreg::ivreg(tfp_log ~ case_incidence_log + agr_share_log + trade_openness_log + cpi_log  + edu_log + rq + factor(Country.Name)|
                                  (hbs_olr_o)  + agr_share_log + trade_openness_log + cpi_log  + edu_log + rq+ factor(Country.Name),
                                data = tfp_data_wgi)
summary(iv_shiftshare_3)

cat("\n==============================\n")
cat("Table A12. Regulatory quality (rq) \n")
cat("==============================\n\n")
print(coeftest(iv_shiftshare_3,vcovHC(iv_shiftshare_3, type = 'HC1', cluster = clust)))

iv_shiftshare_4 <- ivreg::ivreg(tfp_log ~ case_incidence_log + agr_share_log + trade_openness_log + cpi_log  + edu_log + cc + rl + rq + factor(Country.Name)|
                                  (hbs_olr_o)  + agr_share_log + trade_openness_log + cpi_log  + edu_log + cc + rl + rq+ factor(Country.Name),
                                data = tfp_data_wgi)
summary(iv_shiftshare_4)

cat("\n==============================\n")
cat("Table A12. All factors (cc,rl,rq) \n")
cat("==============================\n\n")
print(coeftest(iv_shiftshare_4,vcovHC(iv_shiftshare_4, type = 'HC1', cluster = clust)))



