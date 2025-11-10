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


cat("\n==============================\n")
cat("Table A14. Auxiliary consistency checks \n")
cat("==============================\n\n")


###Children under 5 interaction
under_5_children_file_path <- data_path("children-under-age-5","children-under-age-5.csv")
under_5_children_data <- read.csv(under_5_children_file_path, header = TRUE)


population <- read.csv(
  data_path("population", "API_SP.POP.1564.TO_DS2_en_csv_v2_169.csv"),
  header = TRUE, skip = 4
)
columns_to_remove <- c("Indicator.Code", "Indicator.Name", "X")
population[, columns_to_remove]     <- NULL
names(population) <- gsub("X", "", names(population))
population <- population %>%
  gather("year", "population", -c(Country.Code, Country.Name))
population$year <- as.integer(population$year)

names(under_5_children_data)[2] <- "Country.Code"
names(under_5_children_data)[3] <- "year"
names(under_5_children_data)[4] <- "under_5_children_number"
under_5_children_data[, 5] <- NULL
under_5_children_data[, 1] <- NULL

under_5_children_data$year <- as.integer(under_5_children_data$year)
under_5_children_data$under_5_children_number_log = log(under_5_children_data$under_5_children_number)

under_5_children_interaction_analysis <- under_5_children_data %>%
  inner_join(shift_share_iv, by = c("Country.Code", "year"))


under_5_children_interaction_analysis <- under_5_children_interaction_analysis %>%
  full_join(population,     by = c("Country.Code", "Country.Name", "year")) 
  

under_5_children_interaction_analysis$under_5_share <- under_5_children_interaction_analysis$under_5_children_number/under_5_children_interaction_analysis$population
under_5_children_interaction_analysis$under_5_share_log <- log(under_5_children_interaction_analysis$under_5_share)

under_5_children_interaction_analysis$under_5_child_share_interaction <- under_5_children_interaction_analysis$case_incidence_log*under_5_children_interaction_analysis$under_5_share_log

tfp_data_under_5_children_interaction_analysis <- pdata.frame(under_5_children_interaction_analysis, index = c("Country.Code", "year"))


under_5_children_interaction_model <- plm(
  tfp_log ~ case_incidence_log + under_5_child_share_interaction + 
    agr_share_log + trade_openness_log + cpi_log + edu_log ,
  data = tfp_data_under_5_children_interaction_analysis,
  index = c("Country.Code", "year"),
  model = "within"
)

cat("\n==============================\n")
cat("Table A14. Share of children under 5 \n")
cat("==============================\n\n")
print(summary(under_5_children_interaction_model))



###Children 0-14
children_file_path <- data_path("population_0_14_total_share","API_SP.POP.0014.TO.ZS_DS2_en_csv_v2_506997.csv")

children_0_14 <- read.csv(children_file_path, header = TRUE, skip = 4)

children_0_14[, columns_to_remove] <- NULL
names(children_0_14) <- gsub("X", "", names(children_0_14))
children_0_14 <- gather(children_0_14, key = "year", value = "children_share", -c(Country.Code,Country.Name))

children_0_14$year <- as.integer(children_0_14$year)
children_0_14$children_share_log = log(children_0_14$children_share)

children_interaction_analysis <- children_0_14 %>%
  left_join(shift_share_iv, by = c("Country.Code","Country.Name", "year"))

children_interaction_analysis$child_interaction <- children_interaction_analysis$case_incidence_log*children_interaction_analysis$children_share

setDT(children_interaction_analysis)[, lag_child_interaction := shift(child_interaction, n = 5), by = Country.Name]

tfp_data_child_interaction_analysis <- pdata.frame(children_interaction_analysis, index = c("Country.Code", "year"))


child_interaction_model <- plm(
  tfp_log ~ case_incidence_log + lag_child_interaction +
    agr_share_log + trade_openness_log + cpi_log + edu_log ,
  data = tfp_data_child_interaction_analysis,
  index = c("Country.Code", "year"),
  model = "within"
)

cat("\n==============================\n")
cat("Table A14. Share of children from 0 to 14 years old \n")
cat("==============================\n\n")
print(summary(child_interaction_model))

###Edu interaction
primary_file_path <- data_path("Edu/primary","API_SE.PRM.CMPT.ZS_DS2_en_csv_v2_23374.csv")

primary <- read.csv(primary_file_path, header = TRUE, skip = 4)

primary[, c('Indicator.Name', 'Indicator.Code', 'X')] <- NULL
names(primary) <- gsub("X", "", names(primary))
primary <- gather(primary, key = "year", value = "primary_completion", -c(Country.Code,Country.Name))
primary$year <- as.integer(primary$year)
primary <- primary %>%
  group_by(Country.Code) %>%
  fill(primary_completion, .direction = "downup")
primary$primary_completion_log <- log(primary$primary_completion)
edu_interaction_analysis <- shift_share_iv %>%
  left_join(primary, by = c("Country.Code", "year"))

edu_interaction_analysis$edu_interaction <- edu_interaction_analysis$case_incidence_log*edu_interaction_analysis$primary_completion_log
setDT(edu_interaction_analysis)[, lag_edu_interaction_9 := shift(edu_interaction, n =9), by = Country.Code]

tfp_data_edu_interaction_analysis <- pdata.frame(edu_interaction_analysis, index = c("Country.Code", "year"))

edu_interaction_model_lagged_9 <- plm(
  tfp_log ~ case_incidence_log + lag_edu_interaction_9 +
    agr_share_log + trade_openness_log + cpi_log + primary_completion_log+  edu_log,
  data = tfp_data_edu_interaction_analysis,
  index = c("Country.Code", "year"),
  model = "within"
)


cat("\n==============================\n")
cat("Table A14. Primary school completion rate \n")
cat("==============================\n\n")

print(summary(edu_interaction_model_lagged_9))










