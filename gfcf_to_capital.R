library(dplyr)
library(tidyr)
library(zoo)
library(lmtest)
library(whitestrap)

DATA_DIR <- "Data"    
OUTPUT_DIR <- "Output"    
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, showWarnings = FALSE)
data_path <- function(...) file.path(DATA_DIR, ...)

## Set working directory to the project root
# setwd("path/to/project/root")

gfcf <- read.csv(
  data_path("Capital_constant_2015", "API_NE.GDI.FTOT.KD_DS2_en_csv_v2_11471.csv"),
  header = TRUE, skip = 4
)


columns_to_remove <- c("Indicator.Code","Indicator.Name",'X')

gfcf[, columns_to_remove] <- NULL
names(gfcf) <- gsub("X", "", names(gfcf))
gfcf <- gather(gfcf, key = "year", value = "gfcf", -c(Country.Code,Country.Name))

gfcf <- gfcf %>%
  arrange(Country.Code, year) %>%
  group_by(Country.Code) %>%
  mutate(gfcf_growth = (gfcf - dplyr::lag(gfcf)) / dplyr::lag(gfcf))

# Calculate the initial capital values and the average growth per country
# Ref source
#https://personal.lse.ac.uk/casellif/papers/mpk.pdf
#page 544

initial_capital <- gfcf %>%
  group_by(Country.Code) %>%
  summarize(
    first_gfcf = dplyr::first(na.omit(gfcf)),
    avg_growth = mean(gfcf_growth, na.rm = TRUE),
    year = year[which(!is.na(gfcf))[1]]
  ) %>%
  filter(!is.na(first_gfcf)) %>%
  mutate(initial_capital = first_gfcf / avg_growth)

# Merge the initial capital back into the main data frame
gfcf <- gfcf %>%
  left_join(initial_capital %>% select(Country.Code, capital = initial_capital, year), by = c("Country.Code","year"))


# Function to fill initial_capital for each country dataset
fill_capital <- function(df) {
  for (i in 2:nrow(df)) {
    if (is.na(df$capital[i])) {
      df$capital[i] <- df$gfcf[i-1] + df$capital[i-1]
    }
  }
  return(df)
}

# Split the dataset by country
data_split <- split(gfcf, list(gfcf$Country.Name, gfcf$Country.Code))

# Apply the function to each subset and combine results
data_filled <- do.call(rbind, lapply(data_split, function(df) {
  if (nrow(df) > 0) {
    fill_capital(df)
  } else {
    df   }
}))

#write.csv(data_filled, file = "Data/Capital_constant_2015/gfcf_to_capital.csv", row.names = FALSE)

cat("\nScripts completed successfully.\n")

