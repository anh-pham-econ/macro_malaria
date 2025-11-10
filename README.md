# macro_malaria
Examining the impact of malaria case incidence on total factor productivity leveraging a panel dataset covering 30 malaria-endemic countries from 2000 to 2019

# Replication Package for Estimating the impact of Malaria on Total Factor Productivity

This repository contains the replication code for the malaria–TFP analysis.  
All core results, robustness checks, and mediation analyses can be reproduced from a single master script.

---

## 1. Software Requirements

- **R version:** 4.3.2  
- **IDE:** RStudio

### 1.1. Required R Packages

The main replication scripts use:

```r
library(tidyverse)
library(readxl)
library(plm)
library(ivreg)
library(lmtest)
library(sandwich)
library(data.table)
library(mediation)
library(broom)
```

Additional scripts use:

- **Sickle cell extraction (`sickle_extraction.R`):**
  ```r
  library(tidyverse)
  library(sf)
  library(geodata)
  library(exactextractr)
  library(terra)
  library(raster)
  ```

- **Converting GFCF to capital (`gfcf_to_capital.R`):**
  ```r
  library(dplyr)
  library(tidyr)
  library(zoo)
  library(lmtest)
  library(whitestrap)
  ```

- **Cost–Benefit Analysis (`Appendix_B2_CBR_calculation.R`):**
  ```r
  library(dplyr) 
  library(tidyr) 
  library(purrr)
  library(stringr)
  ```

- **Figures (`macro_malaria_figures.R`):**
  ```r
  library(dplyr)
  library(ggplot2)
  library(gghighlight)
  library(raster)
  library(viridis)
  library(maps)
  ```

You can install any missing packages, for example:

```r
install.packages(c(
  "tidyverse", "readxl", "plm", "ivreg", "lmtest", "sandwich",
  "data.table", "mediation", "broom", "sf", "geodata",
  "exactextractr", "terra", "raster", "tidyr", "zoo",
  "whitestrap", "ggplot2", "gghighlight", "viridis", "maps"
))
```

---

## 2. Directory Structure

The scripts assume the following basic structure:

- `Data/` – raw and processed data files  
- `Output/` – tables, logs, and intermediate outputs  
- `Figure/` – figures  
- `sickle_extraction.R`  
- `gfcf_to_capital.R`  
- `master_macro_malaria.R` – master script 
- `data_processing.R`  
- `main_analysis.R`  
- `A7_robustness_WDI_data_analysis.R`  
- `A8_A15_robustness_output_per_labour_auxiliary_check.R`  
- `A9_robustness_dengue.R`  
- `A10_robustness_HICs.R`  
- `A11_robustness_two_sample_MR.R`  
- `A12_robustness_institutional_quality.R`  
- `A13_causal_mediation_hci.R`  
- `A13_causal_mediation_life_expectancy.R`  
- `A14_auxiliary_check_age_group.R`  
- `Appendix_B2_CBR_calculation.R`
- `macro_malaria_figures.R`
---

## 3. How to Run the Replication

1. **Open R / RStudio.**
2. **Set the working directory** to the project root (where `master_macro_malaria.R` is located), for example:
   ```r
   setwd("path/to/project/root")
   ```
3. **Run the master script:**
   ```r
   source("master_macro_malaria.R")
   ```
   or open `master_macro_malaria.R` in RStudio and run it line by line / as a whole.

The master script calls all core scripts in the correct order and reproduces the main results and robustness checks.

---

## 4. Script Order and Purpose

The master script executes the following scripts in order:

1. **`data_processing.R`**  
   - Cleans and merges raw data, constructs panel variables, and prepares the main analysis dataset.

2. **`main_analysis.R`**  
   - Estimates the baseline malaria–TFP relationship using panel regressions and instrumental variable methods.  
   - Produces main tables and summary statistics used in the paper.

3. **`A7_robustness_WDI_data_analysis.R`**  
   - Replicates the main results using alternative measures or series from the World Development Indicators (WDI).

4. **`A8_A15_robustness_output_per_labour_auxiliary_check.R`**  
   - Robustness checks focusing on output per worker / per labour formulations and related auxiliary tests.

5. **`A9_robustness_dengue.R`**  
   - Uses dengue as an alternative or additional disease exposure to test the specificity of the malaria–TFP relationship.

6. **`A10_robustness_HICs.R`**  
   - Robustness checks restricting or extending the sample to include / exclude high-income countries.

7. **`A11_robustness_two_sample_MR.R`**  
   - Two-sample Mendelian randomization (MR) analysis to assess causality using genetic instruments.

8. **`A12_robustness_institutional_quality.R`**  
   - Robustness checks controlling for or interacting with institutional quality measures.

9. **`A13_causal_mediation_hci.R`**  
   - Causal mediation analysis using the Human Capital Index (HCI) as a mediator between malaria and TFP.

10. **`A13_causal_mediation_life_expectancy.R`**  
   - Causal mediation analysis using life expectancy as a mediator.

11. **`A14_auxiliary_check_age_group.R`**  
   - Auxiliary checks by age group (e.g. focusing on child vs adult malaria exposure) to explore heterogeneity in effects.
     
Each R file (e.g. data_processing.R, main_analysis.R, A7_robustness_WDI_data_analysis.R,...) can also be executed independently outside the master script.
To do so, you need to uncomment and edit the header in each script to load required packages, set up directories, and define paths. This ensures that packages, file paths and dependencies are correctly initialized when running any script individually.
---

## 5. Additional (Standalone) Scripts

These scripts are not strictly required to run the main replication but are included for completeness and to reproduce specific inputs or figures:

- **`sickle_extraction.R`**  
  - Extracts sickle cell (HbS) allele frequency data from raster / geospatial layers and aggregates it to the units used in the analysis (e.g. country-year or similar).  
  - Requires spatial and raster libraries: `sf`, `geodata`, `exactextractr`, `terra`, `raster`, plus `tidyverse`.

- **`gfcf_to_capital.R`**  
  - Converts gross fixed capital formation (GFCF) time series into a capital stock series using a perpetual inventory or related method.  
  - Uses `dplyr`, `tidyr`, `zoo`, `lmtest`, and `whitestrap`.

- **`Appendix_B2_CBR_calculation.R`**  
  - Illustrative cost-benefit calculation linking reductions in malaria incidence to potential productivity and income gains
  - Uses `dplyr`, `tidyr`, `purrr`, and `stringr`.

- **`macro_malaria_figures.R`**  
  - Generates maps and other figures used in the paper (e.g. malaria prevalence maps, HBS maps, TFP patterns).  
  - Relies on `ggplot2`, `gghighlight`, `raster`, `viridis`, `maps`, and `dplyr`.

To run any of these, set the working directory to the project root and then:

```r
source("sickle_extraction.R")
# or
source("gfcf_to_capital.R")
# or
source("macro_malaria_figures.R")
```

---

## 6. Output

- Tables, and intermediate results are written to the `Output/` directory (or as specified in each script).
- Figures are written to the `Figure/` directory (or as specified in each script).
- Filenames correspond closely to the table and figure numbers in the paper and appendix, facilitating cross-referencing.

---

## 7. Contact

For questions about the replication package or issues running the code, please contact the author of the study.
