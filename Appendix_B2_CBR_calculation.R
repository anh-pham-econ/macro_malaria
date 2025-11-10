library(dplyr) 
library(tidyr) 
library(purrr)
library(stringr)

DATA_DIR <- "Data"
OUTPUT_DIR <- "Output"
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, showWarnings = FALSE)
data_path <- function(...) file.path(DATA_DIR, ...)
#Set working directory to the project root
# setwd("path/to/project/root")


tfp_data_path <- data_path("LLMICs_all_data.csv")
tfp_data <- read.csv(tfp_data_path,row.names = NULL)

sample_list_path <- data_path("sample_list.csv")
country_list <- read.csv(sample_list_path, row.names = NULL)

tfp_data_keep <- tfp_data %>%
  semi_join(country_list %>% select(Country.Code), by = "Country.Code")


# ===============================================================
# Malaria CBA — GTS (2016–2030), units in million 2017 USD
# ===============================================================


stopifnot(exists("tfp_data_keep"))
required_cols <- c("Country.Code","Country.Name","year","rgdpna","pop")
missing_cols <- setdiff(required_cols, names(tfp_data_keep))
if (length(missing_cols) > 0) stop("Missing columns in tfp_data_keep: ", paste(missing_cols, collapse=", "))


compute_baseline <- function(df, baseline_years = 2015:2019) {
  df %>%
    mutate(in_window = year %in% baseline_years) %>%
    group_by(Country.Code, Country.Name) %>%
    summarize(
      Y0   = if (any(in_window, na.rm = TRUE)) mean(rgdpna[in_window], na.rm = TRUE)
      else { yrs <- sort(unique(year[!is.na(rgdpna)]), TRUE); mean(rgdpna[year %in% head(yrs, 5)], na.rm = TRUE) },
      Pop0 = if (any(in_window, na.rm = TRUE)) mean(pop[in_window], na.rm = TRUE)
      else { yrs <- sort(unique(year[!is.na(pop)]), TRUE); mean(pop[year %in% head(yrs, 5)], na.rm = TRUE) },
      .groups = "drop"
    ) %>%
    filter(is.finite(Y0), is.finite(Pop0))
}

delta_lnA <- function(theta, r) { stopifnot(r > 0, r < 1); theta * log(1 - r) }

benefit_path <- function(country_tbl, theta, r, k = 1.0,
                         T = 15, path = c("linear","step"), d = 0.05) {
  path <- match.arg(path)
  dlnA <- delta_lnA(theta, r)
  ben_i_full <- country_tbl %>%
    mutate(delta_lnA = dlnA,
           pct_TFP   = 100 * delta_lnA,
           dY_full   = Y0 * delta_lnA * k)   # NOTE: Y0 is in MILLIONS -> dY_full in MILLIONS
  
  mult <- if (path == "step") rep(1, T) else (1:T) / T
  ann  <- tibble(year_index = 1:T, multiplier = mult)
  
  ben_t <- ben_i_full %>%
    tidyr::crossing(ann) %>%
    mutate(dY_it = dY_full * multiplier,
           disc = 1 / (1 + d) ^ year_index,
           PV_it = dY_it * disc)
  
  list(
    per_country = ben_i_full,
    by_year     = ben_t %>% group_by(year_index) %>% summarize(DeltaY = sum(dY_it), PV = sum(PV_it), .groups="drop"),
    NPV_benefit = sum(ben_t$PV_it)
  )
}

cost_from_vector <- function(C_t_millions, d = 0.05) {
  # C_t_millions is already in **million** 2017 USD
  stopifnot(is.numeric(C_t_millions), length(C_t_millions) >= 1)
  tibble(year_index = seq_along(C_t_millions),
         Cost = as.numeric(C_t_millions)) %>%
    mutate(PV = Cost / (1 + d) ^ year_index)
}

run_cba <- function(df_panel,
                    theta = -0.1156, r = 0.9, k = 1.0,
                    T = 15, path = c("linear","step"), d = 0.05,
                    baseline_years = 2015:2019,
                    cost_vector_millions = NULL) {
  
  path <- match.arg(path)
  countries <- compute_baseline(df_panel, baseline_years = baseline_years)
  
  bens <- benefit_path(countries, theta = theta, r = r, k = k, T = T, path = path, d = d)
  
  total_full_effect <- sum(bens$per_country$dY_full, na.rm = TRUE)  # Σ_i ΔY_i (MILLIONS)
  
  if (is.null(cost_vector_millions)) {
    costs <- tibble(year_index = 1:T, Cost = 0, PV = 0)
  } else {
    stopifnot(length(cost_vector_millions) == T)
    costs <- cost_from_vector(cost_vector_millions, d = d)
  }
  
  by_year <- bens$by_year %>%
    left_join(costs, by = "year_index") %>%
    mutate(Net = DeltaY - Cost,
           PV_Net = PV.x - PV.y) %>%
    select(year_index, DeltaY, Cost, Net, PV_benefit = PV.x, PV_cost = PV.y, PV_Net)
  
  out <- list(
    units = "million 2017 USD",
    inputs = list(theta = theta, r = r, k = k, T = T, path = path, d = d, baseline_years = baseline_years),
    per_country_full_effect = bens$per_country %>% arrange(desc(dY_full)),
    annual = by_year,
    NPV_benefit = bens$NPV_benefit,           # MILLIONS
    NPV_cost    = sum(costs$PV),              # MILLIONS
    BCR         = if (sum(costs$PV) > 0) bens$NPV_benefit / sum(costs$PV) else Inf,
    NPV_net     = bens$NPV_benefit - sum(costs$PV),
    DeltaY_full_effect_total = total_full_effect # MILLIONS
  )
  class(out) <- c("cba_result", class(out))
  out
}

print.cba_result <- function(x, ...) {
  cat("---- CBA Summary (units: ", x$units, ") ----\n", sep = "")
  cat("theta =", x$inputs$theta, "| r =", x$inputs$r, "| k =", x$inputs$k,
      "| T =", x$inputs$T, "| path =", x$inputs$path, "| d =", x$inputs$d, "\n")
  cat("Full-effect annual ΔY = Σ_i ΔY_i:", format(round(x$DeltaY_full_effect_total, 3), big.mark=","), x$units, "\n")
  cat("NPV Benefits:", format(round(x$NPV_benefit, 3), big.mark=","), x$units, "\n")
  cat("NPV Costs   :", format(round(x$NPV_cost,    3), big.mark=","), x$units, "\n")
  cat("BCR:", if (is.finite(x$BCR)) round(x$BCR, 3) else "NA", "\n")
  cat("NPV Net:", format(round(x$NPV_net, 3), big.mark=","), x$units, "\n")
  invisible(x)
}

# WHO GTS cost path in **million** USD:
# 7,700M; 8,700M. Start at 0 in 2016 and ramp.
build_gts_cost_millions <- function(years = 2016:2030,
                                    target_2025_billion = 7.7,
                                    target_2030_billion = 8.7,
                                    start_level_million = 0) {
  target_2025_m <- target_2025_billion * 1e3
  target_2030_m <- target_2030_billion * 1e3
  y0 <- min(years)
  cost_m <- sapply(years, function(y) {
    if (y <= 2025) {
      frac <- (y - y0) / (2025 - y0)
      start_level_million + frac * (target_2025_m - start_level_million)
    } else {
      frac <- (y - 2025) / (2030 - 2025)
      target_2025_m + frac * (target_2030_m - target_2025_m)
    }
  })
  tibble(year = years, cost_million = as.numeric(cost_m))
}

apportion_cost <- function(cost_tbl, share = 1.0) {
  cost_tbl %>% mutate(cost_million = cost_million * share)
}

# ---------- RUN ----------
panel <- tfp_data_keep

theta_pref  <- -0.1156
r_target    <- 0.9
k_amp       <- 1.00
disc_rate   <- 0.05
horizon_yrs <- 2016:2030
T_horizon   <- length(horizon_yrs)

gts_cost_tbl_m <- build_gts_cost_millions(years = horizon_yrs,
                                          target_2025_billion = 7.7,
                                          target_2030_billion = 8.7,
                                          start_level_million = 0)
cost_share <- 1.0           # set <1 if your panel is a subset of global burden
gts_cost_tbl_m <- apportion_cost(gts_cost_tbl_m, share = cost_share)
cost_vector_m <- gts_cost_tbl_m$cost_million  # in MILLIONS

res_gts <- run_cba(
  df_panel = panel,
  theta = theta_pref,
  r = r_target,
  k = k_amp,
  T = T_horizon,
  path = "linear",
  d = disc_rate,
  baseline_years = 2015:2019,
  cost_vector_millions = cost_vector_m
)

print(res_gts)

# Sanity check: in linear rollout, last-year ΔY equals Σ_i ΔY_i (both in MILLIONS)
last_year_benefit_m <- res_gts$annual %>% filter(year_index == T_horizon) %>% pull(DeltaY)
cat("\nSanity check — final-year ΔY vs. Σ_i ΔY_i (MILLIONS):\n")
cat("Final-year ΔY: ", round(last_year_benefit_m, 3), "\n",
    "Σ_i ΔY_i:      ", round(res_gts$DeltaY_full_effect_total, 3), "\n", sep = "")

# Annual table with calendar years
annual_tbl <- res_gts$annual %>%
  mutate(calendar_year = min(horizon_yrs) - 1 + year_index) %>%
  select(calendar_year, DeltaY, Cost, Net, PV_benefit, PV_cost, PV_Net)

cat("\n\nAnnual aggregate path (units: million 2017 USD):\n")
print(annual_tbl, row.names = FALSE)

# Implied TFP % at full effect
dlnA_full <- delta_lnA(theta_pref, r_target)
cat(sprintf("\nImplied ΔlnA at full effect: %.3f  (~%0.1f%% TFP)\n", dlnA_full, 100*dlnA_full))


# ===============================================================
# A) Scenario runs
# ===============================================================
scenarios <- tibble(
  name   = c("Conservative", "Average", "High-impact"),
  theta  = c(-0.038, -0.088, -0.116),
  r  = c(0.3, 0.5, 0.9)
)

results <- scenarios %>%
  mutate(
    res = map2(theta, r, ~ run_cba(
      df_panel = panel,
      theta = .x,
      r = .y,
      k = k_amp,
      T = T_horizon,
      path = "linear",
      d = disc_rate,
      baseline_years = 2015:2019,
      cost_vector_millions = cost_vector_m
    ))
  )

# Print summaries
for (i in seq_len(nrow(results))) {
  cat("\n==============================\n")
  cat("Scenario:", results$name[i], "\n")
  print(results$res[[i]])
}

# ===============================================================
# B) Confidence interval for theta
# ===============================================================
theta_hat <- -0.1156
se_theta  <- 0.0416
z_95 <- qnorm(0.975)

ci_lower <- theta_hat - z_95 * se_theta
ci_upper <- theta_hat + z_95 * se_theta

cat("\n95% CI for theta:", round(ci_lower, 3), "to", round(ci_upper, 3), "\n")

# Run CBA at lower and upper bound
res_ci_lower <- run_cba(
  df_panel = panel,
  theta = ci_lower,
  r = r_target,
  k = k_amp,
  T = T_horizon,
  path = "linear",
  d = disc_rate,
  baseline_years = 2015:2019,
  cost_vector_millions = cost_vector_m
)

res_ci_upper <- run_cba(
  df_panel = panel,
  theta = ci_upper,
  r = r_target,
  k = k_amp,
  T = T_horizon,
  path = "linear",
  d = disc_rate,
  baseline_years = 2015:2019,
  cost_vector_millions = cost_vector_m
)

cat("\n--- CI Results ---\n")
cat("Lower-bound theta =", round(ci_lower, 3), "\n")
print(res_ci_lower)

cat("\nUpper-bound theta =", round(ci_upper, 3), "\n")
print(res_ci_upper)

# Report the implied ΔlnA band
dlnA_lower <- delta_lnA(ci_lower, r_target)
dlnA_upper <- delta_lnA(ci_upper, r_target)
cat(sprintf("\nImplied ΔlnA (TFP change) 95%% CI: %.3f (%.1f%%) to %.3f (%.1f%%)\n",
            dlnA_lower, 100*dlnA_lower, dlnA_upper, 100*dlnA_upper))


res_dis_count_3 <- run_cba(
  df_panel = panel,
  theta = theta_pref,
  r = r_target,
  k = k_amp,
  T = T_horizon,
  path = "linear",
  d = 0.05,
  baseline_years = 2015:2019,
  cost_vector_millions = cost_vector_m
)

print(res_dis_count_3)

res_dis_count_10 <- run_cba(
  df_panel = panel,
  theta = theta_pref,
  r = r_target,
  k = k_amp,
  T = T_horizon,
  path = "linear",
  d = 0.1,
  baseline_years = 2015:2019,
  cost_vector_millions = cost_vector_m
)

print(res_dis_count_10)
