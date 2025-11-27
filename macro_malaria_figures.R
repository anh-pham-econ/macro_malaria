library(dplyr)
library(ggplot2)
library(gghighlight)
library(raster)
library(viridis)
library(maps)

DATA_DIR <- "Data"    
OUTPUT_DIR <- "Output"    
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, showWarnings = FALSE)
data_path <- function(...) file.path(DATA_DIR, ...)

##Set working directory to the project root
# setwd("path/to/project/root")

data_file_path <- data_path("LLMICs_all_data.csv")
data <- read.csv(data_file_path, header = TRUE)
data$year = as.integer(data$year)

##Fig 1: cases_by_country
malaria_data <- data %>%
  filter(year %in% 2001:2019) %>%
  group_by(Country.Name) %>%
  filter(all(2001:2019 %in% year)) %>%
  ungroup()

high_incidence_countries <- malaria_data %>%
  group_by(Country.Name) %>%
  summarise(mean_inc = mean(case_incidence, na.rm = TRUE)) %>%
  filter(mean_inc > 38000)

n_colors <- nrow(high_incidence_countries) + 1
vir_cols <- viridis(n_colors, option = "C")

fig1_malaria_cases <- ggplot(malaria_data,
                             aes(x = year,
                                 y = case_incidence,
                                 group = Country.Name,
                                 color = Country.Name)) +
  geom_line(linewidth = 0.8) +
  gghighlight(mean(case_incidence, na.rm = TRUE) > 38000,
              use_direct_label = FALSE,
              unhighlighted_params = list(alpha = 0.2, linewidth = 0.5)) +
  labs(
    x = "Year",
    y = "Cases per 100,000 population at risk",
    color = NULL
  ) +
  #scale_color_viridis_d(option = "C") + 
  scale_color_manual(values = vir_cols) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom"
  )+ guides(col = guide_legend(nrow = 2))
fig1_malaria_cases


# ggsave(
#   filename = "cases_by_country.png",   
#   plot = fig1_malaria_cases,           
#   path = "Figure",
#   width = 10,                        
#   height = 6,                       
#   dpi = 600                         
# )


##Fig 3: blood_disorder_distribution
sickle_raster_path <- data_path(
  "2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed",
  "2013_Sickle_Haemoglobin_HbS_Allele_Freq_Global_5k_Decompressed.geotiff"
)
sickle <- raster(sickle_raster_path)

sickle_df <- as.data.frame(sickle, xy = TRUE)
colnames(sickle_df) <- c("x", "y", "hbs")

world_df <- map_data("world")

fig3_hbs_map <- ggplot() +
  geom_raster(data = sickle_df,
              aes(x = x, y = y, fill = hbs)) +
  geom_path(data = world_df,
            aes(x = long, y = lat, group = group),
            color = "white",
            linewidth = 0.2) +
  scale_fill_viridis(option = "C",
                     name = "HbS allele\nfrequency",
                     na.value = "white"
  ) + 
  labs(
  ) + coord_equal() +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    axis.line  = element_blank(), 
    legend.position = "right"
  )
fig3_hbs_map
# 
# ggsave(
#   filename = "blood_disorder_distribution.png",  
#   plot = fig3_hbs_map,             
#   path = "Figure",
#   width = 10,                      
#   height = 6,                      
#   dpi = 600                        
# )

##Fig 4: malaria_HbS

corr_2019 <- data %>%      
  dplyr::filter(year == 2019)

fig4_corr_hbs <- ggplot(corr_2019,
                        aes(x = case_incidence,
                            y = results_weighted,
                            label = Country.Code)) +
  geom_point(color = "#002395", size = 2) +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "#C00000",
              linewidth = 0.8) +
  geom_text(vjust = -1,
            size = 3,
            color = "#002395",
            check_overlap = TRUE) +
  labs(
    x = "Malaria cases per 100,000 at risk",
    y = "HbS frequency in the population"
  ) +
  theme_classic(base_size = 12) +
  theme(
  )
fig4_corr_hbs
# 
# ggsave(
#   filename = "malaria_HbS.png",  
#   plot = fig4_corr_hbs,          
#   path = "Figure",
#   width = 10,                    
#   height = 6,                    
#   dpi = 600                      
# )

##Fig 5: time_series_OLR_malaria

yearly_avg <- data %>%
  group_by(year) %>%
  summarise(
    mean_case_incidence_log = mean(case_incidence_log, na.rm = TRUE),
    mean_log_olr_o          = mean(log(olr_o), na.rm = TRUE),
    .groups = "drop"
  )

fig5_time_series <- ggplot(yearly_avg, aes(x = year)) +
  geom_line(aes(y = mean_case_incidence_log,
                linetype = "Malaria incidence (log)",
                color    = "Malaria incidence (log)"),color = "#C00000",
            linewidth = 0.8) +
  geom_line(aes(y = mean_log_olr_o,
                linetype = "OLR (log)",
                color    = "OLR (log)"), color = "#002395",
            linewidth = 0.8) +
  labs(
    x = "Year",
    y = "Log value (cross-country mean)",
    color = NULL,
    linetype = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom"
  )
fig5_time_series
# 
# ggsave(
#   filename = "time_series_OLR_malaria.png", 
#   plot = fig5_time_series,            
#   path = "Figure",
#   width = 10,                       
#   height = 6,                     
#   dpi = 600                       
# )

##Fig 6: Correlation_malaria_OLR

fig6_corr_olr <- ggplot(yearly_avg,
                        aes(x = mean_case_incidence_log,
                            y = mean_log_olr_o)) +
  geom_point(color = "#002395", size = 2) +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "#C00000",
              linewidth = 0.8) +
  labs(
    x = "Mean malaria incidence (log)",
    y = "Mean OLR (log)"
  ) +
  theme_classic(base_size = 12) +
  theme(
  )
fig6_corr_olr
# 
# ggsave(
#   filename = "Correlation_malaria_OLR.png",  
#   plot = fig6_corr_olr,           
#   path = "Figure",
#   width = 10,                      
#   height = 6,                      
#   dpi = 600                        
# )

##Fig 7: added-variable-plot
reg_y <- lm(case_incidence_log ~ agr_share_log + trade_openness_log +
              cpi_log + edu_log + factor(Country.Name),
            data = data)

reg_x <- lm(hbs_olr_o ~ agr_share_log + trade_openness_log +
              cpi_log + edu_log + factor(Country.Name),
            data = data)

av_df <- data.frame(
  hbs_res  = resid(reg_x),
  case_res = resid(reg_y)
)

fig7_avplot <- ggplot(av_df,
                      aes(x = hbs_res,
                          y = case_res)) +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             linewidth = 0.4) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             linewidth = 0.4) +
  geom_point(alpha = 0.5, size = 2,color = "#002395") +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "#C00000",
              linewidth = 0.8) +
  labs(
    x = "Shift–share IV residual",
    y = "Malaria incidence residual"
  ) +
  theme_classic(base_size = 12) +
  theme()
fig7_avplot

# ggsave(
#   filename = "added-variable-plot.png",  
#   plot = fig7_avplot,            
#   path = "Figure",
#   width = 10,                    
#   height = 6,                    
#   dpi = 600                      
# )

##Fig 8: correlation_malaria_shift_share_by_country

fig8_corr_iv <- ggplot(data,
                       aes(x = case_incidence_log,
                           y = hbs_olr_o)) +
  geom_point(aes(color = Country.Code),
             alpha = 0.6,
             size = 2) +
  geom_smooth(method = "lm",
              se = FALSE,
              color = "#C00000",
              linewidth = 0.8) +
  scale_color_viridis_d(option = "C") +
  labs(
    x = "Malaria incidence (log)",
    y = "HbS × OLR (shift–share IV)",
    color = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "bottom") + guides(col = guide_legend(nrow = 3))
fig8_corr_iv


# ggsave(
#   filename = "correlation_malaria_shift_share_by_country.png", 
#   plot = fig8_corr_iv,            
#   path = "Figure",
#   width = 10,                     
#   height = 6,                     
#   dpi = 600                       
# )


fig10_tfp_index <- ggplot(data,
                          aes(x = year)) +
  geom_line(aes(y = rtfpna,
                linetype = "Feenstra et al. (2015)",
                color    = "Feenstra et al. (2015)"),color = "#002395",
            linewidth = 0.7) +
  geom_line(aes(y = tfp_index,
                linetype = "Estimated",
                color    = "Estimated"),color = "#C00000",
            linewidth = 0.7) +
  facet_wrap(~ Country.Code) +
  labs(
    y = "TFP index",
    x = "Year",
    color = NULL,
    linetype = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "bottom",
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )
fig10_tfp_index


# ggsave(
#   filename = "TFP_index.png",  
#   plot = fig10_tfp_index,      
#   path = "Figure",
#   width = 10,                  
#   height = 6,                  
#   dpi = 600                    
# )

##Fig 11: tfps
tfps <- data.frame(
  factors = c("Malaria", "Years of Schooling", "Life Expectancy",
              "Undernourishment", "Lack Access Water",
              "Air Pollution", "HIV/AIDS"),
  values = c(-0.116, 0.087, 0.114, -0.140, -0.039, -0.073, -0.146),
  se = c(0.0416, 0.020, 0.036, 0.051, 0.012, 0.025, 0.051)
)

tfps <- tfps %>%
  mutate(
    lower_ci = values - 1.96 * se,
    upper_ci = values + 1.96 * se,
    is_malaria = factors == "Malaria"
  )

fig11_tfps <- ggplot(tfps,
                     aes(y = factors,
                         x = values)) +
  geom_vline(xintercept = 0,
             color = "black",
             linewidth = 0.4) +
  geom_errorbar(aes(xmin = lower_ci,
                    xmax = upper_ci,
                    color = is_malaria),   # 🔹 map color here
                width = 0.2,
                linewidth = 0.6) +
  geom_point(aes(color = is_malaria),
             size = 3) +
  scale_color_manual(values = c("TRUE" = "#C00000",
                                "FALSE" = "#002395"),
                     guide = "none") +
  coord_flip() +
  labs(
    x = "% change in TFP",
    y = "Factors"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 10)
  )

fig11_tfps

# ggsave(
#   filename = "tfps.png",  
#   plot = fig11_tfps,      
#   path = "Figure",
#   width = 10,                        
#   height = 6,                       
#   dpi = 600                        
# )

cat("\nScripts completed successfully.\n")








