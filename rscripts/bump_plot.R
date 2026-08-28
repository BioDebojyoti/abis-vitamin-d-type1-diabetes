rm(list = ls())

library(dplyr)
library(tidyr)
library(utils)
library(ggplot2)

library(survival)
library(survminer)

data_dir = "/Users/debda22/OneDrive - Linköpings universitet/T1D_vitamin_D_results/cox_analysis/concentrations/"

df_cox_full <- read.table(paste0(data_dir,"age20/long_format_cox_age_full.csv"), header = TRUE) %>%
  dplyr::mutate(scenario= 4)
df_cox_age_5 <- read.table(paste0(data_dir,"age5/long_format_cox_age_5.csv"), header = TRUE) %>%
  dplyr::mutate(scenario= 1)
df_cox_age_10 <- read.table(paste0(data_dir,"age10/long_format_cox_age_10.csv"), header = TRUE) %>%
  dplyr::mutate(scenario= 2)
df_cox_age_15 <- read.table(paste0(data_dir,"age15/long_format_cox_age_15.csv"), header = TRUE) %>%
  dplyr::mutate(scenario= 3)

# Reshape data to long format
df_long <- rbind(
  df_cox_full,
  df_cox_age_5,
  df_cox_age_10,
  df_cox_age_15
) %>%
  dplyr::mutate(
    ABISnr = factor(ABISnr),
    event = factor(event),
    origin=0
  )

# Order scenarios
df_long$scenario <- factor(
  df_long$scenario,
  levels = c(1,2,3,4),
  labels = c("~5 age", "~10 age", "~15 age", "~full follow-up")
  )

censoring_design <- df_long %>%
  ggplot(aes(y = ABISnr)) +
  # all participants baseline in light blue
  geom_segment(aes(x = origin, xend = exit_age, yend = ABISnr),
               color = "skyblue1", alpha = 0.2, linewidth = 0.1) +

  # overlay red for events
  geom_segment(data = df_long %>% filter(event == 1),
               aes(x = origin, xend = exit_age, y = ABISnr, yend = ABISnr),
               color = "red", alpha = 0.6, linewidth = 0.3) +

  # start point
  geom_point(aes(x = origin), color = "gray70", size = 0.3) +

  # end point
  geom_point(aes(x = exit_age, color = factor(event)), size = 1.2) +

  facet_wrap(~ scenario,
             ncol = 1,
             scales = "free_y",
             labeller = as_labeller(
               c(
                 "~full follow-up" = "Full follow-up; events = 83; censored = 296",
                 "~15 age"   = "Censored at 15 yrs; events = 60; censored = 319",
                 "~10 age"   = "Censored at 10 yrs; events = 35; censored = 344",
                 "~5 age"    = "Censored at  5 yrs; events = 12; censored = 367"
               )
             )) +
  scale_color_manual(values = c("0" = "gray40", "1" = "red"),
                     labels = c("Censored", "Event"),
                     name = "Status") +
  labs(x = "Follow-up time (years)",
       y = "Participant (ABISnr)"#,
       # title = "Follow-up intervals under different censoring strategies"
       ) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank())


ggplot2::ggsave(
  "/Users/debda22/OneDrive - Linköpings universitet/T1D_vitamin_D_results/cox_analysis/concentrations/censoring_design.pdf",
  censoring_design,
  height = 11.69,
  width = 8.27,
  units = "in",
  dpi = 300,
  font = "Times"
)



model1 <- coxph(
  Surv(exit_age, event) ~ conc_25OHD,
  data = df_cox_age_5
  )

model2 <- coxph(
  Surv(exit_age, event) ~ conc_25OHD,
  data = df_cox_age_10
  )

model3 <- coxph(
  Surv(exit_age, event) ~ conc_25OHD,
  data = df_cox_age_15
)

model4 <- coxph(
  Surv(exit_age, event) ~ conc_25OHD,
  data = df_cox_full
)

med_val <- sapply(c("~5 age", "~10 age", "~15 age", "~full follow-up"), function(x) {
  median((df_long %>% dplyr::filter(scenario == x & event == 0))$conc_25OHD, na.rm = TRUE)
})

survfits <- list(
  "~5 age"  = survfit(model1, newdata = data.frame(conc_25OHD = med_val[[1]])),
  "~10 age" = survfit(model2, newdata = data.frame(conc_25OHD = med_val[[2]])),
  "~15 age" = survfit(model3, newdata = data.frame(conc_25OHD = med_val[[3]])),
  "~full follow-up" = survfit(model4, newdata = data.frame(conc_25OHD = med_val[[4]]))
)

combined_plot <- ggsurvplot_combine(
  fit = survfits,
  legend.title = "Scenario",
  legend.labs = names(survfits),
  palette = c("#1b9e77", "#d95f02", "#7570b3", "#e7298a"),
  risk.table = FALSE
)
