### Summary stats for the cross sectional analysis
# Load required libraries
library(readxl)
library(writexl)
library(here)
library(dplyr)
library(tidyr)

# Load excel file
file_path <- here::here("data", "megamastersheet.xlsx")
data <- readxl::read_excel(file_path)
print(dim(data))
# [1] 10802    36

# only keep the columns indicated in column keep
data <- data[data$keep == 1, ]
print(dim(data))
# [1] 9768   36

# Convert categorical variables to factors
data <- data %>%
    dplyr::mutate(
        subject_id = factor(subject_id),
        mri_field = factor(mri_field),
        session_id = factor(session_id),
        gender = factor(gender),
        diagnosis = factor(diagnosis, levels = c("CN", "MCI", "AD"))
    )

#### Get baseline data
baseline_data <- data %>%
    dplyr::filter(session_id == "M000")
print(dim(baseline_data))
# [1] 2330   36

# Subset baseline data
baseline_data_subset <- baseline_data %>%
    dplyr::select(
        subject_id,
        starts_with("pad_"),
        chron_age,
        education_level,
        diagnosis,
        gender,
        aes,
        gm,
        icv,
        gm_icv,
        ADNI_MEM
    )

# get number of nans
nans <- baseline_data_subset %>%
    dplyr::summarise_all(~ sum(is.na(.)))
print(nans)

# Remove rows with missing values
baseline_data_subset <- na.omit(baseline_data_subset)
print(dim(baseline_data_subset))
# [1] 2330   15

# store data
store_data_path <- here::here("data", "cs_baseline.xlsx")
writexl::write_xlsx(baseline_data_subset, path = store_data_path)

### get means and confidence intervals per group for each pad_* by fitting a lm without covariates
# Subset the dataset to keep only columns starting with "pad_" + gm_icv + diagnosis
df_subset <- baseline_data_subset %>% dplyr::select(starts_with("pad_"), gm_icv, diagnosis)

## compute the mean and SD for each diagnostic group
summary_stats <- df_subset %>%
  group_by(diagnosis) %>%
  summarise(across(where(is.numeric), 
                   list(mean = ~ round(mean(.), 2), sd = ~ round(sd(.), 2)), 
                   .names = "{.col}_{.fn}"))

# Write a function to get the coefficients and confidence intervals
get_lm_stats <- function(x) {
  model <- lm(x ~ 1)
  ci <- confint(model)
  stats <- c(mean = coef(model), lower_ci = ci[1], upper_ci = ci[2])
  names(stats) <- c("mean", "lower_ci", "upper_ci")
  return(stats)
}

summary_stats_lm <- df_subset %>%
  group_by(diagnosis) %>%
  summarise(across(where(is.numeric), 
                   list(mean = ~ get_lm_stats(.x)["mean"], lower_ci = ~ get_lm_stats(.x)["lower_ci"], upper_ci = ~ get_lm_stats(.x)["upper_ci"]), 
                   .names = "{.col}_{.fn}"))

summary_stats_lm <- summary_stats_lm %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))

# compute the std for each variable as well
sds <- df_subset %>%
  group_by(diagnosis) %>%
  summarise(across(where(is.numeric), 
                   list(sd = ~ round(sd(.x), 2)), 
                   .names = "{.col}_{.fn}"))

# Combine variance with the summary statistics
summary_stats_lm_sd <- summary_stats_lm %>%
  left_join(sds, by = "diagnosis")

# Compute the variance for brainageR
brainageR_sd <- sd(df_subset$pad_brainageR[df_subset$diagnosis == "CN"], na.rm = TRUE)
brainageR_mean <- mean(df_subset$pad_brainageR[df_subset$diagnosis == "CN"], na.rm = TRUE)
print(paste("SD of brainageR:", round(brainageR_sd, 2)))
print(paste("Mean of brainageR:", round(brainageR_mean, 2)))

# View results
print(t(summary_stats_lm_sd))

# Store results
store_summary_path <- here::here("results", "cs_summary.xlsx")
writexl::write_xlsx(summary_stats_lm_sd, path = store_summary_path)
