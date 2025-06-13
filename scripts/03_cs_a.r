# Analysis a) to see whether an increased predicted age deviation (PAD) at baseline can differentiate
# between the clinical groups: Cognitive Normal (CN), Mild Cognitive Impaired (MCI), and Alzheimer's Disease (AD).
# This script runs an independent sample-t-test (or ANCOVA with covariates) to test for significant differences
# in PAD between the groups (CN v.s. MCI, CN v.s. AD, MCI v.s. AD)

# Load required libraries
library(readxl)
library(writexl)
library(here)
library(stats)
library(dplyr)
library(effsize)
library(data.table)


#### Functions and Definitions ####

#' Calculate Cohen's d effect size
#'
#' @param data A dataframe containing the data
#' @param group1 The name of the first group
#' @param group2 The name of the second group
#' @return A dataframe with effect size and confidence interval


calculate_cohens_d <- function(data, group1, group2) {
    data_group1 <- data %>%
        dplyr::filter(.data$diagnosis == group1) %>%
        dplyr::pull(.data$pad_model)

    data_group2 <- data %>%
        dplyr::filter(.data$diagnosis == group2) %>%
        dplyr::pull(.data$pad_model)

    result <- effsize::cohen.d(data_group1, data_group2)

    # eff-size diff between vs group 0 and vs group 1
    # regress out the effect of the co-vars, compar the residuals between vs 0 and vs 1
    # mtcars$res <- lm(data = mtcars, hp ~ qsec + drat)$residuals
    # effsize::cohen.d(mtcars$res[mtcars$vs == 1], mtcars$res[mtcars$vs == 0])

    return(
        data.frame(
            effect_size = result$estimate,
            lower_ci = result$conf.int[1],
            upper_ci = result$conf.int[2]
        )
    )
}

#' Run statistical test excluding a specific diagnosis
#'
#' @param data A dataframe containing the data
#' @param excluded_diagnosis The diagnosis to exclude from the analysis
#' @return A vector of test results for the diagnosis coefficient
run_statistical_test <- function(data, excluded_diagnosis, formula = "pad_model ~ diagnosis + aes + gender + chron_age") {
    filtered_data <- data %>%
        dplyr::filter(.data$diagnosis != excluded_diagnosis)

    formula <- as.formula(formula)
    model <- lm(formula, data = filtered_data)
    model_summary <- summary(model)

    print(coef(model_summary))

    test_result <- model_summary$coefficients[2, ]

    return(test_result)
}


#### Script ####

# Load excel file
file_path <- here::here("data", "cs_baseline.xlsx")
data <- readxl::read_excel(file_path)
print("Number of subjects:")
print(dim(data))

# Convert categorical variables to factors
baseline_data <- data %>%
    dplyr::mutate(
        subject_id = factor(subject_id),
        gender = factor(gender),
        diagnosis = factor(diagnosis, levels = c("CN", "MCI", "AD"))
    )

# Select columns
dvs <- c(
    "pad_brainageR",
    "pad_DeepBrainNet",
    "pad_brainage",
    "pad_enigma",
    "pad_pyment",
    "pad_mccqrnn",
    "gm_icv"
)
idx_dvs <- which(colnames(baseline_data) %in% dvs)

# Initialize lists to store summary and results
results_list <- list()
plot_list <- list()

for (i in seq_along(dvs)) {
    message(sprintf("Running model %s", dvs[i]))

    # get index of current model
    idx <- idx_dvs[i]

    # rename the column of this index so that we can run the same lines of code for all models
    baseline_data_model <- baseline_data %>%
        dplyr::select(pad_model = all_of(idx), diagnosis, chron_age, aes, gender)

    # Compute effect sizes between groups
    effect_size_cn_mci <- calculate_cohens_d(baseline_data_model, "CN", "MCI")
    effect_size_cn_ad <- calculate_cohens_d(baseline_data_model, "CN", "AD")
    effect_size_mci_ad <- calculate_cohens_d(baseline_data_model, "MCI", "AD")

    # Combine into dataframe
    print("Effect sizes:")
    effect_sizes <- dplyr::bind_rows(
        effect_size_cn_mci %>% dplyr::mutate(comparison = "CN_MCI", model = dvs[i]),
        effect_size_cn_ad %>% dplyr::mutate(comparison = "CN_AD", model = dvs[i]),
        effect_size_mci_ad %>% dplyr::mutate(comparison = "MCI_AD", model = dvs[i])
    )
    print(effect_sizes)

    # run tests for sanity check without covariates:coefficients
    print("Tests without covariates:")
    test_cn_ad_no_cov <- run_statistical_test(baseline_data_model, "AD", formula = "pad_model ~ diagnosis")
    test_cn_mci_no_cov <- run_statistical_test(baseline_data_model, "MCI", formula = "pad_model ~ diagnosis")
    test_mci_ad_no_cov <- run_statistical_test(baseline_data_model, "CN", formula = "pad_model ~ diagnosis")

    tests_no_cov <- dplyr::bind_rows(test_cn_ad_no_cov, test_cn_mci_no_cov, test_mci_ad_no_cov)

    # Run tests for significant differences in means between groups
    print("Tests with covariates:")
    test_cn_mci <- run_statistical_test(baseline_data_model, "AD")
    test_cn_ad <- run_statistical_test(baseline_data_model, "MCI")
    test_mci_ad <- run_statistical_test(baseline_data_model, "CN")

    tests <- dplyr::bind_rows(test_cn_mci, test_cn_ad, test_mci_ad)

    results <- dplyr::bind_cols(effect_sizes, tests) %>%
        dplyr::select(
            model,
            comparison,
            effect_size,
            lower_ci,
            upper_ci,
            Estimate,
            `Std. Error`,
            `t value`,
            `Pr(>|t|)`
        )

    results_list[[i]] <- results
}

# Combine stats
total_result <- dplyr::bind_rows(results_list)
total_result_rounded <- total_result %>% dplyr::mutate(across(where(is.double), \(x) round(x, digits = 2)))
total_result_rounded["Pr(>|t|)"] <- total_result["Pr(>|t|)"]

# Store results
store_results_path <- here::here("results", "cs_a_results.xlsx")
writexl::write_xlsx(total_result_rounded, path = store_results_path)
