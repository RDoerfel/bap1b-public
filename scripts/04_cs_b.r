# Analysis b) an increased PAD is associated with decreased cognition at baseline
# Load required libraries
library(readxl)
library(writexl)
library(here)
library(stats)
library(dplyr)
library(lgr)
library(data.table)
library(ppcor)
library(LMMstar)

# Function to create structured log entry with multi-line content
log_structured_multiline <- function(logger, content, type = "model_summary") {
    output <- capture.output(content)
    
    # Create a properly formatted multi-line string
    formatted_content <- paste(output, collapse = "\n")
  
    # Log with the content as a separate field
    logger$info(
        sprintf("=== %s ===", type),
        content = sprintf("\n%s\n", formatted_content)
    )
}

# Initialize logger if not already initialized
if (!exists("logger")) {
    logger <- lgr::get_logger("logger")
    if (length(logger$appenders) == 0) {
        logger$set_appenders(AppenderConsole$new())
    }
}

#### Functions ####
run_statistical_test <- function(data, included_diagnosis, formula = "ADNI_MEM ~ pad_model + aes + gender + chron_age + education_level") {
    logger$info(glue::glue("Running model for {formula}"))

    model <- lm(formula, data = data)

    return(model)
}

compute_partial_corellation <- function(test_result) {
    # Compute the partial correlation coefficient
    # eq. 3 from Lipsitz, T. et al. (2001)

    # z_k = beta_k / se_k
    z_k <- test_result['Estimate'] / test_result['Std. Error']
    df <- test_result['df']
    p_k <- z_k / sqrt(df) / sqrt(1 + z_k^2 / df)

    # it should get the name pc_k
    names(p_k) <- "pc_k"

    return(p_k)
}


compute_pck <- function(data, covariates) {
    data_df <- as.data.frame(data)
    # Convert gender to binary
    data_df$gender <- ifelse(data_df$gender == "Male", 1, 0)
    pcor_result <- ppcor::pcor.test(x = data_df$pad_model,
                                    y = data_df$ADNI_MEM,
                                    z = data_df[, covariates],
                                    method = "pearson")
    results <- data.frame(
        pcc = pcor_result$estimate,
        pcc_p.value = pcor_result$p.value
    )
    return(results)
} 


format_model_summary <- function(model) {
    model_summary <- summary(model)
    test_result <- model_summary$coefficients[2, ]

    # add model$df to the test_result 
    test_result <- c(test_result, df = as.integer(model$df.residual))
    names(test_result) <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)", "df")

    return(test_result)
}

format_par_cor_results <- function(parcor) {
    # Extract the relevant information from the partial correlation result
    test_result <- data.frame(
        parcor = parcor$estimate,
        parcor_se = parcor$se,
        parcor_df = parcor$df,
        parcor_lower = parcor$lower,
        parcor_upper = parcor$upper,
        parcor_p.value = parcor$p.value
    )

    return(test_result)
}


# Function to perform ANOVA and extract F and p-values
perform_anova <- function(model1, model2, logger) {
    anova_result <- anova(model1, model2)
    log_structured_multiline(logger, anova_result, type = "anova_result")

    f_value <- anova_result$`F`[2]
    p_value <- anova_result$`Pr(>F)`[2]
    
    # Log the ANOVA result
    logger$info(glue::glue("ANOVA Result: F value = {f_value}, Pr(>F) = {p_value}"))

    return(list(f_value = f_value, p_value = p_value))
}

#### Scripts ####
file_path <- here::here("data", "cs_baseline.xlsx")
data <- readxl::read_excel(file_path)
logger$info(glue::glue("Number of subjects: {dim(data)}"))

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
results_list <- list()

total_preds_cn <- data.table()
total_preds_mci <- data.table()
total_preds_ad <- data.table()

for (i in seq_along(dvs)) {
    logger$info(glue::glue("#### Running model {dvs[i]} #### \n"))

    # get index of current model
    idx <- idx_dvs[i]

    # rename the column of this index so that we can run the same lines of code for all models
    baseline_data_model <- baseline_data %>%
        dplyr::select(pad_model = all_of(idx), diagnosis, chron_age, aes, gender, ADNI_MEM, education_level)


    # Cognitive normal
    logger$info("=== Run for CN === \n")

    data_cn <- baseline_data_model %>% dplyr::filter(baseline_data_model$diagnosis == "CN")
    test_cn <- run_statistical_test(data_cn, "CN", "ADNI_MEM ~ pad_model + aes + gender + chron_age")
    test_cn_edu <- run_statistical_test(data_cn, "CN", "ADNI_MEM ~ pad_model + aes + gender + chron_age + education_level")

    test_results_cn <- dplyr::bind_rows(
        format_model_summary(test_cn),
        format_model_summary(test_cn_edu)
    )

    # compute partial correlation
    test_results_cn <- test_results_cn %>%
        dplyr::mutate(compute_partial_corellation(test_results_cn))

    # use the ppcor package to compute the partial correlation coefficient
    pcc_cn <- compute_pck(data_cn, c("aes", "gender", "chron_age" ))
    pcc_cn_edu <- compute_pck(data_cn, c("aes", "gender", "chron_age", "education_level"))

    pcc_cn <- dplyr::bind_rows(pcc_cn, pcc_cn_edu)

    test_results_cn <- dplyr::bind_cols(test_results_cn, pcc_cn)

    # Use LMMStar with CIs
    par_cor_cn <- partialCor(c(pad_model, ADNI_MEM) ~ aes + gender + chron_age, data = data_cn)
    par_cor_cn_edu <- partialCor(c(pad_model, ADNI_MEM) ~ aes + gender + chron_age + education_level, data = data_cn)

    # bind rows
    par_cor_cn <- dplyr::bind_rows(
        format_par_cor_results(par_cor_cn),
        format_par_cor_results(par_cor_cn_edu)
    )

    # add the partial correlation to the test results
    test_results_cn <- dplyr::bind_cols(test_results_cn, par_cor_cn)  

    # runn anova
    anova_cn <- perform_anova(test_cn, test_cn_edu, logger)

    # add to test results 
    test_results_cn <- test_results_cn %>% dplyr::mutate(
        "F value" = c("-", anova_cn$f_value),
        "Pr(>|F|)" = c("-", anova_cn$p_value)
    )

    log_structured_multiline(logger, test_results_cn, type = "test_result")

    # get the ci for the prediction using inbuild functions from R. fix the covariates to their means, iterate over min/max of pad_model

    data_test_cn <- data.frame(
        pad_model = seq(min(data_cn$pad_model), max(data_cn$pad_model), length.out = 100),
        aes = mean(data_cn$aes),
        gender = "Male",
        chron_age = mean(data_cn$chron_age)
    )

    pred_ci_cn <- predict(test_cn, newdata = data_test_cn, interval = "confidence", level = 0.95)
    colnames(pred_ci_cn) <- paste0(dvs[i], "_", c("fit", "lwr", "upr"))
    total_preds_cn <- cbind(total_preds_cn, pred_ci_cn)

    # Mild cognitive impairment
    logger$info("=== Run for MCI === \n")
    data_mci <- baseline_data_model %>% dplyr::filter(baseline_data_model$diagnosis == "MCI")
    test_mci <- run_statistical_test(data_mci, "MCI", "ADNI_MEM ~ pad_model + aes + gender + chron_age")
    test_mci_edu <- run_statistical_test(data_mci, "MCI", "ADNI_MEM ~ pad_model + aes + gender + chron_age + education_level")

    test_results_mci <- dplyr::bind_rows(
        format_model_summary(test_mci),
        format_model_summary(test_mci_edu)
    )

    # compute partial correlation
    test_results_mci <- test_results_mci %>%
        dplyr::mutate(compute_partial_corellation(test_results_mci))

    # use the ppcor package to compute the partial correlation coefficient
    pcc_mci <- compute_pck(data_mci, c("aes", "gender", "chron_age" ))
    pcc_mci_edu <- compute_pck(data_mci, c("aes", "gender", "chron_age", "education_level"))

    pcc_mci <- dplyr::bind_rows(pcc_mci, pcc_mci_edu)

    test_results_mci <- dplyr::bind_cols(test_results_mci, pcc_mci)

    # Use LMMStar with CIs
    par_cor_mci <- partialCor(c(pad_model, ADNI_MEM) ~ aes + gender + chron_age, data = data_mci)
    par_cor_mci_edu <- partialCor(c(pad_model, ADNI_MEM) ~ aes + gender + chron_age + education_level, data = data_mci)

    # bind rows
    par_cor_mci <- dplyr::bind_rows(
        format_par_cor_results(par_cor_mci),
        format_par_cor_results(par_cor_mci_edu)
    )

    # add the partial correlation to the test results
    test_results_mci <- dplyr::bind_cols(test_results_mci, par_cor_mci)
    
    # runn anova
    anova_mci <- perform_anova(test_mci, test_mci_edu, logger)

    # add to test results
    test_results_mci <- test_results_mci %>% dplyr::mutate(
        "F value" = c("-", anova_mci$f_value),
        "Pr(>|F|)" = c("-", anova_mci$p_value)
    )

    log_structured_multiline(logger, test_results_mci, type = "test_result")

    # get the ci for the prediction using inbuild functions from R. fix the covariates to their means, iterate over min/max of pad_model

    data_test_mci <- data.frame(
        pad_model = seq(min(data_mci$pad_model), max(data_mci$pad_model), length.out = 100),
        aes = mean(data_mci$aes),
        gender = "Male",
        chron_age = mean(data_mci$chron_age)
    )

    pred_ci_mci <- predict(test_mci, newdata = data_test_mci, interval = "confidence", level = 0.95)
    colnames(pred_ci_mci) <- paste0(dvs[i], "_", c("fit", "lwr", "upr"))
    total_preds_mci <- cbind(total_preds_mci, pred_ci_mci)

    # Alzheimer's disease
    logger$info("=== Run for AD === \n")
    data_ad <- baseline_data_model %>% dplyr::filter(baseline_data_model$diagnosis == "AD")
    test_ad <- run_statistical_test(data_ad, "AD", "ADNI_MEM ~ pad_model + aes + gender + chron_age ")
    test_ad_edu <- run_statistical_test(data_ad, "AD", "ADNI_MEM ~ pad_model + aes + gender + chron_age + education_level")

    test_results_ad <- dplyr::bind_rows(
        format_model_summary(test_ad),
        format_model_summary(test_ad_edu)
    )

    # compute partial correlation
    test_results_ad <- test_results_ad %>%
        dplyr::mutate(compute_partial_corellation(test_results_ad))

    # use the ppcor package to compute the partial correlation coefficient
    pcc_ad <- compute_pck(data_ad, c("aes", "gender", "chron_age" ))
    pcc_ad_edu <- compute_pck(data_ad, c("aes", "gender", "chron_age", "education_level"))
    
    pcc_ad <- dplyr::bind_rows(pcc_ad, pcc_ad_edu)
    test_results_ad <- dplyr::bind_cols(test_results_ad, pcc_ad)

    # Use LMMStar with CIs
    par_cor_ad <- partialCor(c(pad_model, ADNI_MEM) ~ aes + gender + chron_age, data = data_ad)
    par_cor_ad_edu <- partialCor(c(pad_model, ADNI_MEM) ~ aes + gender + chron_age + education_level, data = data_ad)

    # bind rows
    par_cor_ad <- dplyr::bind_rows(
        format_par_cor_results(par_cor_ad),
        format_par_cor_results(par_cor_ad_edu)
    )

    # add the partial correlation to the test results
    test_results_ad <- dplyr::bind_cols(test_results_ad, par_cor_ad)

    # runn anova
    anova_ad <- perform_anova(test_ad, test_ad_edu, logger)

    # add to test results
    test_results_ad <- test_results_ad %>% dplyr::mutate(
        "F value" = c("-", anova_ad$f_value),
        "Pr(>|F|)" = c("-", anova_ad$p_value)
    )

    log_structured_multiline(logger, test_results_ad, type = "test_result")

    data_test_ad <- data.frame(
        pad_model = seq(min(data_ad$pad_model), max(data_ad$pad_model), length.out = 100),
        aes = mean(data_ad$aes),
        gender = "Male",
        chron_age = mean(data_ad$chron_age)
    )

    pred_ci_ad <- predict(test_ad, newdata = data_test_ad, interval = "confidence", level = 0.95)
    colnames(pred_ci_ad) <- paste0(dvs[i], "_", c("fit", "lwr", "upr"))
    total_preds_ad <- cbind(total_preds_ad, pred_ci_ad)

    # append results with diagnosis and model
    test_results_cn <- test_results_cn %>% dplyr::mutate(pad = dvs[i], group = "CN", type = c("model", "model+edu"))
    test_results_mci <- test_results_mci %>% dplyr::mutate(pad = dvs[i], group = "MCI", type = c("model", "model+edu"))
    test_results_ad <- test_results_ad %>% dplyr::mutate(pad = dvs[i], group = "AD", type = c("model", "model+edu"))

    # combine rows
    test_results_list <- list(test_results_cn, test_results_mci, test_results_ad)

    results_list[[i]] <- test_results_list
}

# combine results
results_cs_b <- dplyr::bind_rows(results_list)

# reorder columns
results_cs_b <- results_cs_b %>% dplyr::select(pad, type, group, everything())

# combine predictions pred_ci_cn, pred_ci_mci, pred_ci_ad. Add a column diagnosis to each one first
total_preds_cn[, diagnosis := "CN"]
total_preds_mci[, diagnosis := "MCI"]
total_preds_ad[, diagnosis := "AD"]

total_preds <- rbind(total_preds_cn, total_preds_mci, total_preds_ad)


# Store results
store_summary_path <- here::here("results", "cs_b_results.xlsx")
writexl::write_xlsx(results_cs_b, path = store_summary_path)
logger$info(glue::glue("Results stored at {store_summary_path}"))

# store predictions
store_preds_path <- here::here("results", "cs_b_predictions.xlsx")
writexl::write_xlsx(total_preds, path = store_preds_path)
logger$info(glue::glue("Predictions stored at {store_preds_path}"))
