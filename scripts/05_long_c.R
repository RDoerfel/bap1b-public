# Analysis c) to investigate whether an increased PAD at baseline is associated with future conversion from CN/MCI to MCI/AD.
# Load required libraries
library(mvtnorm)
library(data.table)
library(LMMstar)
library(mets)
library(riskRegression)
library(dplyr)
library(lava)

### Functions
calc_var_lp <- function(vec_values, vcov) {
    # compute the variance of the linear predictor
    # lp = b0 + b1*x1 + b2*x2 + ... + bn*x3 = t(b)x
    # var(t(b)x) = t(x) var(b) x, with x being a column vector
    return(t(vec_values) %*% vcov %*% vec_values)
}

pred_lp <- function(vec_coeffs, vec_values) {
    # compute linear combination of values and coefficients to get the linear predictor
    # t(b)X = sum(b_i * X_i)
    return(t(vec_coeffs) %*% vec_values)
}

calc_ci_95 <- function(prediction_lp, var_lp) {
    # compute the confidence interval for the prediction
    lower <- plogis(prediction_lp - 1.96 * sqrt(var_lp))
    upper <- plogis(prediction_lp + 1.96 * sqrt(var_lp))
    return(c(lower, upper))
}

predict_and_ci95 <- function(vec_coeffs, vcov, vec_values) {
    var_lp <- calc_var_lp(vec_values, vcov)
    pred_lp <- pred_lp(vec_coeffs, vec_values)
    ci <- calc_ci_95(pred_lp, var_lp)
    prob <- plogis(pred_lp)

    # return a vector with the prediction and the confidence interval and names
    result <- c(prob, ci)
    names(result) <- c("prob", "ci_lower", "ci_upper")

    return(result)
}

### Script

# Load excel file
file_path <- here::here("data", "megamastersheet.xlsx")
data <- readxl::read_excel(file_path)
print(dim(data))
# [1] 10802    36

# only take data that is indicated by column keep
data_keep <- data[data$keep == 1, ]
print(dim(data_keep))
# [1] 9768   36

# Convert categorical variables to factors
data_keep <- data_keep %>%
    dplyr::mutate(
        subject_id = factor(subject_id),
        session_id = factor(session_id),
        mri_field = factor(mri_field),
        session_id = factor(session_id),
        gender = factor(gender),
        diagnosis = factor(diagnosis, levels = c("CN", "MCI", "AD"))
    )

# make it  a table
data_keep <- data.table(data_keep)

# select columns
data_subset <- data_keep %>% dplyr::select(
    subject_id,
    session_id,
    gender,
    diagnosis,
    aes,
    gm_icv,
    starts_with("pad_"),
    chron_age,
    time_from_baseline
)

# only until 4 years
data_subset <- data_subset %>% dplyr::filter(time_from_baseline <= 4)
print("Number of scans until 4 years:")
print(dim(data_subset))
# [1] 8067   15

# count nas
na_counts <- data_subset %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), ~ sum(is.na(.))))

# Print the result
print("Number of missing values:")
print(na_counts)

## Only look at CN that potentially transition later on
data_subset <- data_subset %>%
    dplyr::arrange(subject_id, time_from_baseline) %>%
    dplyr::group_by(subject_id) %>%
    dplyr::filter(dplyr::first(diagnosis) == "CN") %>%
    dplyr::ungroup()

print("Number of scans after filtering for CN at baseline:")
print(dim(data_subset))
# [1] 2637   15

# Get baseline values
data_subset_base <- data_subset %>%
    dplyr::group_by(subject_id) %>%
    dplyr::arrange(subject_id, session_id) %>%
    dplyr::mutate(
        pad_brainageR_base = dplyr::first(pad_brainageR),
        pad_DeepBrainNet_base = dplyr::first(pad_DeepBrainNet),
        pad_brainage_base = dplyr::first(pad_brainage),
        pad_enigma_base = dplyr::first(pad_enigma),
        pad_pyment_base = dplyr::first(pad_pyment),
        pad_mccqrnn_base = dplyr::first(pad_mccqrnn),
        gm_icv_base = dplyr::first(gm_icv),
        aes_base = dplyr::first(aes),
        age_base = dplyr::first(chron_age),
        diagnosis_base = dplyr::first(diagnosis),
        diagnosis_change = diagnosis_base != diagnosis,
        change_occurred = sum(diagnosis_change) > 0,
        last_false = cumsum(diagnosis_change == FALSE)
    ) %>%
    dplyr::ungroup()

# subset to only take last available entry
data_subset_base_survive <- data_subset_base %>%
    dplyr::group_by(subject_id) %>%
    dplyr::filter(
        (any(change_occurred) & cumsum(diagnosis_change) == 1) | (!any(change_occurred) & row_number() == n())
    )

print("Number of scans entering the analysis")
print(dim(data_subset_base_survive))
# [1] 875  29

# compute "prevelance" and number of subjects with change
print("Prevelance of change from CN to MCI/AD")
print(mean(data_subset_base_survive$diagnosis_change))
# [1] 0.07428571

print("Number of subjects with change from CN to MCI/AD")
print(sum(data_subset_base_survive$diagnosis_change))
# [1] 65

# if it is the last entry and change did not occur yet, set time to shortly after visit
data_subset_base_survive_timing <- data_subset_base_survive %>%
    dplyr::mutate(survive_timing = dplyr::case_when(
        session_id == "M048" & change_occurred == FALSE ~ time_from_baseline + 1e-5,
        TRUE ~ time_from_baseline
    ))

# if a subject drops out at the first visit, set time to shortly after visit
data_subset_base_survive_timing <- data_subset_base_survive_timing %>%
    dplyr::mutate(survive_timing = dplyr::case_when(
        survive_timing == 0 ~ 1e-5,
        TRUE ~ survive_timing
    ))

# convert to table
data_subset_base_survive_timing <- data.table(data_subset_base_survive_timing)

# save to excel
store_path <- here::here("data", "data_long_c.xlsx")
writexl::write_xlsx(data_subset_base_survive_timing, path = store_path)

# Select columns
dvs <- c(
    "pad_brainageR_base",
    "pad_DeepBrainNet_base",
    "pad_brainage_base",
    "pad_enigma_base",
    "pad_pyment_base",
    "pad_mccqrnn_base",
    "gm_icv_base"
)

idx_dvs <- which(colnames(data_subset_base_survive_timing) %in% dvs)

# add columns model, coef, intercept, p_value coef, p_value intercept
total_result <- data.table(
    model = character(),
    odds = numeric(),
    odds_lower = numeric(),
    odds_upper = numeric(),
    odds_p_value = numeric(),
    prob_pad_minus_sd = numeric(),
    ci_lower_pad_minus_sd = numeric(),
    ci_upper_pad_minus_sd = numeric(),
    prob_pad_mu = numeric(),
    ci_lower_pad_mu = numeric(),
    ci_upper_pad_mu = numeric(),
    prob_pad_plus_sd = numeric(),
    ci_lower_pad_plus_sd = numeric(),
    ci_upper_pad_plus_sd = numeric(),
    auc = numeric(),
    auc_lower = numeric(),
    auc_upper = numeric(),
    brier = numeric(),
    brier_lower = numeric(),
    brier_upper = numeric()
)

total_preds <- data.table()

wglm_models <- list()

for (i in seq_along(dvs)) {
    message(sprintf("Running model %s", dvs[i]))

    # get index of current model
    idx <- idx_dvs[i]

    # rename the column of this index so that we can run the same lines of code for all models
    data_model <- data_subset_base_survive_timing %>%
        dplyr::select(subject_id,
            session_id,
            gender,
            age_base,
            aes_base,
            diagnosis_base,
            diagnosis_change,
            survive_timing,
            pad_base = all_of(idx),
        )

    # conver diagnosis change to 1 or 0
    data_model$diagnosis_change <- as.numeric(data_model$diagnosis_change)
    # center data
    data_model$pad_baseC <- data_model$pad_base - mean(data_model$pad_base)
    data_model$aes_baseC <- data_model$aes_base - mean(data_model$aes_base)
    data_model$age_baseC <- data_model$age_base - mean(data_model$age_base)

    # data_model$pad_baseC <- data_model$pad_base
    # data_model$aes_baseC <- data_model$aes_base
    # data_model$age_baseC <- data_model$age_base

    ### Analysis
    ## First, fake the complete case (no dropouts)
    print("Fake complete case analysis:")
    data_model_nc <- data_model[survive_timing >= 4.0]
    print("Number of scans:")
    print(dim(data_model_nc))
    # [1] 170  12

    print("Prevelance of change from CN to MCI/AD")
    print(mean(data_model_nc$diagnosis_change))
    # [1] 0.05294118

    # naive / glm model. should give the same results as wglm for complete case
    print("Naive model")
    nc_glm <- glm(diagnosis_change ~ pad_baseC + age_baseC + aes_baseC + gender,
        family = binomial(link = "logit"),
        data = data_model_nc
    )

    print(summary(nc_glm)$coef)

    print(plogis(cumsum(coef(nc_glm))))

    # wglm (without censoring, should have the same result)
    print("wglm")
    nc_wglm <- wglm(Surv(survive_timing, diagnosis_change) ~ pad_baseC + age_baseC + aes_baseC + gender,
        times = 4,
        data = data_model_nc,
        formula.censor = ~ pad_baseC + age_baseC + aes_baseC + gender,
        product.limit = TRUE,
        fitter = "coxph",
    )

    print(summary(nc_wglm)$coef)

    print(plogis(cumsum(coef(nc_wglm))))

    print(quantile(weights(nc_wglm)))

    print(model.tables(nc_wglm))

    hist(weights(nc_wglm))

    # AUC
    print("AUC")
    score_nc <- Score(
        list(glm = nc_glm, wglm = nc_wglm),
        formula = Surv(survive_timing, diagnosis_change) ~ pad_baseC + age_baseC + aes_baseC + gender,
        data = data_model_nc,
        metrics = c("AUC", "Brier"),
        plot = c("ROC"),
        case = 1,
        times = 4
    )

    # auc
    print(score_nc$AUC)
    print(score_nc$Brier)
    plotROC(score_nc)

    ## Prevelance of change from CN to MCI/AD
    print("Actual analysis accounting for censoring:")
    print(mean(data_model$diagnosis_change))
    # [1] 0.07411631

    ## Standard Analysis (naive)
    c_glm <- glm(diagnosis_change ~ pad_baseC + age_baseC + aes_baseC + gender,
        family = binomial(link = "logit"),
        data = data_model
    )

    print(summary(c_glm)$coef)
    print(plogis(cumsum(coef(c_glm))))
    print(confint(c_glm))

    ## censoring wglm
    c_wglm <- wglm(Surv(survive_timing, diagnosis_change) ~ pad_baseC + age_baseC + aes_baseC + gender, 
        times = 4,
        data = data_model,
        formula.censor = ~pad_baseC + age_baseC + aes_baseC + strat(gender),
        product.limit = TRUE,
        fitter = "coxph",
    )

    # add to list
    wglm_models[[dvs[i]]] <- c_wglm

    print(model.tables(c_wglm))
    print(quantile(weights(c_wglm)))
    print(confint(c_wglm))
    hist(weights(c_wglm))

    print(plogis(cumsum(coef(c_wglm))))

    # AUC
    print("AUC")
    score <- Score(
        list(glm = c_glm, wglm = c_wglm),
        times = 4,
        formula = Surv(survive_timing, diagnosis_change) ~ pad_baseC + age_baseC + aes_baseC + gender,
        data = data_model,
        metrics = c("AUC", "Brier"),
        plot = c("ROC")
    )

    print(score$AUC)
    print(score$Brier)
    plotROC(score)

    # get variance and covariance for PAD = mu_pad 
    c_wglm_vcov <- crossprod(lava::iid(c_wglm)[])
    vec_coeffs <- coef(c_wglm)
    gender <- "F"
    vec_values <- c(1, mean(data_model$pad_baseC), mean(data_model$age_baseC), mean(data_model$aes_baseC), gender == gender)

    # compute the prediction and confidence interval
    predict_and_ci_mu <- predict_and_ci95(vec_coeffs, c_wglm_vcov, vec_values)

    # compute prob of change with pad_baseC = + sd_pad
    # if dvs is gm_icv, we need to adjust the values
    vec_values[2] <- sd(data_model$pad_baseC)

    predict_and_ci_plus_sd <- predict_and_ci95(vec_coeffs, c_wglm_vcov, vec_values)

    # print the results
    print("Predictions and confidence intervals")
    print(predict_and_ci_plus_sd)

    # predict compute prob of change with pad_baseC = - sd_pad
    # if dvs is gm_icv, we need to adjust the values
    vec_values[2] <- - sd(data_model$pad_baseC)

    predict_and_ci_minus_sd <- predict_and_ci95(vec_coeffs, c_wglm_vcov, vec_values)

    # print the results
    print("Predictions and confidence intervals")
    print(predict_and_ci_minus_sd)

    # get odds
    results_wglm <- model.tables(c_wglm)
    esimate_wglm <- results_wglm[2,c("estimate", "lower", "upper", "p.value")]
    odds <- exp(esimate_wglm[,])
    odds["p.value"] <- esimate_wglm["p.value"]

    # store in table total_result
    total_result <- rbind(
        total_result,
        data.table(
            model = dvs[i],
            odds = odds[,"estimate"],
            odds_lower = odds[,"lower"],
            odds_upper = odds[,"upper"],
            odds_p_value = odds[, "p.value"],
            prob_pad_minus_sd = predict_and_ci_minus_sd["prob"],
            ci_lower_pad_minus_sd = predict_and_ci_minus_sd["ci_lower"],
            ci_upper_pad_minus_sd = predict_and_ci_minus_sd["ci_upper"],
            prob_pad_mu = predict_and_ci_mu["prob"],
            ci_lower_pad_mu = predict_and_ci_mu["ci_lower"],
            ci_upper_pad_mu = predict_and_ci_mu["ci_upper"],
            prob_pad_plus_sd = predict_and_ci_plus_sd["prob"],
            ci_lower_pad_plus_sd = predict_and_ci_plus_sd["ci_lower"],
            ci_upper_pad_plus_sd = predict_and_ci_plus_sd["ci_upper"],
            auc = score$AUC$score[[2, "AUC"]],
            auc_lower = score$AUC$score[[2, "lower"]],
            auc_upper = score$AUC$score[[2, "upper"]],
            brier = score$Brier$score[[3, "Brier"]],
            brier_lower = score$Brier$score[[3, "lower"]],
            brier_upper = score$Brier$score[[3, "upper"]]
        )
    )

    # compute ci and prob for all values in pad_baseC_seq for plotting purposes
    # this should be transformed into a matrix multiplication. nobody likes ugly for loops

    wiggle = (max(data_model$pad_baseC) - min(data_model$pad_baseC)) / 10

    pad_baseC_seq <- seq(min(data_model$pad_baseC) - wiggle, max(data_model$pad_baseC) + wiggle, length.out = 100)
    pred_ci <- matrix(NA, nrow = length(pad_baseC_seq), ncol = 3)
    for (j in seq_along(pad_baseC_seq)) {
        vec_values[2] <- pad_baseC_seq[j]
        pred_ci[j, ] <- predict_and_ci95(vec_coeffs, c_wglm_vcov, vec_values)
    }

    # provide column names with currend dv name e.g.
    colnames(pred_ci) <- c(paste0(dvs[i], "_pred"), paste0(dvs[i], "_ci_lower"), paste0(dvs[i], "_ci_upper"))

    # merge with total preds
    total_preds <- cbind(total_preds, pred_ci)
}

# Combine stats
# remove _base from model names
total_result$model <- gsub("_base", "", total_result$model)

# write to excel
store_results_path <- here::here("results", "long_c_results.xlsx")
writexl::write_xlsx(total_result, path = store_results_path)

# write total_preds to excel
store_preds_path <- here::here("results", "long_c_preds.xlsx")
writexl::write_xlsx(total_preds, path = store_preds_path)


# plot roc for all but first make stupid baseline
glm_base <- glm(diagnosis_change ~ age_baseC + aes_baseC + gender,
    family = binomial(link = "logit"),
    data = data_model
)

score <- Score(
    list(glm = glm_base, 
         brainageR = wglm_models$pad_brainageR_base,
         DeepBrainNet = wglm_models$pad_DeepBrainNet_base,
         brainage = wglm_models$pad_brainage_base,
         enigma = wglm_models$pad_enigma_base,
         pyment = wglm_models$pad_pyment_base,
         mccqrnn = wglm_models$pad_mccqrnn_base,
         gm_icv = wglm_models$gm_icv_base),
    times = 4,
    formula = Surv(survive_timing, diagnosis_change) ~ age_baseC + aes_base + gender,
    data = data_model,
    metrics = c("AUC", "Brier"),
    plot = c("ROC")
)

plotROC(score)

# store data
roc_data <- score$ROC$plotframe
store_roc_path <- here::here("results", "long_c_roc_data.xlsx")
writexl::write_xlsx(roc_data, path = store_roc_path)
