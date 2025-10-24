# Analysis c) to investigate whether an increased PAD at baseline is associated with future conversion from CN/MCI to MCI/AD.
# Load required libraries
library(mvtnorm)
library(data.table)
library(LMMstar)
library(mets)
library(riskRegression)
library(dplyr)
library(lava)

### Script

# Load excel file
file_path <- here::here("data", "megamastersheet.xlsx")
data <- readxl::read_excel(file_path)
print(dim(data))
# [1] 10802    42

# only take data that is indicated by column keep
data_keep <- data[data$keep == 1, ]
print(dim(data_keep))
# [1] 9768   42

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
    time_from_baseline,
    hippocampus_icv
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
    dplyr::filter(dplyr::first(session_id) == "M000" & dplyr::first(diagnosis) == "CN") %>%
    dplyr::ungroup()

print("Number of scans after filtering for CN at baseline:")
print(dim(data_subset))
# [1] 2600   15

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
# [1] 1083   28

# compute "prevelance" and number of subjects with change
print("Prevelance of change from CN to MCI/AD")
print(mean(data_subset_base_survive$diagnosis_change))
# [1] 0.2539243

print("Number of subjects with change from CN to MCI/AD")
print(sum(data_subset_base_survive$diagnosis_change))
# [1] 275

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
store_path <- here::here("data", "data_long_c_cox.xlsx")
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
    hazard = numeric(),
    hazard_lower = numeric(),
    hazard_upper = numeric(),
    hazard_p_value = numeric(),
    survival_minus_sd = numeric(),
    ci_lower_minus_sd = numeric(),
    ci_upper_minus_sd = numeric(),
    survival_mu = numeric(),
    ci_lower_mu = numeric(),
    ci_upper_mu = numeric(),
    survival_plus_sd = numeric(),
    ci_lower_plus_sd = numeric(),
    ci_upper_plus_sd = numeric(),
    auc = numeric(),
    auc_lower = numeric(),
    auc_upper = numeric(),
    brier = numeric(),
    brier_lower = numeric(),
    brier_upper = numeric()
)

total_preds <- data.table()

cox_models <- list()

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

    ## Prevelance of change from CN to MCI/AD
    print("Actual analysis accounting for censoring:")
    print(mean(data_model$diagnosis_change))
    # [1] 0.07411631

    ## censoring cox
    ## censoring wglm (cox)
    e_cox <- coxph(Surv(survive_timing, diagnosis_change) ~ pad_baseC + age_baseC + aes_baseC + gender,
        data = data_model,
        x = TRUE
    )

    # add to list
    cox_models[[dvs[i]]] <- e_cox


    print(summary(e_cox)$coef)
    print(exp(cumsum(coef(e_cox))))

    # AUC
    print("AUC")
    score <- Score(
        list(cox = e_cox),
        times = 4,
        formula = Surv(survive_timing, diagnosis_change) ~ pad_baseC + age_baseC + aes_baseC + gender,
        data = data_model,
        metrics = c("AUC", "Brier"),
        plot = c("ROC")
    )

    print(score$AUC)
    print(score$Brier)
    plotROC(score)

    # compute survial probabilities at time 4 for mean, mean - sd, mean + sd of pad_baseC
    cox_probs <- predictCox(e_cox, times = c(4), newdata = data.frame(gender="Female",
        pad_baseC=c(mean(data_model$pad_baseC) - sd(data_model$pad_baseC), mean(data_model$pad_baseC), mean(data_model$pad_baseC) + sd(data_model$pad_baseC)),
        age_baseC=mean(data_model$age_baseC),
        aes_baseC=mean(data_model$aes_baseC)),
        product.limit = TRUE,
        se = TRUE)

    # get hazard ratios
    exp_coeffs <- exp(coef(e_cox))
    exp_confint <- exp(confint(e_cox))
    p_values <- summary(e_cox)$coefficients[, "Pr(>|z|)"]

    hazard <- exp_coeffs['pad_baseC']
    hazard_lower <- exp_confint['pad_baseC', 1]
    hazard_upper <- exp_confint['pad_baseC', 2]
    hazard_p <- p_values['pad_baseC']

    # store in table total_result
    total_result <- rbind(
        total_result,
        data.table(
            model = dvs[i],
            hazard = hazard,
            hazard_lower = hazard_lower,
            hazard_upper = hazard_upper,
            hazard_p_value = hazard_p,
            survival_minus_sd = cox_probs$survival[1, ],
            ci_lower_minus_sd = cox_probs$survival.lower[1, ],
            ci_upper_minus_sd = cox_probs$survival.upper[1, ],
            survival_mu = cox_probs$survival[2, ],
            ci_lower_mu = cox_probs$survival.lower[2, ],
            ci_upper_mu = cox_probs$survival.upper[2, ],
            survival_plus_sd = cox_probs$survival[3, ],
            ci_lower_plus_sd = cox_probs$survival.lower[3, ],
            ci_upper_plus_sd = cox_probs$survival.upper[3, ],
            auc = score$AUC$score[[1, "AUC"]],
            auc_lower = score$AUC$score[[1, "lower"]],
            auc_upper = score$AUC$score[[1, "upper"]],
            brier = score$Brier$score[[2, "Brier"]],
            brier_lower = score$Brier$score[[2, "lower"]],
            brier_upper = score$Brier$score[[2, "upper"]]
        )
    )
}

# Combine stats
# remove _base from model names
total_result$model <- gsub("_base", "", total_result$model)

# write to excel
store_results_path <- here::here("results", "long_c_results_cox.xlsx")
writexl::write_xlsx(total_result, path = store_results_path)

# plot roc for all but first make stupid baseline
glm_base <- glm(diagnosis_change ~ age_baseC + aes_baseC + gender,
    family = binomial(link = "logit"),
    data = data_model
)

score <- Score(
    list(glm = glm_base, 
         brainageR = cox_models$pad_brainageR_base,
         DeepBrainNet = cox_models$pad_DeepBrainNet_base,
         brainage = cox_models$pad_brainage_base,
         enigma = cox_models$pad_enigma_base,
         pyment = cox_models$pad_pyment_base,
         mccqrnn = cox_models$pad_mccqrnn_base,
         gm_icv = cox_models$gm_icv_base),
    times = 4,
    formula = Surv(survive_timing, diagnosis_change) ~ age_baseC + aes_base + gender,
    data = data_model,
    metrics = c("AUC", "Brier"),
    plot = c("ROC")
)

plotROC(score)

# store data
roc_data <- score$ROC$plotframe
store_roc_path <- here::here("results", "long_c_roc_data_cox.xlsx")
writexl::write_xlsx(roc_data, path = store_roc_path)
