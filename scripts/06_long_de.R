# Analyses d) and e): to see whether an increased PAD at baseline is associated with a future increase in atrophy (d) and cognition (e)

# Load required libraries
library(mvtnorm)
library(data.table)
library(LMMstar)
library(mets)
library(riskRegression)
library(dplyr)
library(boot)
library(writexl)
library(mmrm)

### Functions
get_times <- function(sessions) {
    time <- as.numeric(gsub("M", "", levels(sessions), fixed = TRUE))
    return(time)
}

get_contrast <- function(time, do_special) {
    # do_special is only true if we look at the same measure at baseline vs cahnge from baseline.
    # e.g. gm_v vs change in gm_v
    n_times <- length(time)
    weight <- (time - mean(time)) / sum((time - mean(time))^2) # (time is centered and scaled)

    if(do_special) {
        contrast <- rbind(c(1, rep(0, n_times-1)), weight)
    } else {
        contrast <- rbind(c(1, rep(0, n_times)), c(0, weight))
    }
    return(contrast)
}

select_outcome_vatriable <- function(data, outcome_variable) {
    data_subset <- data %>%
        dplyr::select("subject_id", "session_id", "time_from_baseline", "pad_base", all_of(outcome_variable))
    data_renamed <- data_subset %>% dplyr::rename(target = all_of(outcome_variable))
    data_table <- setDT(data_renamed)
    return(data_table)
}

strategy_1 <- function(data_long) {
    slopes <- data_long %>%
        dplyr::group_by(subject_id) %>%
        dplyr::arrange(time_from_baseline) %>%
        dplyr::summarise(
            slope_cogn = coef(lm(target ~ time_from_baseline))["time_from_baseline"],
            PAD = pad_base[1]
        )

    cor <- cor.test(slopes$PAD, slopes$slope_cogn)
    table_cor <- data.table(estimate = cor$estimate, lower = cor$conf.int[1], upper = cor$conf.int[2], p.value = cor$p.value)

    return(table_cor)
}

prepare_data_for_strat2 <- function (data_model_target) {
    data_semi_long <- melt(data_model_target,
                            id.vars = c("subject_id", "session_id", "time_from_baseline"),
                            measure.vars = c("pad_base", "target"),
                            variable.name = "variable",
                            value.name = "value")

    data_long <- data_semi_long %>% dplyr::filter(!(variable == "pad_base" & time_from_baseline != 0.0))

    # add session indicator to variable name (e.g., ADNI_MEM_M000)
    data_long$variable <- paste(data_long$variable, data_long$session_id, sep = "_")
    data_long$variable <- factor(data_long$variable)
    data_long$variable <- relevel(data_long$variable, "pad_base_M000")

    return(data_long)
}

strategy_2 <- function(data_long, contrast) {
    unstructured_lmm <- lmm(value ~ variable, repetition = ~ variable | subject_id, structure = "UN", data = data_long)
    cor <- lava::estimate(unstructured_lmm, function(p){
        sigma_vy <- sigma(unstructured_lmm, p = p)
        sigma_vb <- contrast %*% sigma_vy %*% t(contrast)
        cov2cor(sigma_vb)[1, 2]
    })
    # only keep coluns estimate, lower, upper, p.value
    cor <- cor[, c("estimate", "lower", "upper", "p.value")]
    return(cor)
}

strategy_2_mmrm_lmm <- function(data_longlong, contrast) {
    # first, estimate the residual variance-covariance matrix using MMRM (has no convergence issues)
    strat2_mmrm <- mmrm::mmrm(value ~ variable + us(variable | subject_id), data = data_longlong)
    sigma_mmrm <- VarCorr(strat2_mmrm)
    # second, use the estimated residual variance-covariance matrix as initial values for the LMM
    strat2_lmm_mmrm <- lmm(
        value ~ variable,
        repetition = ~ variable | subject_id,
        structure = "UN",
        data = data_longlong,
        control = list(init = sigma_mmrm)
    )
    # we do this to use the delta method to get confidence intervals for the estimate.
    cor_lmm_mmrm <- lava::estimate(strat2_lmm_mmrm, function(p) {
        sigma_vy <- sigma(strat2_lmm_mmrm, p = p)
        sigma_vb <- contrast %*% sigma_vy %*% t(contrast)
        cov2cor(sigma_vb)[1, 2]
    })
    cor_lmm_mmrm <- cor_lmm_mmrm[, c("estimate", "lower", "upper", "p.value")]
    return(cor_lmm_mmrm)
}

prepare_data_for_strat3 <- function(data_long) {
    # make it wide for strategy 3
    data_long_no_time <- data_long %>% dplyr::select(-time_from_baseline)
    data_wide <- dcast(data_long_no_time,
                        subject_id ~ variable,
                        value.var = "value")
    return(data_wide)
}

get_sigma_vb_strat_3 <- function(data, contrast) {
    sigma_vy <- cov(data, use = "pairwise.complete.obs")
    sigma_vb <- contrast %*% sigma_vy %*% t(contrast)
    return(sigma_vb)
}

bs_strat3 <- function(data, indices, contrast) {
    data_boot <- data[indices, ]
    sigma_vb_strat3 <- get_sigma_vb_strat_3(data_boot, contrast)
    cor <- cov2cor(sigma_vb_strat3)[1, 2]
    return(cor)
}

strategy_3 <- function(data, contrast) {
    bootstrap <- boot(data = data, stype = "i", statistic = bs_strat3, R = 1000, contrast = contrast)
    sd_bs <- sd(bootstrap$t)
    mean_bs = mean(bootstrap$t)
    ci_l <- bootstrap$t0 - 1.96 * sd_bs
    ci_u <- bootstrap$t0 + 1.96 * sd_bs
    # weird thing: mean_bs is not the same as what we get from summary(bootstrap)
    original <- bootstrap$t0
    results <- cbind(original, ci_l, ci_u, NA)
    # change names of columns
    colnames(results) <- c("estimate", "lower", "upper", "p.value")
    # make it a data.table
    results <- data.table(results)
    plot(bootstrap)
    return(results)
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
    starts_with("pad_"),
    gm_icv,
    chron_age,
    time_from_baseline,
    ADNI_MEM
)

# only until 4 years
data_4years <- data_subset %>% dplyr::filter(time_from_baseline <= 4)
print("Number of scans until 4 years:")
print(dim(data_4years))
# [1] 8067   16

# Only keep yearly visits
data_4years_yearly <- data_4years[data_4years$session_id %in% c("M006", "M018", "M030", "M042") == FALSE, ]
data_4years_yearly <- droplevels(data_4years_yearly)
print(table(data_4years_yearly$session_id))

# Only keep subjects with more than 2 visits
data_4years_yearly_min2 <- data_4years_yearly %>% dplyr::group_by(subject_id) %>% dplyr::filter(n() > 2) %>% dplyr::ungroup()
print("Number of scans until with more than 2 visits:")
print(dim(data_4years_yearly_min2))

# store results in list
results <- list()

# iterate through all diagnosis (for now, only CN and MCI)
for (diag in c("CN", "MCI", "AD")) {
    # only CN in first scan TODO: only check for first scan!
    data_diagnosis <- data_4years_yearly_min2 %>%
        dplyr::arrange(subject_id, time_from_baseline) %>%
        dplyr::group_by(subject_id) %>%
        dplyr::filter(dplyr::first(diagnosis) == diag) %>%
        dplyr::ungroup()

    # if diag == "AD", drop scans M036 and M048
    if (diag == "AD") {
        data_diagnosis <- data_diagnosis %>% dplyr::filter(session_id != "M036" & session_id != "M048")
        data_diagnosis <- droplevels(data_diagnosis)
    }
    
    print(sprintf("Number of %s scans until 4 years:", diag))
    print(dim(data_diagnosis))
    print(table(data_diagnosis$session_id))

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

    # this is dangerous! expects that the dvs are actually in that order in the inital dataframe
    idx_dvs <- which(colnames(data_diagnosis) %in% dvs)

    for (i in 1:length(dvs)) {

        # get index of current model
        idx <- idx_dvs[i]

        # rename the column of this index so that we can run the same lines of code for all models
        data_model <- data_diagnosis %>%
            dplyr::select(subject_id,
                          session_id,
                          time_from_baseline,
                          pad_model = all_of(idx),
                          diagnosis,
                          chron_age,
                          aes,
                          gender,
                          gm_icv,
                          ADNI_MEM,
                          gm_icv = gm_icv)

        # add pad at baseline to data
        data_model <- data_model %>%
            dplyr::group_by(subject_id) %>%
            arrange(subject_id, session_id) %>%
            mutate(pad_base = first(pad_model)) %>%
            ungroup()

        # iterate through target variables ADNI_MEM and gm_icv
        results_target <- list()
        targets <- c("ADNI_MEM", "gm_icv")
        for (target in targets) {
            message(sprintf("Diagnosis %s, model %s, target %s", diag, dvs[i], target))

            # get times
            time <- get_times(data_model$session_id)
            do_special <- (target == "gm_icv") && (dvs[i] == "gm_icv")
            contrast <- get_contrast(time, do_special)

            # subset to only actually used data
            # prepare data to be in formats that are needed for each strategy
            data_model_target <- select_outcome_vatriable(data_model, target)

            # prepare data for strategy 2
            data_long_target <- prepare_data_for_strat2(data_model_target)

            # if we look at gm_icv, we need to remove the pad_base_M000 entry, as it is the same as target_M000
            if (do_special) {
                data_long_target <- data_long_target %>% dplyr::filter(!(variable == "pad_base_M000"))
                data_long_target <- droplevels(data_long_target)
            }

            # prepare data for strategy 3
            data_wide_target <- prepare_data_for_strat3(data_long_target)

            ## Strategy 1: Compute slopes per subject.
            strat1 <- strategy_1(data_model_target)
            print(strat1)

            ## Strdategy 2: estimate covariance matrix using lmm
            strat2 <- strategy_2(data_long_target, contrast)
            print(strat2)

            ## Strategy 3: estimate covariance matrix using cov
            # run strategy 3
            strat3 <- strategy_3(data_wide_target[, 2 : (1 + dim(contrast)[2])], contrast)
            print(strat3)

            # add results to results_all list
            strat1$diagnosis <- diag
            strat1$target <- target
            strat1$strategy <- "1"
            strat1$model <- dvs[i]
            strat2$diagnosis <- diag
            strat2$target <- target
            strat2$strategy <- "2"
            strat2$model <- dvs[i]
            strat3$diagnosis <- diag
            strat3$target <- target
            strat3$strategy <- "3"
            strat3$model <- dvs[i]
            strategies_results_target <- list(strat1, strat2, strat3)
            results_target[[length(results_target) + 1]] <- rbindlist(strategies_results_target)
        }
        results[[length(results) + 1]] <- rbindlist(results_target)
    }
}
# combine all results
results_table <- rbindlist(results)

# store results in excel
store_results_path <- here::here("results", "long_de_results.xlsx")
writexl::write_xlsx(results_table, store_results_path)

