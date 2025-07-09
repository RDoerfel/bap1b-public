### Summary stats for the cross sectional analysis
# Load required libraries
library(readxl)
library(writexl)
library(here)
library(dplyr)

# Load excel file
file_path <- here::here("data", "megamastersheet_simulated.xlsx")
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

#### Prepare data for Table 1
# Table with 4 rows: CN, MCI, AD, Total. Columns: N M000, N M012, N M024, N M036, N M048, Age mean (SD), Sex [F/M], Years of EDU mean (SD), ADNI_MEM mean (SD), AES mean (SD)

# Create a summary table
summary_table <- data %>%
    dplyr::group_by(diagnosis) %>%
    dplyr::summarise(
        M000 = sum(session_id == "M000"),
        M012 = sum(session_id == "M012"),
        M024 = sum(session_id == "M024"),
        M036 = sum(session_id == "M036"),
        M048 = sum(session_id == "M048"),
        Age_mean_sd = paste0(round(mean(chron_age[session_id == "M000"], na.rm = TRUE), 2), " (", round(sd(chron_age[session_id == "M000"], na.rm = TRUE), 2), ")"),
        Sex_FM = paste0(sum(gender[session_id == "M000"] == "Female"), "/", sum(gender[session_id == "M000"] == "Male")),
        EDU_mean_sd = paste0(round(mean(education_level[session_id == "M000"], na.rm = TRUE), 2), " (", round(sd(education_level[session_id == "M000"], na.rm = TRUE), 2), ")"),
        ADNI_MEM_mean_sd = paste0(round(mean(ADNI_MEM[session_id == "M000"], na.rm = TRUE), 2), " (", round(sd(ADNI_MEM[session_id == "M000"], na.rm = TRUE), 2), ")"),
        AES_mean_sd = paste0(round(mean(aes[session_id == "M000"], na.rm = TRUE), 2), " (", round(sd(aes[session_id == "M000"], na.rm = TRUE), 2), ")")
    )

# Add a total row
total_row <- data %>%
    dplyr::summarise(
        diagnosis = "Total",
        M000 = sum(session_id == "M000"),
        M012 = sum(session_id == "M012"),
        M024 = sum(session_id == "M024"),
        M036 = sum(session_id == "M036"),
        M048 = sum(session_id == "M048"),
        Age_mean_sd = paste0(round(mean(chron_age[session_id == "M000"], na.rm = TRUE), 2), " (", round(sd(chron_age[session_id == "M000"], na.rm = TRUE), 2), ")"),
        Sex_FM = paste0(sum(gender[session_id == "M000"] == "Female"), "/", sum(gender[session_id == "M000"] == "Male")),
        EDU_mean_sd = paste0(round(mean(education_level[session_id == "M000"], na.rm = TRUE), 2), " (", round(sd(education_level[session_id == "M000"], na.rm = TRUE), 2), ")"),
        ADNI_MEM_mean_sd = paste0(round(mean(ADNI_MEM[session_id == "M000"], na.rm = TRUE), 2), " (", round(sd(ADNI_MEM[session_id == "M000"], na.rm = TRUE), 2), ")"),
        AES_mean_sd = paste0(round(mean(aes[session_id == "M000"], na.rm = TRUE), 2), " (", round(sd(aes[session_id == "M000"], na.rm = TRUE), 2), ")")
    )

# Combine the summary table and the total row
summary_table <- dplyr::bind_rows(summary_table, total_row)

# Print the summary table
print(summary_table)

# Save summary table
summary_table_path <- here::here("results", "summary_table.xlsx")
writexl::write_xlsx(summary_table, path = summary_table_path)