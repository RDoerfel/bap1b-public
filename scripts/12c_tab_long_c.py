# %% ### import modules
from pathlib import Path
import pandas as pd

# %% read results table
results_dir = Path(__file__).resolve().parents[1] / "results"

# %% data
df = pd.read_excel(results_dir / "long_c_results_cox.xlsx")

# %% round intercept and coef to 3 decimal places. same for intercept_ste and coef_ste
df["hazard"] = df["hazard"].round(2)
df["hazard_lower"] = df["hazard_lower"].round(2)
df["hazard_upper"] = df["hazard_upper"].round(2)
df["survival_minus_sd"] = df["survival_minus_sd"].round(2)
df["ci_lower_minus_sd"] = df["ci_lower_minus_sd"].round(2)
df["ci_upper_minus_sd"] = df["ci_upper_minus_sd"].round(2)
df["survival_mu"] = df["survival_mu"].round(2)
df["ci_lower_mu"] = df["ci_lower_mu"].round(2)
df["ci_upper_mu"] = df["ci_upper_mu"].round(2)
df["survival_plus_sd"] = df["survival_plus_sd"].round(2)
df["ci_lower_plus_sd"] = df["ci_lower_plus_sd"].round(2)
df["ci_upper_plus_sd"] = df["ci_upper_plus_sd"].round(2)
df["auc"] = df["auc"].round(2)
df["auc_lower"] = df["auc_lower"].round(2)
df["auc_upper"] = df["auc_upper"].round(2)
df["brier"] = df["brier"].round(2)
df["brier_lower"] = df["brier_lower"].round(2)
df["brier_upper"] = df["brier_upper"].round(2)

# %% combine hazard and hazard_lower and hazard_upper to hazard (hazard_lower, hazard_upper)
df["Hazard (95%CI)"] = (
    df["hazard"].astype(str) + " (" + df["hazard_lower"].astype(str) + "," + df["hazard_upper"].astype(str) + ")"
)

# %% if the p-value is below 0.01 write < 0.01, otherwise write the value
df["p-value (Hazard)"] = df["hazard_p_value"].apply(lambda x: "<0.001" if x < 0.001 else round(x, 3))

# %% in model remove pad_
df["Model"] = df["model"].str.replace("pad_", "")

# %% combine probpad_minus5 and cipad_minus5 to P(PAD=-5) (95%CI)
df["P(x=-SD) (95%CI)"] = (
    df["survival_minus_sd"].astype(str)
    + " ("
    + df["ci_lower_minus_sd"].astype(str)
    + ","
    + df["ci_upper_minus_sd"].astype(str)
    + ")"
)

# %% combine survival0 and ci_pad0 to P(PAD=0) (95%CI)
df["P(x=mu) (95%CI)"] = (
    df["survival_mu"].astype(str) + " (" + df["ci_lower_mu"].astype(str) + "," + df["ci_upper_mu"].astype(str) + ")"
)

# %% combine survival5 and ci_pad5 to P(PAD=5) (95%CI)
df["P(x=SD) (95%CI)"] = (
    df["survival_plus_sd"].astype(str)
    + " ("
    + df["ci_lower_plus_sd"].astype(str)
    + ","
    + df["ci_upper_plus_sd"].astype(str)
    + ")"
)

# %% AUC
df["AUC (95%CI)"] = df["auc"].astype(str) + " (" + df["auc_lower"].astype(str) + "," + df["auc_upper"].astype(str) + ")"

# %% Brier
df["Brier (95%CI)"] = (
    df["brier"].astype(str) + " (" + df["brier_lower"].astype(str) + "," + df["brier_upper"].astype(str) + ")"
)

# %% in model replace gm_icv with GM/ICV
df["Model"] = df["Model"].str.replace("gm_icv", "GM/ICV")

# %% save table
df_to_save = df[
    [
        "Model",
        "Hazard (95%CI)",
        "p-value (Hazard)",
        "P(x=-SD) (95%CI)",
        "P(x=mu) (95%CI)",
        "P(x=SD) (95%CI)",
        "AUC (95%CI)",
        "Brier (95%CI)",
    ]
]
df_to_save.to_excel(results_dir / "long_c_results_formatted_cox.xlsx", index=False)
# %%
