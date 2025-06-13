# %% ### import modules
from pathlib import Path
import pandas as pd

# %% read results table
results_dir = Path(__file__).resolve().parents[1] / "results"

# %% data
df = pd.read_excel(results_dir / "long_c_results.xlsx")

# %% round intercept and coef to 3 decimal places. same for intercept_ste and coef_ste
df["odds"] = df["odds"].round(2)
df["odds_lower"] = df["odds_lower"].round(2)
df["odds_upper"] = df["odds_upper"].round(2)
df["odds_p_value"] = df["odds_p_value"].round(3)
df["prob_pad_minus_sd"] = df["prob_pad_minus_sd"].round(2)
df["ci_lower_pad_minus_sd"] = df["ci_lower_pad_minus_sd"].round(2)
df["ci_upper_pad_minus_sd"] = df["ci_upper_pad_minus_sd"].round(2)
df["prob_pad_mu"] = df["prob_pad_mu"].round(2)
df["ci_lower_pad_mu"] = df["ci_lower_pad_mu"].round(2)
df["ci_upper_pad_mu"] = df["ci_upper_pad_mu"].round(2)
df["prob_pad_plus_sd"] = df["prob_pad_plus_sd"].round(2)
df["ci_lower_pad_plus_sd"] = df["ci_lower_pad_plus_sd"].round(2)
df["ci_upper_pad_plus_sd"] = df["ci_upper_pad_plus_sd"].round(2)
df["auc"] = df["auc"].round(2)
df["auc_lower"] = df["auc_lower"].round(2)
df["auc_upper"] = df["auc_upper"].round(2)
df["brier"] = df["brier"].round(2)
df["brier_lower"] = df["brier_lower"].round(2)
df["brier_upper"] = df["brier_upper"].round(2)

# %% combine odds and odds_lower and odds_upper to odds (odds_lower, odds_upper)
df["Odds (95%CI)"] = (
    df["odds"].astype(str) + " (" + df["odds_lower"].astype(str) + "," + df["odds_upper"].astype(str) + ")"
)

# %% if the p-value is below 0.01 write < 0.01, otherwise write the value
df["p-value (Odds)"] = df["odds_p_value"].apply(lambda x: "<0.001" if x < 0.001 else round(x, 3))

# %% in model remove pad_
df["Model"] = df["model"].str.replace("pad_", "")

# %% combine probpad_minus5 and cipad_minus5 to P(PAD=-5) (95%CI)
df["P(x=-SD) (95%CI)"] = (
    df["prob_pad_minus_sd"].astype(str)
    + " ("
    + df["ci_lower_pad_minus_sd"].astype(str)
    + ","
    + df["ci_upper_pad_minus_sd"].astype(str)
    + ")"
)

# %% combine prob_pad0 and ci_pad0 to P(PAD=0) (95%CI)
df["P(x=mu) (95%CI)"] = (
    df["prob_pad_mu"].astype(str)
    + " ("
    + df["ci_lower_pad_mu"].astype(str)
    + ","
    + df["ci_upper_pad_mu"].astype(str)
    + ")"
)

# %% combine prob_pad5 and ci_pad5 to P(PAD=5) (95%CI)
df["P(x=SD) (95%CI)"] = (
    df["prob_pad_plus_sd"].astype(str)
    + " ("
    + df["ci_lower_pad_plus_sd"].astype(str)
    + ","
    + df["ci_upper_pad_plus_sd"].astype(str)
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
        "Odds (95%CI)",
        "p-value (Odds)",
        "P(x=-SD) (95%CI)",
        "P(x=mu) (95%CI)",
        "P(x=SD) (95%CI)",
        "AUC (95%CI)",
        "Brier (95%CI)",
    ]
]
df_to_save.to_excel(results_dir / "long_c_results_formatted.xlsx", index=False)
# %%
