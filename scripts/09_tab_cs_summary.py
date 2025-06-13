# %% ### import modules
from pathlib import Path
import pandas as pd

# %% read results table
results_dir = Path(__file__).resolve().parents[1] / "results"
df = pd.read_excel(results_dir / "cs_summary.xlsx")

# %%
# Define the new column names
metrics = ["pad_brainageR", "pad_DeepBrainNet", "pad_brainage", "pad_enigma", "pad_pyment", "pad_mccqrnn", "gm_icv"]

df_mean = df.copy()
# Reformat columns
for metric in metrics:
    df_mean[metric] = df_mean.apply(
        lambda row: f"{row[f'{metric}_mean']} ({row[f'{metric}_lower_ci']}, {row[f'{metric}_upper_ci']})",
        axis=1,
    )

# Keep only the formatted columns
df_mean = df_mean[["diagnosis"] + metrics]

# remove pad_ from column names
df_mean.columns = ["diagnosis"] + [col.replace("pad_", "") for col in df_mean.columns if col != "diagnosis"]
df_mean.rename(columns={"gm_icv": "Gray Matter"}, inplace=True)

df_mean = df_mean.set_index("diagnosis").T

# %% now for variance
df_var = df.copy()
for metric in metrics:
    df_var[metric] = df_var.apply(
        lambda row: f"{row[f'{metric}_sd']}",
        axis=1,
    )

# Keep only the formatted columns
df_var = df_var[["diagnosis"] + metrics]

# remove pad_ from column names
df_var.columns = ["diagnosis"] + [col.replace("pad_", "") for col in df_var.columns if col != "diagnosis"]
df_var.rename(columns={"gm_icv": "Gray Matter"}, inplace=True)

df_var = df_var.set_index("diagnosis").T
# %% combine
df_combined = pd.concat([df_mean, df_var], axis=1)
df_combined.columns = pd.MultiIndex.from_product([["Mean (95%CI)", "SD"], ["CN", "MCI", "AD"]])

# %%
print(df_combined)

# %%
df_combined.to_excel(results_dir / "cs_summary_formatted.xlsx")

# %%
