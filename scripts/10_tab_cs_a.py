# %% ### import modules
from pathlib import Path
import seaborn as sns
import pandas as pd

from matplotlib import pyplot as plt

import bap1b.figures as figures

# %% read results table
results_dir = Path(__file__).resolve().parents[1] / "results"

df = pd.read_excel(results_dir / "cs_a_results.xlsx")

# %% round all numeric values to 2 decimals columns ["Estimate", "SE", "effect_size", "lower_ci", "upper_ci"]
df[["Estimate", "Std. Error", "effect_size", "lower_ci", "upper_ci"]] = df[
    ["Estimate", "Std. Error", "effect_size", "lower_ci", "upper_ci"]
].round(2)

# %% combine effect_size, lower_ci, upper_ci to effect_size (lower_ci, upper_ci)
df["Cohens d (95%CI)"] = (
    df["effect_size"].astype(str) + " (" + df["lower_ci"].astype(str) + ", " + df["upper_ci"].astype(str) + ")"
)

# %% combine Estimate and SE to Estimate (SE)
df["Beta (SE)"] = df["Estimate"].astype(str) + " (" + df["Std. Error"].astype(str) + ")"

# %% if the p-value is below 0.01 write < 0.01, otherwise write the value
df["p-value"] = df["Pr(>|t|)"].apply(lambda x: "<0.001" if x < 0.001 else round(x, 3))

# %% in comparison replace _ with vs
df["Comparison"] = df["comparison"].str.replace("_", " vs ")

# %% in model remove pad_
df["Model"] = df["model"].str.replace("pad_", "")

# %% in model replace gm_icv with GM/ICV
df["Model"] = df["Model"].str.replace("gm_icv", "GM/ICV")

# %% pivot the table based on comparison
df = df.pivot(index="Model", columns="Comparison", values=["Beta (SE)", "Cohens d (95%CI)", "p-value"])
df.columns = df.columns.swaplevel(0, 1)
df = df.sort_index(axis=1, level=0)
df = df.reindex(["brainageR", "DeepBrainNet", "brainage", "enigma", "pyment", "mccqrnn", "GM/ICV"])
df

# %% sort the table based on the list: brainageR, DeepBrainNet, brainage, enigma, pyment, mccqrnn

# %% save table
df.to_excel(results_dir / "cs_a_results_formatted.xlsx", index=True)
