# %% ### import modules
from pathlib import Path
import pandas as pd

# %% read results table
results_dir = Path(__file__).resolve().parents[1] / "results"

# %% data
df = pd.read_excel(results_dir / "long_fg_results.xlsx")

# %% round values
df["estimate"] = df["estimate"].round(2)
df["lower"] = df["lower"].round(2)
df["upper"] = df["upper"].round(2)
df["p.value"] = df["p.value"].round(3)

# %% combine columns (we need the a in ar for sorting. stupid, but it is ok for now)
df["ar (95%)"] = df["estimate"].astype(str) + " (" + df["lower"].astype(str) + "," + df["upper"].astype(str) + ")"

# rename model to Model
df = df.rename(columns={"model": "Model", "p.value": "p-value"})

# remove pad_
df["Model"] = df["Model"].str.replace("pad_", "")

# rename gm_icv to GM/ICV in target
df["target"] = df["target"].str.replace("gm_icv", "GM/ICV")

# replace p-value < 0.001 with <0.001
df["p-value"] = df["p-value"].apply(lambda x: "<0.001" if x < 0.001 else x)

# %% select
# only keep straggy = 2
df = df[df["strategy"] == 2]
# df = df[df["diagnosis"] == "CN"]

# %% make table wide based
df = df.pivot(index=["Model", "diagnosis"], columns="target", values=["ar (95%)", "p-value"])

# %%
df.columns = df.columns.swaplevel(0, 1)
df = df.sort_index(axis=1, level=0)

# %% make it even wider by target
df = df.unstack(level=1)
# %%
df.columns = df.columns.swaplevel(1, 2)
# %%
df = df.sort_index(axis=1, level=0)
# %%
# rename column from ar to r (95%)
df.columns = df.columns.set_levels(["r (95%)", "p-value"], level=2)

# remove pad_ from model (index)
models = ["brainageR", "DeepBrainNet", "brainage", "enigma", "pyment", "mccqrnn", "gm_icv"]

# sort MultiIndex by models
df = df.reindex(models, level=0)

# %% save
df.to_excel(results_dir / "long_fg_results_formatted.xlsx", index=True)

# %%
