# %%
from pathlib import Path
import pandas as pd

# %%
df = pd.read_excel(Path(__file__).parent.parent / "data" / "megamastersheet.xlsx")
df_fs = pd.read_csv(Path(__file__).parent.parent / "data" / "LandRVolumes.csv")

# %%
mm3_to_liter = 1e-6
df_fs["hippocampus"] = (df_fs["Lhippo"] + df_fs["Rhippo"]) * mm3_to_liter

# %% rename SubjID to mr_id_x
df_fs = df_fs.rename(columns={"SubjID": "mr_id_x"})

# %%
df = df.merge(df_fs[["mr_id_x", "hippocampus"]], on="mr_id_x", how="left")

# %% normalize by icv
df["hippocampus_icv"] = df["hippocampus"] / df["icv"] * 100

# %% save
df.to_excel(Path(__file__).parent.parent / "data" / "megamastersheet.xlsx", index=False)

# %%
