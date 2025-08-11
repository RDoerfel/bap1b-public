# %%
from pathlib import Path
import pandas as pd
import numpy as np

# %%
df = pd.read_excel(Path(__file__).parent.parent / "data" / "megamastersheet.xlsx")
# %%
df.head()
# %% keep: subject_id, session_id, gender, chron_age, education_level, diagnosis, ADNI_MEM, pad_brainageR, pad_brainage, pad_DeepBrainNet, pad_pyment, pad_enigma, pad_mccqrnn, keep, gm_icv
df = df[
    [
        "subject_id",
        "session_id",
        "gender",
        "time_from_baseline",
        "mri_field",
        "chron_age",
        "education_level",
        "diagnosis",
        "ADNI_MEM",
        "pad_brainageR",
        "pad_brainage",
        "pad_DeepBrainNet",
        "pad_pyment",
        "pad_enigma",
        "pad_mccqrnn",
        "keep",
        "gm",
        "icv",
        "gm_icv",
        "aes",
        "hippocampus_icv",
    ]
]
# %%
for col in [
    "chron_age",
    "ADNI_MEM",
    "pad_brainageR",
    "pad_brainage",
    "pad_DeepBrainNet",
    "pad_pyment",
    "pad_enigma",
    "pad_mccqrnn",
    "gm",
    "icv",
    "gm_icv",
    "aes",
    "hippocampus_icv",
]:
    mean = df[col].mean()
    std = df[col].std()
    df[col] = np.random.normal(loc=mean, scale=std, size=len(df))
# %% education level should be a random integer between 1 and 20
df["education_level"] = np.random.randint(1, 21, size=len(df))

# %% diagnosis should be a random choice between "CN", "MCI", "AD"
df["diagnosis"] = np.random.choice(["CN", "MCI", "AD"], size=len(df), replace=True)

# %% gender should be randomly sampled from "M" and "F", but stay the same for each subject
for i, subject in enumerate(df["subject_id"].unique()):
    df.loc[df["subject_id"] == subject, "gender"] = np.random.choice(["M", "F"])
    df.loc[df["subject_id"] == subject, "subject_id"] = i
    df.loc[df["subject_id"] == subject, "mri_field"] = np.random.choice([1.5, 3], p=[0.3, 0.7])

# %% randomly drop 20% of the data
df = df.sample(frac=0.8, random_state=42).reset_index(drop=True)

# %% write to file
df.to_excel(Path(__file__).parent.parent / "data" / "megamastersheet_simulated.xlsx", index=False)

# %% write to csv
