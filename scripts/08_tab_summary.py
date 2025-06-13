# %% ### import modules
from pathlib import Path
import pandas as pd

# %% read results table
results_dir = Path(__file__).resolve().parents[1] / "results"
df = pd.read_excel(results_dir / "summary_table.xlsx")

# %% rename and structure columns:
new_columns = [
    ("", "Group"),  # diagnosis -> Group
    ("Scans", "M000"),  # N_M000 -> M000 under Scans
    ("Scans", "M012"),  # N_M012 -> M012 under Scans
    ("Scans", "M024"),  # N_M024 -> M024 under Scans
    ("Scans", "M036"),  # N_M036 -> M036 under Scans
    ("Scans", "M048"),  # N_M048 -> M048 under Scans
    ("Baseline", "Age [years] mean (SD)"),  # Age_mean_sd -> Age mean (SD)
    ("Baseline", "Sex F/M"),  # Sex_FM -> Sex F/M
    ("Baseline", "Education [years] mean (SD)"),  # EDU_mean_sd -> Education [years] mean (SD)
    ("Baseline", "ADNI MEM mean (SD)"),  # ADNI_MEM_mean_sd -> ADNI MEM mean (SD)
    ("Baseline", "AES mean (SD)"),  # AES_mean_sd -> AES mean (SD)
]

# Assign the new MultiIndex to the DataFrame columns
df.columns = pd.MultiIndex.from_tuples(new_columns)

print(df)
# %% save
df.to_excel(results_dir / "summary_table_formatted.xlsx")

# %%
