# %% ### import modules
from pathlib import Path
import pandas as pd

# %% read results table
results_dir = Path(__file__).resolve().parents[1] / "results"

df = pd.read_excel(results_dir / "cs_b_results.xlsx")

# combine Estimate and Std. Error to Estimate (Std. Error)
# round to 5 decimal places
df["Estimate"] = df["Estimate"].round(3)
df["Std. Error"] = df["Std. Error"].round(3)

df["Beta (SE)"] = df["Estimate"].astype(str) + " (" + df["Std. Error"].astype(str) + ")"

df["PCC (95%CI)"] = (
    df["parcor"].round(2).astype(str)
    + " ("
    + df["parcor_lower"].round(2).astype(str)
    + ","
    + df["parcor_upper"].round(2).astype(str)
    + ")"
)
df["PCC p-value"] = df["parcor_p.value"].apply(lambda x: "<0.001" if x < 0.001 else round(x, 3))

# if the p-value is below 0.01 write < 0.01, otherwise write the value
df["p-value"] = df["Pr(>|t|)"].apply(lambda x: "<0.001" if x < 0.001 else round(x, 3))

# in model remove pad_
df["Model"] = df["pad"].str.replace("pad_", "")

# in model replace gm_icv with GM/ICV
df["Model"] = df["Model"].str.replace("gm_icv", "GM/ICV")

# rename group to Group
df.rename(columns={"group": "Group"}, inplace=True)

# only take type == model+edu
df_edu_suppl = df[df["type"] == "model+edu"]
df = df[df["type"] == "model"]


def pivot_table(df, index, columns, values):
    df = df.pivot(index=index, columns=columns, values=values)
    df.columns = df.columns.swaplevel(0, 1)
    df = df.sort_index(axis=1, level=0)
    return df


def sort_rows(df, index):
    df = df.reindex(index)
    return df


def sort_columns(df, columns):
    df = df.reindex(columns, axis=1, level=0)
    return df


def arrange_table(df):
    df_pivot_edu_suppl = pivot_table(
        df,
        "Model",
        "Group",
        [
            "Beta (SE)",
            "p-value",
            "PCC (95%CI)",
            "PCC p-value",
        ],
    )
    df_pivot_edu_suppl = sort_rows(
        df_pivot_edu_suppl, ["brainageR", "DeepBrainNet", "brainage", "enigma", "pyment", "mccqrnn", "GM/ICV"]
    )
    df_pivot_edu_suppl = sort_columns(df_pivot_edu_suppl, ["CN", "MCI", "AD"])
    return df_pivot_edu_suppl


df_arranged_edu_suppl = arrange_table(df_edu_suppl)
df_arranged = arrange_table(df)

# %% drop Beta (SE)
df_arranged_edu_suppl_beta = df_arranged_edu_suppl.drop(
    columns=[
        "PCC (95%CI)",
        "PCC p-value",
    ],
    level=1,
)
df_arranged_edu_suppl_pcc = df_arranged_edu_suppl.drop(columns=["Beta (SE)", "p-value"], level=1)

# %% drop Beta (SE) for main table
df_arranged_beta = df_arranged.drop(
    columns=[
        "PCC (95%CI)",
        "PCC p-value",
    ],
    level=1,
)

df_arranged_pcc = df_arranged.drop(columns=["Beta (SE)", "p-value"], level=1)
# %% save table
df_arranged_edu_suppl_beta.to_excel(results_dir / "cs_b_results_formatted_edu_suppl_beta.xlsx", index=True)
df_arranged_edu_suppl_pcc.to_excel(results_dir / "cs_b_results_formatted_edu_suppl_pcc.xlsx", index=True)
df_arranged_beta.to_excel(results_dir / "cs_b_results_formatted_beta.xlsx", index=True)
df_arranged_pcc.to_excel(results_dir / "cs_b_results_formatted_pcc.xlsx", index=True)

# %%
