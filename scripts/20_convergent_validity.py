# %%
from pathlib import Path
import pandas as pd
from scipy.stats import pearsonr
import statsmodels.api as sm
import numpy as _np

# %% load data
data_dir = Path(__file__).resolve().parents[1] / "data"
results_dir = Path(__file__).resolve().parents[1] / "results"
df = pd.read_excel(data_dir / "cs_baseline.xlsx")

# %% pairplot of models and gm_icv
# filter on CN
df["gm_icv"] = df["gm"] / df["icv"]
models = ["brainageR", "DeepBrainNet", "brainage", "enigma", "pyment", "mccqrnn"]
pad_models = ["pad_" + model for model in models]
pad_models_gm_icv = pad_models + ["gm_icv"]
keep = pad_models_gm_icv + ["diagnosis"]
data = df[keep]

# %% compute the correlation and ci between each package using scipy
df_results = pd.DataFrame(index=pad_models_gm_icv, columns=pad_models_gm_icv)

# %%
df_cn = data[df["diagnosis"] == "CN"].copy()


# %% compute pearson's r between "enigma", "pyment"
col1 = "pad_enigma"
col2 = "pad_brainage"


def compute_pearsonr(data, col1, col2):
    # center
    data[col1] = data[col1] - data[col1].mean()
    data[col2] = data[col2] - data[col2].mean()
    result = pearsonr(data[col1], data[col2])
    corr = round(result[0], 2)
    pval = result[1]
    return corr, pval


corr_enigma_brainage, pval_enigma_brainage = compute_pearsonr(df_cn, col1, col2)


#  fit linear regression line between "enigma", "pyment"
def fit_linear_regression(data, col1, col2, col1_new=10):
    X = data[col1]
    y = data[col2]
    mu_x = X.mean()
    mu_y = y.mean()
    sd_x = X.std()
    sd_y = y.std()

    # center
    X = X - mu_x
    y = y - mu_y

    # unit variance
    X = X / sd_x
    y = y / sd_y

    model = sm.OLS(y, X).fit()
    slope = model.params[col1]

    # prediction interval
    X_new = _np.array([col1_new])
    X_new = X_new - mu_x
    X_new = X_new / sd_x
    pred = model.get_prediction(X_new)
    pred_summary = pred.summary_frame(alpha=0.05)
    pred_summary["mean"][0] = pred_summary["mean"][0] * sd_y + mu_y
    pred_summary["obs_ci_lower"][0] = pred_summary["obs_ci_lower"][0] * sd_y + mu_y
    pred_summary["obs_ci_upper"][0] = pred_summary["obs_ci_upper"][0] * sd_y + mu_y

    return slope, pred_summary


def get_pred_interval(data, col1, col2, col1_new=10):
    slope, pred_summary = fit_linear_regression(data, col1, col2, col1_new)
    print(f"Pearson's r between {col1} and {col2}: {slope} ")
    print(f"95% prediction interval at {col1} x={col1_new} years:")
    print(
        f"prediction {pred_summary['mean'][0]:.2f} ({pred_summary['obs_ci_lower'][0]:.2f}, {pred_summary['obs_ci_upper'][0]:.2f})"
    )


get_pred_interval(df_cn, col1="pad_enigma", col2="pad_brainage", col1_new=5)
get_pred_interval(df_cn, col1="pad_enigma", col2="pad_pyment", col1_new=5)


# %% manual
