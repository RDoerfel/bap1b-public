# %%
from pathlib import Path
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import pingouin as pg
import numpy as _np

import bap1b.figures as figures

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


# lambda for p-value stars if p < 0.001 then ***, if p < 0.01 then **, if p < 0.05 then *
def pvalue_stars(p):
    return " ***" if p < 0.001 else " **" if p < 0.01 else " *" if p < 0.05 else ""


def compute_correlations(data, models):
    df_results = pd.DataFrame(index=models, columns=models)
    for model1 in models:
        for model2 in models:
            result = pearsonr(data[model1], data[model2])
            corr = round(result[0], 2)
            ci_lower, ci_upper = result.confidence_interval()
            ci_lower = round(ci_lower, 2)
            ci_upper = round(ci_upper, 2)
            pstars = pvalue_stars(result[1])
            df_results.loc[model1, model2] = f"{corr}\n({ci_lower}, {ci_upper})"
    return df_results


def compute_partial_correlations(data, models):
    df_results = pd.DataFrame(index=models, columns=models)
    for model1 in models:
        for model2 in models:
            if model1 == model2:
                df_results.loc[model1, model2] = "1.0\n(1.0, 1.0)"
                continue
            # create dummy variables for diagnosis using pandas and keep as covariates list
            data = data.copy()
            diag_dummies = pd.get_dummies(data["diagnosis"], prefix="diagnosis", drop_first=True)
            data_dummy = pd.concat([data, diag_dummies], axis=1)
            covars = list(diag_dummies.columns)
            data_subset = data_dummy[[model1, model2] + covars]
            pcorr = pg.partial_corr(data=data_subset, x=model1, y=model2, covar=covars, method="pearson")
            corr = round(pcorr["r"].values[0], 2)
            ci_lower = round(pcorr["CI95%"][0][0], 2)
            ci_upper = round(pcorr["CI95%"][0][1], 2)
            df_results.loc[model1, model2] = f"{corr}\n({ci_lower}, {ci_upper})"
    return df_results


df_corrs = compute_correlations(df, pad_models_gm_icv)
df_partial_corrs = compute_partial_correlations(df, pad_models_gm_icv)


# %%
def custom_pairplot(data, outcome_columns, figsize=(10, 10), cmap="coolwarm", partial=False):
    """
    Create a custom pairplot with correlation heatmaps in lower triangle
    and scatter plots in upper triangle.

    Parameters:
    -----------
    data : pandas.DataFrame
        Input dataframe with numerical columns
    figsize : tuple, optional (default=(10,10))
        Size of the entire figure
    cmap : str, optional (default='coolwarm')
        Colormap for the correlation heatmap

    Returns:
    --------
    matplotlib.figure.Figure
        The created pairplot figure
    """
    # colors:
    ki_color_names = ["Plum", "Orange", "Blue", "Black", "Green", "Yellow"]
    ki_colors_dark = ["#4F0433", "#B84145", "#002C34", "#000000", "#094334", "#F59A00"]
    ki_colors_normal = ["#870052", "#FF876F", "#4DB5BC", "#666666", "#54B986", "#FFC66D"]
    ki_colors_light = ["#EDDBE4", "#FFDDD6", "#CCEBED", "#F1F1F1", "#C7ECDC", "#FFE7C2"]

    # convert to color palettes
    ki_palette_dark = sns.color_palette(ki_colors_dark)
    ki_palette_normal = sns.color_palette(ki_colors_normal)
    ki_palette_light = sns.color_palette(ki_colors_light)

    color_dict_lines = {"CN": ki_palette_normal[4], "MCI": ki_palette_normal[5], "AD": ki_palette_normal[1]}
    color_dict_shade = {"CN": ki_palette_light[4], "MCI": ki_palette_light[5], "AD": ki_palette_light[1]}

    # Compute correlation matrix
    if partial:
        corr_labels = compute_partial_correlations(data, outcome_columns)
    else:
        corr_labels = compute_correlations(data, outcome_columns)

    # Determine number of columns
    n = len(outcome_columns)

    # remove pad_ from column names
    column_names = [name.replace("pad_", "") for name in outcome_columns]

    # replace gm_icv with GM/ICV
    column_names = [name.replace("gm_icv", "Gray Matter") for name in column_names]

    # Create figure and axes
    fig, axes = plt.subplots(n, n, figsize=figsize)

    # ensure axes is 2D (handles n==1)
    axes = _np.atleast_2d(axes)

    # Iterate through all subplot combinations
    for i in range(n):
        for j in range(n):
            ax = axes[i, j]

            # remove ticks
            ax.set_xticks([])
            ax.set_yticks([])

            # Diagonal: histogram
            if i == j:
                for diag in ["CN", "MCI", "AD"]:
                    data_diag = data[data["diagnosis"] == diag]
                    ax.hist(data_diag[outcome_columns[i]].dropna(), bins=30, color=color_dict_lines[diag], alpha=0.8)
                ax.spines[["left", "top", "right"]].set_visible(False)
                ax.set_yticks([])

            # upper triangle: correlation heatmap (i < j)
            elif i < j:
                text = corr_labels.iloc[i, j]
                corr = float(text.split("\n")[0])
                sns.heatmap(
                    [[corr]],
                    ax=ax,
                    cbar=False,
                    cmap=cmap,
                    vmax=1,
                    vmin=-1,
                )
                ax.axis("off")

                # add text to heatmap
                ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=8, color="black")

                ax.set_yticklabels([])
                ax.set_xticklabels([])

            # Lower triangle: scatter plot (i < j)
            else:
                for diag in ["CN", "MCI", "AD"]:
                    data_diag = data[data["diagnosis"] == diag]
                    ax.scatter(
                        data_diag[outcome_columns[j]],
                        data_diag[outcome_columns[i]],
                        alpha=0.6,
                        s=1,
                        color=color_dict_lines[diag],
                    )
                ax.spines[["top", "right"]].set_visible(False)

            if j == 0:
                ax.set_ylabel(column_names[i], rotation=90, fontsize=8)
            if i == n - 1:
                ax.set_xlabel(column_names[j], rotation=0, fontsize=8)

    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            ax.set_yticklabels([])
            ax.set_xticklabels([])

    plt.tight_layout()
    return fig


# set rc params small for better visibility (seperate for labels and ticks)
figures.set_rc_params(fontfamily="Arial", small=5, medium=6, big=7)

inch_to_cm = 2.54
figwidth = 17 / inch_to_cm
fig = custom_pairplot(data, outcome_columns=pad_models_gm_icv, figsize=(figwidth, figwidth), partial=True)

fig.savefig(results_dir / "correlation_plot.png", dpi=300)
fig.savefig(results_dir / "correlation_plot.pdf", dpi=300)

# %% generate three different group results including CI for data['diagnosis'] == {CN, MCI, AD}
for diag in ["CN", "MCI", "AD"]:
    data_diag = data[data["diagnosis"] == diag]
    df_results_diag = compute_correlations(data_diag, pad_models_gm_icv)
    df_results_diag.to_excel(results_dir / f"correlation_results_{diag}.xlsx")
    fig_diag = custom_pairplot(
        data_diag, outcome_columns=pad_models_gm_icv, figsize=(figwidth, figwidth), partial=False
    )
    fig_diag.savefig(results_dir / f"correlation_plot_{diag}.png", dpi=300)
    fig_diag.savefig(results_dir / f"correlation_plot_{diag}.pdf", dpi=300)

# %%
fig_pcc = custom_pairplot(data, outcome_columns=pad_models_gm_icv, figsize=(figwidth, figwidth), partial=True)
fig_pcc.savefig(results_dir / "correlation_plot_partial.png", dpi=300)
fig_pcc.savefig(results_dir / "correlation_plot_partial.pdf", dpi=300)

# %%
