# %%
from pathlib import Path
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

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


df_results = compute_correlations(df, pad_models_gm_icv)


# %%
def custom_pairplot(data, outcome_columns, figsize=(10, 10), cmap="coolwarm"):
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
    # Compute correlation matrix
    corr_matrix = data[outcome_columns].corr()
    corr_labels = compute_correlations(data, outcome_columns)

    # Determine number of columns
    n = len(outcome_columns)

    # remove pad_ from column names
    column_names = [name.replace("pad_", "") for name in outcome_columns]

    # replace gm_icv with GM/ICV
    column_names = [name.replace("gm_icv", "Gray Matter") for name in column_names]

    # Create figure and axes
    fig, axes = plt.subplots(n, n, figsize=figsize)

    # Iterate through all subplot combinations
    for i in range(n):
        for j in range(n):
            ax = axes[i, j]

            # Diagonal: show column names
            if i == j:
                # histplot using plt.hist to avoid labels
                ax.hist(data[outcome_columns[i]], bins=30, color="tab:blue", alpha=0.9)

                # remove left, top and right spines
                ax.spines[["left", "top", "right"]].set_visible(False)
                ax.set_yticks([0])

            # Lower triangle: correlation heatmap
            elif i > j:
                ax.scatter(data[outcome_columns[j]], data[outcome_columns[i]], alpha=0.9, s=0.05, color="tab:blue")

            # Upper triangle: scatter plot
            elif i < j:
                # get text from corr_labels
                sns.heatmap(
                    [[corr_matrix.iloc[i, j]]],
                    ax=ax,
                    cbar=False,
                    cmap=cmap,
                    vmax=1,
                    vmin=-1,
                )
                ax.axis("off")

                # add text to heatmap
                text = corr_labels.iloc[i, j]
                ax.text(0.5, 0.5, text, ha="center", va="center", fontsize=8, color="black")

            if i == (n - 1):
                ax.set_xlabel(column_names[j], rotation=0)

            if j == 0:
                ax.set_ylabel(column_names[i], rotation=90)

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
fig = custom_pairplot(data, outcome_columns=pad_models_gm_icv, figsize=(figwidth, figwidth))

fig.savefig(results_dir / "correlation_plot.png", dpi=300)
fig.savefig(results_dir / "correlation_plot.pdf", dpi=300)

# %%
