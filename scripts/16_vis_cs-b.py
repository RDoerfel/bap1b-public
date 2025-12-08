# %% ### import modules
from pathlib import Path
import seaborn as sns
import pandas as pd
import numpy as np

import bap1b.figures as figures

# %% ### set paths
data_dir = Path(__file__).resolve().parents[1] / "data"
results_dir = Path(__file__).resolve().parents[1] / "results"

# %% ### load data
df = pd.read_excel(data_dir / "cs_baseline.xlsx")
df_predictions = pd.read_excel(results_dir / "cs_b_predictions.xlsx")

# replace CN in diagnosis with NC
df["diagnosis"] = df["diagnosis"].replace({"CN": "NC"})
df_predictions["diagnosis"] = df_predictions["diagnosis"].replace({"CN": "NC"})
# %%
models = ["brainageR", "DeepBrainNet", "brainage", "enigma", "pyment", "mccqrnn", "gm_icv"]

# Farben festlegen
ki_color_names = ["Plum", "Orange", "Blue", "Black", "Green", "Yellow"]
ki_colors_dark = ["#4F0433", "#B84145", "#002C34", "#000000", "#094334", "#F59A00"]
ki_colors_normal = ["#870052", "#FF876F", "#4DB5BC", "#666666", "#54B986", "#FFC66D"]
ki_colors_light = ["#EDDBE4", "#FFDDD6", "#CCEBED", "#F1F1F1", "#C7ECDC", "#FFE7C2"]

# convert to color palettes
ki_palette_dark = sns.color_palette(ki_colors_dark)
ki_palette_normal = sns.color_palette(ki_colors_normal)
ki_palette_light = sns.color_palette(ki_colors_light)

color_dict_lines = {"NC": ki_palette_normal[4], "MCI": ki_palette_normal[5], "AD": ki_palette_normal[1]}
color_dict_shade = {"NC": ki_palette_light[4], "MCI": ki_palette_light[5], "AD": ki_palette_light[1]}

scatter_color = ki_palette_dark[3]  # Black

# Abbildung und Achsen erstellen
figures.set_rc_params(fontfamily="Arial", small=5, medium=6, big=7)

width = 18.3 - 1
height = 12

### Supplement Figure
fig, axes = figures.get_figures(rows=3, cols=7, unit="cm", figwidth=width, figheight=height, sharex=False, sharey=True)
for d, diag in enumerate(["NC", "MCI", "AD"]):
    for m, model in enumerate(models):
        ax = axes[d, m]
        ylim = (-3, 3)
        yticks = [-2, -1, 0, 1, 2]
        y_col = "ADNI_MEM"
        if model == "gm_icv":
            x_col = model
            x_lable = "GMV/ICV [%]"
            xlim = (20, 60)
            xtiks = [30, 40, 50]
            model = "Gray Matter"
        else:
            x_col = "pad_" + model
            x_lable = "PAD [years]"
            xlim = (-35, 35)
            xtiks = [-20, 0, 20]

        df_model = df[["subject_id", "diagnosis", y_col, x_col]]

        df_diagnosis = df_model[df_model["diagnosis"] == diag]
        df_pred_diag = df_predictions[df_predictions["diagnosis"] == diag]

        # regression line with ci
        x_regression = np.linspace(min(df_diagnosis[x_col]), max(df_diagnosis[x_col]), 100)
        y_regression = df_pred_diag[x_col + "_fit"]
        ci_lwr = df_pred_diag[x_col + "_lwr"]
        ci_upr = df_pred_diag[x_col + "_upr"]

        # plot regression line (predictions) + confidence intervals
        ax.plot(x_regression, y_regression, color=color_dict_lines[diag], lw=1)
        ax.fill_between(x_regression, ci_lwr, ci_upr, color=color_dict_shade[diag], alpha=0.9)
        # scatter plot
        sns.scatterplot(x=x_col, y=y_col, data=df_diagnosis, ax=ax, color=scatter_color, legend=False, s=15, alpha=0.3)
        if (m == 0) & (d == 0):
            figures.set_axes_description(
                ax,
                xlim=xlim,
                ylim=ylim,
                xticks=xtiks,
                xlable=x_lable,
                yticks=yticks,
                ylable="NC\nADNI-Mem",
                title=model,
            )
        elif (m == 0) & (d == 1):
            figures.set_axes_description(
                ax,
                xlim=xlim,
                ylim=ylim,
                xticks=xtiks,
                xlable=x_lable,
                yticks=yticks,
                ylable="MCI\nADNI-Mem",
                title=model,
            )
        elif (m == 0) & (d == 2):
            figures.set_axes_description(
                ax,
                xlim=xlim,
                ylim=ylim,
                xticks=xtiks,
                xlable=x_lable,
                yticks=yticks,
                ylable="AD\nADNI-Mem",
                title=model,
            )
        else:
            figures.set_axes_description(
                ax, xlim=xlim, ylim=ylim, xticks=xtiks, xlable=x_lable, yticks=yticks, title=model
            )
        if d != 2:
            # remove xlabel
            ax.set_xlabel(None)
        if d != 0:
            # remove title
            ax.set_title(None)

fig.set_style(spinewidth=0.9)
fig.save("figure_cs-b_suppl.png", path=results_dir)
fig.save("figure_cs-b_suppl.pdf", path=results_dir)
fig.save("figure_cs-b_suppl.tiff", path=results_dir)

# %% Main Article
figures.set_rc_params(fontfamily="Arial", small=5, medium=6, big=7)
width = 18.3 - 1
height = 5

fig, axes = figures.get_figures(rows=1, cols=7, unit="cm", figwidth=width, figheight=height, sharex=False, sharey=True)
for m, model in enumerate(models):
    ax = axes[m]
    ylim = (-2, 3.5)
    yticks = [-1.5, 0, 1.5, 3]
    y_col = "ADNI_MEM"
    if model == "gm_icv":
        x_col = model
        x_lable = "GMV/ICV [%]"
        xlim = (20, 60)
        xtiks = [30, 40, 50]
        model = "Gray Matter"
    else:
        x_col = "pad_" + model
        x_lable = "PAD [years]"
        xlim = (-35, 35)
        xtiks = [-20, 0, 20]

    df_model = df[["subject_id", "diagnosis", y_col, x_col]]

    df_diagnosis = df_model[df_model["diagnosis"] == "NC"]
    df_pred_diag = df_predictions[df_predictions["diagnosis"] == "NC"]

    # regression line with ci
    x_regression = np.linspace(min(df_diagnosis[x_col]), max(df_diagnosis[x_col]), 100)
    y_regression = df_pred_diag[x_col + "_fit"]
    ci_lwr = df_pred_diag[x_col + "_lwr"]
    ci_upr = df_pred_diag[x_col + "_upr"]

    # plot regression line (predictions) + confidence intervals
    ax.plot(x_regression, y_regression, color=color_dict_lines["NC"], lw=1)
    ax.fill_between(x_regression, ci_lwr, ci_upr, color=color_dict_shade["NC"], alpha=0.9)
    # scatter plot
    sns.scatterplot(x=x_col, y=y_col, data=df_diagnosis, ax=ax, color=scatter_color, legend=False, s=10, alpha=0.3)
    if m == 0:
        figures.set_axes_description(
            ax, xlim=xlim, ylim=ylim, xticks=xtiks, xlable=x_lable, yticks=yticks, ylable="ADNI-Mem", title=model
        )
    else:
        figures.set_axes_description(ax, xlim=xlim, ylim=ylim, xticks=xtiks, xlable=x_lable, yticks=yticks, title=model)

# Apply the updated GridSpec to the figure
fig.set_style(spinewidth=0.9)
fig.save("figure_cs-b.pdf", results_dir)
fig.save("figure_cs-b.png", results_dir)
fig.save("figure_cs-b.tiff", results_dir)

# %%
