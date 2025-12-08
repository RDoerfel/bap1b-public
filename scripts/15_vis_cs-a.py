# %% ### import modules
from pathlib import Path
import seaborn as sns
import pandas as pd

import bap1b.figures as figures

# %% ### set paths
data_dir = Path(__file__).resolve().parents[1] / "data"
results_dir = Path(__file__).resolve().parents[1] / "results"

# %% ### load data
df = pd.read_excel(data_dir / "cs_baseline.xlsx")

# replace CN in diagnosis with NC
df["diagnosis"] = df["diagnosis"].replace({"CN": "NC"})

# %%
models = [
    "brainageR",
    "DeepBrainNet",
    "brainage",
    "enigma",
    "pyment",
    "mccqrnn",
    "gm_icv",
]

# %% Supplementary plot

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
figures.set_rc_params(fontfamily="Arial", small=6, medium=7, big=7)

## Extended Figure 1
width = 18.3 - 1
height = 21 - 1
fig, axes = figures.get_figures(rows=7, cols=4, unit="cm", figwidth=width, figheight=height, sharex=False, sharey=False)

for m, model in enumerate(models):
    xlim = (50, 95)
    xticks = [50, 75, 100]
    if model == "gm_icv":
        y_col = model
        y_lable = "GMV/ICV [%]"
        ylim = (0, 70)
    else:
        y_col = "pad_" + model
        y_lable = f"{model}\nPAD [years]"
        ylim = (-30, 30)

    df_subset = df[["subject_id", "diagnosis", "chron_age", y_col]]

    # Scatterplots mit Regressionslinien
    for i, diagnosis in enumerate(["NC", "MCI", "AD"]):
        ax = axes[m, i]
        subset = df_subset[df_subset["diagnosis"] == diagnosis]
        ax = figures.scatter_plot(
            ax=ax,
            data=subset,
            x_col="chron_age",
            y_col=y_col,
            scatter_color=scatter_color,
            line_color=color_dict_lines[diagnosis],
            label=diagnosis,
        )
        if (i == 0) & (m == 0):
            ax = figures.set_axes_description(ax, xlim=xlim, ylim=ylim, ylable=y_lable, xticks=xticks, title=diagnosis)
        elif (i != 0) & (m == 0):
            ax = figures.set_axes_description(ax, xlim=xlim, ylim=ylim, xticks=xticks, title=diagnosis)
        elif (i == 0) & (m != 0):
            ax = figures.set_axes_description(ax, xlim=xlim, ylim=ylim, xticks=xticks, ylable=y_lable)
        elif (i == 1) & (m == 6):
            ax = figures.set_axes_description(
                ax,
                xlim=xlim,
                ylim=ylim,
                xlable="Chronological Age [years]",
                xticks=xticks,
            )
        else:
            ax = figures.set_axes_description(ax, xlim=xlim, ylim=ylim, xticks=xticks)

    # Boxplots und Violinplot
    ax = axes[m, 3]
    ax = figures.violinplot(ax=ax, data=df_subset, x_col="diagnosis", y_col=y_col, color_dict=color_dict_lines)
    if m == 6:
        ax = figures.set_axes_description(ax, ylim=ylim, xlable="Diagnosis")
    else:
        ax = figures.set_axes_description(ax, ylim=ylim)

fig.set_style(spinewidth=0.9)
fig.save("figure_cs-a_suppl.pdf", results_dir)
fig.save("figure_cs-a_suppl.png", results_dir)
fig.save("figure_cs-a_suppl.tiff", results_dir)

# %% Main Figure
figures.set_rc_params(fontfamily="Arial", small=5, medium=6, big=7)
width = 18.3 - 1
height = 5
fig, axes = figures.get_figures(rows=1, cols=7, unit="cm", figwidth=width, figheight=height, sharex=False, sharey=False)
for m, model in enumerate(models):
    if model == "gm_icv":
        y_col = model
        title = "Gray Matter"
        y_lable = "GMV/ICV [%]"
        ylim = (10, 70)
        yticks = [20, 40, 60]
    else:
        y_col = "pad_" + model
        title = model
        y_lable = "PAD [years]"
        ylim = (-30, 30)
        yticks = [-20, 0, 20]

    df_subset = df[["subject_id", "diagnosis", "chron_age", y_col]]

    # Scatterplots mit Regressionslinien

    # Boxplots und Violinplot
    ax = axes[m]
    ax = figures.violinplot(ax=ax, data=df_subset, x_col="diagnosis", y_col=y_col, color_dict=color_dict_lines)
    if m == 0 or m == 6:
        ax = figures.set_axes_description(ax, ylim=ylim, ylable=y_lable, title=title)
    elif m == 3:
        ax = figures.set_axes_description(ax, ylim=ylim, xlable="Diagnosis", title=title)
    else:
        ax = figures.set_axes_description(ax, ylim=ylim, title=title)

    ax.set_yticks(yticks)

# Apply the updated GridSpec to the figure
fig.set_style(spinewidth=0.9)
fig.save("figure_cs-a.pdf", results_dir)
fig.save("figure_cs-a.png", results_dir)
fig.save("figure_cs-a.tiff", results_dir)

# %%
