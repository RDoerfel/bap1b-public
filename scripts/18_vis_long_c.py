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
df_results = pd.read_excel(results_dir / "long_c_results.xlsx")

df_preds = pd.read_excel(results_dir / "long_c_preds.xlsx")

df_data = pd.read_excel(data_dir / "data_long_c.xlsx")

# replce CN with NC
df_data["diagnosis"] = df_data["diagnosis"].replace({"CN": "NC"})

models = ["brainageR", "DeepBrainNet", "brainage", "enigma", "pyment", "mccqrnn", "gm_icv"]

# %% plot probability usin logisticc regression. df_results has a column model, intercept and coef.
# I want to plot the probability of the model for the range of -40 to 20.

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

# set overall parameters
figures.set_rc_params(fontfamily="Arial", small=5, medium=6, big=7)

#
max_width = 18.3
width = max_width - 1
height = 5
print("height", height, "width", width)
fig, axes = figures.get_figures(rows=1, cols=7, unit="cm", figwidth=width, figheight=height, sharex=False, sharey=True)

for i, model in enumerate(models):
    print(model)
    if model == "gm_icv":
        model = "gm_icv"
        title = "Gray Matter"
        xlabel = "GMV/ICV [%]"
    else:
        model = "pad_" + model
        title = model[4:]
        xlabel = "PAD [years]"

    data_col = model + "_base"

    # set x based on data range
    x = np.linspace(df_data[data_col].min(), df_data[data_col].max(), 100)
    ax = axes[i]

    y = df_preds[data_col + "_pred"]
    y_cinf_lower = df_preds[data_col + "_ci_lower"]
    y_cinf_upper = df_preds[data_col + "_ci_upper"]

    # plot expit
    ax.plot(x, y, color=color_dict_lines["NC"])

    # plot confidence intervall
    ax.fill_between(x, y_cinf_lower, y_cinf_upper, color=color_dict_shade["NC"], alpha=0.9)

    # plot data centered data
    sns.rugplot(df_data[data_col], ax=ax, color="black", height=0.05, lw=0.6, alpha=0.1)

    ax.set_ylabel("")
    ax.set_xlabel("")
    if i == 0:
        ax.set_ylabel("Probability of change")
    elif i == 3 or i == 6:
        ax.set_xlabel(xlabel)

    ax.set_yticks([0, 0.1, 0.2, 0.3, 0.4])
    # ax.set_xlim(xlim)
    ax.set_ylim((0, 0.4))
    ax.set_title(title)

fig.set_style(spinewidth=0.9)
fig.save("figure_long-c.png", results_dir)
fig.save("figure_long-c.pdf", results_dir)
fig.save("figure_long-c.tiff", results_dir)


# %% ROC curve
roc_data = pd.read_excel(results_dir / "long_c_roc_data.xlsx")

# data is in vide format. plot TPR vs FPR for each model. hue is the model
fig, ax = figures.get_figures(rows=1, cols=1, unit="cm", figwidth=7, figheight=7, sharex=False, sharey=True)
sns.lineplot(data=roc_data, x="FPR", y="TPR", hue="model", linewidth=0.9, ax=ax)
ax.plot([0, 1], [0, 1], color="black", linestyle="--", linewidth=0.5)
fig.set_style(spinewidth=0.9)
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_xticks([0, 0.5, 1])
ax.set_yticks([0, 0.5, 1])
# add legend to the right
ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))

fig.save("figure_long-c_roc.png", results_dir)
fig.save("figure_long-c_roc.pdf", results_dir)
fig.save("figure_long-c_roc.tiff", results_dir)

# %%
