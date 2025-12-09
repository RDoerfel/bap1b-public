# %% ### import modules
from pathlib import Path
import seaborn as sns
import pandas as pd

import bap1b.figures as figures

# %% set path
data_dir = Path(__file__).resolve().parents[1] / "data"
results_dir = Path(__file__).resolve().parents[1] / "results"

# %% Load data
for adjustment in ["adjusted", "unadjusted"]:
    for analysis in ["de", "fg"]:
        # skip de for adjusted
        if analysis == "de" and adjustment == "adjusted":
            continue
        print(f"Processing {analysis} {adjustment}")

        df_results = pd.read_excel(results_dir / f"long_{analysis}_results.xlsx")

        # rename CN to NC
        df_results["diagnosis"] = df_results["diagnosis"].replace({"CN": "NC"})

        # ensure integer-like strategy values become plain strings (handles NaN)
        df_results["strategy"] = df_results["strategy"].apply(lambda x: str(x) if pd.notna(x) else "")

        # subset for strategy 2
        if adjustment == "adjusted":
            df_results = df_results[df_results["strategy"] == "2_adjusted"]
        else:
            df_results = df_results[df_results["strategy"] == "2"]

        # settings for plots
        models = ["brainageR", "DeepBrainNet", "brainage", "enigma", "pyment", "mccqrnn", "gm_icv"]
        model_labels = ["brainageR", "DeepBrainNet", "brainage", "enigma", "pyment", "mccqrnn", "Gray Matter"]

        # colors
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

        # plot
        max_width = 8.9
        width = max_width - 1
        height = 5
        fig, axes = figures.get_figures(
            rows=1, cols=2, unit="cm", figwidth=width, figheight=height, sharex=True, sharey=True
        )
        index = len(models) * 4
        revert_index = range(index)[::-1]
        for d, diag in enumerate(["NC", "MCI", "AD"]):
            targets = ["ADNI_MEM", "gm_icv"]
            for i, target in enumerate(targets):
                df_target = df_results[df_results["target"] == target]
                df = df_target[df_target["diagnosis"] == diag]
                ax = axes[i]
                indices = revert_index[d::4]
                if target == "gm_icv":
                    indices = indices[:-1]
                    df = df[df["model"] != "gm_icv"]
                ax.errorbar(
                    df["estimate"],
                    indices,
                    xerr=[df["estimate"] - df["lower"], df["upper"] - df["estimate"]],
                    fmt="o",
                    color=color_dict_lines[diag],
                    ms=2,
                    capsize=2,
                    lw=1,
                    label=diag,
                )
                ax.axvline(0, color="black", lw=0.7, ls="--")
                ax.set_yticks(revert_index[1::4])
                ax.set_xlim(-0.65, 0.65)
                ax.set_xticks([-0.6, -0.3, 0, 0.3, 0.6])

        # set y title and labels
        axes[0].set_ylabel("Model")
        axes[1].set_ylabel("")
        axes[0].set_yticklabels(model_labels)

        # set x title
        axes[0].set_xlabel("Pearson's r (ADNI-Mem)")
        axes[1].set_xlabel("Pearson's r (GMV/ICV)")

        # add legend to first plot, bottom left outside
        axes[0].legend(
            bbox_to_anchor=(-0.3, -0.1),
            ncol=1,
            title="Diagnosis",
            title_fontsize=6,
            fontsize=4,
            frameon=False,
        )

        # set style
        fig.set_style(spinewidth=0.8)

        # save
        fig.save(f"figure_long-{analysis}_{adjustment}.png", results_dir)
        fig.save(f"figure_long-{analysis}_{adjustment}.pdf", results_dir)
        fig.save(f"figure_long-{analysis}_{adjustment}.tiff", results_dir)

# %% Supplement figure (Hippocampal Volume)
target = "hippocampus_icv"
max_width = 8.9
width = max_width - 1
height = 5

for adjustment in ["adjusted", "unadjusted"]:
    for analysis in ["de", "fg"]:
        # skip de for adjusted
        if analysis == "de" and adjustment == "adjusted":
            continue
        print(f"Processing {analysis} {adjustment}")
        fig, ax = figures.get_figures(
            rows=1, cols=1, unit="cm", figwidth=width / 1.5, figheight=height, sharex=True, sharey=True
        )

        df_results = pd.read_excel(results_dir / f"long_{analysis}_results.xlsx")
        df_results["diagnosis"] = df_results["diagnosis"].replace({"CN": "NC"})
        # ensure integer-like strategy values become plain strings (handles NaN)
        df_results["strategy"] = df_results["strategy"].apply(lambda x: str(x) if pd.notna(x) else "")
        # subset for strategy 2
        if adjustment == "adjusted":
            df_results = df_results[df_results["strategy"] == "2_adjusted"]
        else:
            df_results = df_results[df_results["strategy"] == "2"]

        # settings for plots
        models = ["brainageR", "DeepBrainNet", "brainage", "enigma", "pyment", "mccqrnn", "gm_icv"]
        model_labels = ["brainageR", "DeepBrainNet", "brainage", "enigma", "pyment", "mccqrnn", "Gray Matter"]

        # colors
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

        index = len(models) * 4
        revert_index = range(index)[::-1]
        for d, diag in enumerate(["NC", "MCI", "AD"]):
            df_target = df_results[df_results["target"] == target]
            df = df_target[df_target["diagnosis"] == diag]
            ax.errorbar(
                df["estimate"],
                revert_index[d::4],
                xerr=[df["estimate"] - df["lower"], df["upper"] - df["estimate"]],
                fmt="o",
                color=color_dict_lines[diag],
                ms=2,
                capsize=2,
                lw=1,
                label=diag,
            )
            ax.axvline(0, color="black", lw=0.7, ls="--")
            ax.set_yticks(revert_index[1::4])
            ax.set_xlim(-0.65, 0.65)
            ax.set_xticks([-0.6, -0.3, 0, 0.3, 0.6])

        # set y title and labels
        ax.set_ylabel("Model")
        ax.set_yticklabels(model_labels)

        # set x title
        ax.set_xlabel("Pearson's r (HGMV/ICV)")

        # add legend to first plot, bottom left outside
        ax.legend(
            bbox_to_anchor=(-0.3, -0.1),
            ncol=1,
            title="Diagnosis",
            title_fontsize=6,
            fontsize=4,
            frameon=False,
        )

        # set style
        fig.set_style(spinewidth=0.8)

        # save
        fig.save(f"figure_long-{analysis}_{adjustment}_hippocampus_suppl.png", results_dir)
        fig.save(f"figure_long-{analysis}_{adjustment}_hippocampus_suppl.pdf", results_dir)
        fig.save(f"figure_long-{analysis}_{adjustment}_hippocampus_suppl.tiff", results_dir)

# %%
