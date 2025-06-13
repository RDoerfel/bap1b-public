# %%
from pathlib import Path
import seaborn as sns
import pandas as pd

import bap1b.figures as figures

# %% ### set paths
data_dir = Path(__file__).resolve().parents[1] / "data"
results_dir = Path(__file__).resolve().parents[1] / "results"

# %% ### load data
df = pd.read_excel(data_dir / "cs_baseline.xlsx")

# %% plot ADNI_MEM vs age, ADNI_mem vs educaion_level

# Farben festlegen
colors = sns.color_palette("Dark2", 3)
colors_hue = sns.color_palette("Set2", 3)
color_dict = {"CN": colors[0], "MCI": colors[1], "AD": colors[2]}
scatter_color = "black"  # sns.color_palette("Set2", 8)[7]

# Abbildung und Achsen erstellen
figures.set_rc_params(fontfamily="Arial", small=6, medium=8, big=8)

# figures
fig, axes = figures.get_figures(rows=1, cols=2, unit="cm", figwidth=7, figheight=4, sharex=False, sharey=True)

# age
ax = axes[0]
sns.scatterplot(
    data=df, x="chron_age", y="ADNI_MEM", hue="diagnosis", ax=ax, legend=False, palette=colors_hue, size=0.5
)
for diagnosis in ["CN", "MCI", "AD"]:
    df_diagnosis = df[df["diagnosis"] == diagnosis]
    sns.regplot(data=df_diagnosis, x="chron_age", y="ADNI_MEM", ax=ax, scatter=False, color=color_dict[diagnosis])

ax.set_xlim([45, 95])
ax.set_xticks([50, 60, 70, 80, 90])
ax.set_yticks([-3, -1.5, 0, 1.5, 3])
# education
ax = axes[1]
sns.scatterplot(
    data=df, x="education_level", y="ADNI_MEM", hue="diagnosis", ax=ax, legend=False, palette=colors_hue, size=0.5
)
for diagnosis in ["CN", "MCI", "AD"]:
    df_diagnosis = df[df["diagnosis"] == diagnosis]
    sns.regplot(data=df_diagnosis, x="education_level", y="ADNI_MEM", ax=ax, scatter=False, color=color_dict[diagnosis])

ax.set_xlim([0, 25])
ax.set_xticks([5, 10, 15, 20])


fig.set_style()

fig.save("cs_Adni_mem_age_edu.jpg", results_dir)

# %% education_level vs pad_pyment
fig, axes = figures.get_figures(rows=1, cols=1, unit="cm", figwidth=4.5, figheight=4.5, sharex=False, sharey=True)
sns.scatterplot(
    data=df, x="education_level", y="pad_pyment", ax=axes, legend=False, color=scatter_color, size=0.5, alpha=0.3
)
sns.regplot(data=df_diagnosis, x="education_level", y="pad_pyment", ax=axes, scatter=False, color=scatter_color)
fig.set_style()


fig
fig.save("cs_pad-pyment_edu.jpg", results_dir)
