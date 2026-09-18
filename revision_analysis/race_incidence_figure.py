# -*- coding: utf-8 -*-
"""Figure: annual age-standardized D1 incidence by race/ethnicity (small multiples).

Study: appendiceal goblet cell adenocarcinoma (GCA) case-definition research
in SEER-17, 2000-2023 (revision phase).

Reads the CSV outputs produced by race_incidence_analysis.py:
  outputs/revision_analysis/b3_d1_asr_annual_by_race.csv
  outputs/revision_analysis/b1_d1_asr_by_race_2000_2023.csv

Writes the figure to outputs/revision_analysis/b3_d1_asr_annual_by_race.png
(excluded from version control).

Dependencies: Python 3.x, pandas, matplotlib.
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

OUT = Path("outputs") / "revision_analysis"
ann = pd.read_csv(OUT / "b3_d1_asr_annual_by_race.csv")
b1 = pd.read_csv(OUT / "b1_d1_asr_by_race_2000_2023.csv")

RACE_ORDER = ["NHW", "NHB", "NHAIAN", "NHAPI", "Hispanic"]
LABEL = {"NHW": "Non-Hispanic White", "NHB": "Non-Hispanic Black",
         "NHAIAN": "Non-Hispanic American Indian/\nAlaska Native",
         "NHAPI": "Non-Hispanic Asian or\nPacific Islander",
         "Hispanic": "Hispanic (all races)"}
COLORS = {"NHW": "#1f5fa8", "NHB": "#c0392b", "NHAIAN": "#7f8c8d",
          "NHAPI": "#16866a", "Hispanic": "#d98e04"}

fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.2))
axes = axes.ravel()
for ax in axes[1:5]:
    ax.sharex(axes[0])
for i, race in enumerate(RACE_ORDER):
    ax = axes[i]
    d = ann[ann["race"] == race].sort_values("year")
    ls = "--" if race == "NHAIAN" else "-"
    ax.fill_between(d["year"], d["lcl95"], d["ucl95"], color=COLORS[race], alpha=0.15, linewidth=0)
    ax.plot(d["year"], d["asr_per_100k"], ls, color=COLORS[race], lw=1.8,
            marker="o", ms=3.5, markerfacecolor="white", markeredgewidth=1.1)
    pooled = b1.loc[b1["race"] == race, "asr_per_100k"].iloc[0]
    n = int(b1.loc[b1["race"] == race, "cases"].iloc[0])
    ax.axhline(pooled, color=COLORS[race], lw=0.9, ls=":", alpha=0.8)
    ax.set_title(LABEL[race], fontsize=11, loc="left", fontweight="bold")
    ax.text(0.98, 0.95, f"n = {n}", transform=ax.transAxes, ha="right", va="top", fontsize=9.5, color="#333")
    if race == "NHAIAN":
        ax.text(0.98, 0.82, "unstable (<16 cases):\ndescriptive only", transform=ax.transAxes,
                ha="right", va="top", fontsize=8.5, color="#a04000")
    ax.set_xlim(1999.5, 2023.5)
    ax.set_xticks([2000, 2005, 2010, 2015, 2020])
    ax.tick_params(labelsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", alpha=0.25, lw=0.6)
    ax.set_ylim(bottom=0)

# 6th panel: mini-summary of pooled ASR with CI
ax = axes[5]
d = b1.set_index("race").loc[RACE_ORDER]
ypos = np.arange(len(RACE_ORDER))[::-1]
ax.errorbar(d["asr_per_100k"], ypos,
            xerr=[d["asr_per_100k"] - d["lcl95"], d["ucl95"] - d["asr_per_100k"]],
            fmt="s", color="#333", ecolor="#777", elinewidth=1.2, capsize=3, ms=5)
ax.set_yticks(ypos)
ax.set_yticklabels(["NHW", "NHB", "NHAIAN", "NHAPI", "Hispanic"], fontsize=9.5)
ax.set_title("Pooled 2000–2023 ASR (95% CI)", fontsize=11, loc="left", fontweight="bold")
ax.set_xlabel("ASR per 100,000 person-years", fontsize=9.5)
ax.tick_params(labelsize=9)
ax.spines[["top", "right"]].set_visible(False)
ax.grid(axis="x", alpha=0.25, lw=0.6)
ax.set_xlim(left=0)

for ax in axes[3:5]:
    ax.set_xlabel("Year of diagnosis", fontsize=10)
for ax in [axes[0], axes[3]]:
    ax.set_ylabel("ASR per 100,000 person-years", fontsize=10)

fig.suptitle("Age-standardized incidence of appendiceal goblet cell adenocarcinoma (ICD-O-3 8243/3)\n"
             "by race/ethnicity, SEER-17, 2000–2023 (2000 US Standard Population; shaded band = 95% CI)",
             fontsize=12.5, fontweight="bold", y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.94])
png = OUT / "b3_d1_asr_annual_by_race.png"
fig.savefig(png, dpi=200, bbox_inches="tight")
print("saved", png)
