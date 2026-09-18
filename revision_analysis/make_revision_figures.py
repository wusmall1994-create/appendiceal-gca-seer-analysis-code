# -*- coding: utf-8 -*-
"""Revision figures: redraw Figure 1 (aligned x-axes, 2019/2021 reference lines)
and create Figure S4 (subgroup forest plot).

Study: appendiceal goblet cell adenocarcinoma (GCA) case-definition research
in SEER-17, 2000-2023 (revision phase).

Required inputs (all produced by earlier scripts in this repository):
  outputs/formal_next_phase/formal_incidence_annual_fay_feuer.csv
      (R/01_main_analysis.R)
  outputs/formal_next_phase/coding_composition_by_period_formal.csv
      (R/01_main_analysis.R)
  outputs/formal_next_phase/formal_incidence_period_definition_effect.csv
      (R/01_main_analysis.R)
  outputs/revision_analysis/b2_d1_asr_by_race_periods_rr.csv
      (revision_analysis/race_incidence_analysis.py)
  outputs/revision_analysis/task3_rate_ratios_2020_2023_vs_2000_2004.csv
      (revision_analysis/run_phase1.py)

Survival hazard ratios and equivalence assessments shown in Figure S4 are the
values reported in the revision (from revision_analysis/run_phase1.py and
R/01_main_analysis.R); they are entered as plotting constants here.

Outputs are written to outputs/revision_analysis/ (excluded from version
control).

Dependencies: Python 3.x, numpy, matplotlib.
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

FORMAL = Path("outputs") / "formal_next_phase"
REV = Path("outputs") / "revision_analysis"
REV.mkdir(parents=True, exist_ok=True)

C_D1 = "#1F5C8B"   # dark blue
C_D2 = "#DD8452"   # orange
C_D3 = "#8172B3"   # purple

DEF_LABEL = {
    "D1_strict_GCA": "Strict 8243",
    "D2_legacy_GCC_adenocarcinoid": "8243 + 8245",
    "D3_historical_expanded_spectrum": "8243 + 8244 + 8245",
}
DEF_COLOR = {"D1_strict_GCA": C_D1, "D2_legacy_GCC_adenocarcinoid": C_D2, "D3_historical_expanded_spectrum": C_D3}

PERIODS = ["2000-2004", "2005-2009", "2010-2014", "2015-2019", "2020-2023"]
MID = {"2000-2004": 2002, "2005-2009": 2007, "2010-2014": 2012, "2015-2019": 2017, "2020-2023": 2021.5}
WID = {"2000-2004": 4.6, "2005-2009": 4.6, "2010-2014": 4.6, "2015-2019": 4.6, "2020-2023": 3.6}

XLIM = (1999.3, 2023.8)
XTICKS = [2000, 2005, 2010, 2015, 2020]


def add_ref_lines(ax, y_base):
    """Draw 2019/2021 reference lines with small rotated tags near the axis base."""
    ax.axvline(2019, color="grey", linestyle="--", linewidth=1.1, zorder=1)
    ax.axvline(2021, color="grey", linestyle=":", linewidth=1.4, zorder=1)
    ax.text(2019, y_base, " 2019", color="grey", fontsize=8, ha="center", va="bottom", rotation=90)
    ax.text(2021, y_base, " 2021", color="grey", fontsize=8, ha="center", va="bottom", rotation=90)


# ---------------- Figure 1 ----------------
annual = {}
with (FORMAL / "formal_incidence_annual_fay_feuer.csv").open(encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        annual.setdefault(row["definition"], []).append(
            (int(row["year"]), float(row["age_adjusted_rate_per_100k"]),
             float(row["lcl_95"]), float(row["ucl_95"])))

comp = {p: {} for p in PERIODS}
with (FORMAL / "coding_composition_by_period_formal.csv").open(encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        comp[row["period"]][row["histology"]] = float(row["share_percent"])

effect = {}
with (FORMAL / "formal_incidence_period_definition_effect.csv").open(encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        if row["period"] == "2000-2023" or row["definition"] == "D1_strict_GCA":
            continue
        effect.setdefault(row["definition"], []).append(
            (MID[row["period"]], float(row["rate_inflation_vs_D1_percent"]),
             float(row["inflation_lcl_95"]), float(row["inflation_ucl_95"])))

fig = plt.figure(figsize=(11.5, 9.2), dpi=300)
gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], hspace=0.30, wspace=0.24,
                      left=0.075, right=0.975, top=0.93, bottom=0.09)
axa = fig.add_subplot(gs[0, :])
axb = fig.add_subplot(gs[1, 0])
axc = fig.add_subplot(gs[1, 1])

# Panel a: annual ASR
for d in ("D1_strict_GCA", "D2_legacy_GCC_adenocarcinoid", "D3_historical_expanded_spectrum"):
    ys = np.array([r[0] for r in annual[d]])
    rate = np.array([r[1] for r in annual[d]])
    lo = np.array([r[2] for r in annual[d]])
    hi = np.array([r[3] for r in annual[d]])
    axa.fill_between(ys, lo, hi, color=DEF_COLOR[d], alpha=0.16, linewidth=0)
    axa.plot(ys, rate, color=DEF_COLOR[d], linewidth=2.2, label=DEF_LABEL[d])
axa.set_title("a  Reported incidence depends on case definition", loc="left", fontsize=12, fontweight="bold")
axa.set_ylabel("Age-adjusted rate per 100,000", fontsize=10.5)
axa.set_xlim(*XLIM)
axa.set_xticks(XTICKS)
axa.set_ylim(0, 0.33)
add_ref_lines(axa, 0.004)
axa.legend(fontsize=9, loc="upper left", frameon=False, ncol=3)
axa.text(2000.2, 0.245, "Dashed line: D1 joinpoint (2019)\nDotted line: ICD-O-3.2 preferred terms (2021)",
         fontsize=8.5, color="grey", va="top")
axa.tick_params(labelsize=9.5)
for s in ("top", "right"):
    axa.spines[s].set_visible(False)

# Panel b: coding composition (stacked bars on continuous axis)
bottoms = np.zeros(len(PERIODS))
for hist, color, label in (("8243", C_D1, "8243 strict GCA"),
                           ("8244", C_D3, "8244 mixed"),
                           ("8245", C_D2, "8245 adenocarcinoid")):
    vals = np.array([comp[p].get(hist, 0.0) for p in PERIODS])
    axb.bar([MID[p] for p in PERIODS], vals, width=[WID[p] for p in PERIODS],
            bottom=bottoms, color=color, edgecolor="white", linewidth=0.6, label=label)
    bottoms += vals
axb.set_title("b  Coding composition of the expanded spectrum", loc="left", fontsize=12, fontweight="bold")
axb.set_ylabel("Share of expanded spectrum (%)", fontsize=10.5)
axb.set_xlabel("Year of diagnosis", fontsize=10.5)
axb.set_xlim(*XLIM)
axb.set_xticks(XTICKS)
axb.set_ylim(0, 100)
add_ref_lines(axb, 1.5)
axb.legend(fontsize=8.5, loc="center left", frameon=True, facecolor="white", framealpha=0.95,
           edgecolor="none")
axb.tick_params(labelsize=9.5)
for s in ("top", "right"):
    axb.spines[s].set_visible(False)

# Panel c: definition effect over time
for d, label in (("D2_legacy_GCC_adenocarcinoid", "8243 + 8245"),
                 ("D3_historical_expanded_spectrum", "Expanded spectrum")):
    pts = effect[d]
    x = np.array([p[0] for p in pts])
    y = np.array([p[1] for p in pts])
    lo = np.array([p[2] for p in pts])
    hi = np.array([p[3] for p in pts])
    axc.errorbar(x, y, yerr=[y - lo, hi - y], color=DEF_COLOR[d], marker="o", markersize=5,
                 linewidth=1.8, capsize=3, label=label)
axc.axhline(0, color="black", linewidth=0.8)
axc.set_title("c  Definition effect over time", loc="left", fontsize=12, fontweight="bold")
axc.set_ylabel("Rate change vs strict 8243 (%)", fontsize=10.5)
axc.set_xlabel("Year of diagnosis", fontsize=10.5)
axc.set_xlim(*XLIM)
axc.set_xticks(XTICKS)
axc.set_ylim(-8, 145)
add_ref_lines(axc, -6)
axc.legend(fontsize=8.5, loc="lower left", frameon=False)
axc.tick_params(labelsize=9.5)
for s in ("top", "right"):
    axc.spines[s].set_visible(False)

out1 = REV / "Figure1_case_definition_and_coding_migration.png"
fig.savefig(out1, dpi=300, facecolor="white")
plt.close(fig)
print("Figure1:", out1)

# ---------------- Figure S4: forest plot ----------------
rows = []  # (section_title_or_None, label, est, lo, hi, note)
rows.append(("Age-standardized rate ratio by race/ethnicity (2010–2023 vs 2000–2009)", None, None, None, None, None))
with (REV / "b2_d1_asr_by_race_periods_rr.csv").open(encoding="utf-8-sig") as f:
    race_lbl = {"NHW": "Non-Hispanic White", "NHB": "Non-Hispanic Black",
                "NHAPI": "Non-Hispanic Asian/Pacific Islander", "Hispanic": "Hispanic (all races)"}
    for row in csv.DictReader(f):
        if row["race"] in race_lbl:
            rows.append((None, race_lbl[row["race"]], float(row["rr_p2_vs_p1"]),
                         float(row["rr_lcl95"]), float(row["rr_ucl95"]), None))
rows.append(("Age-specific rate ratio (2020–2023 vs 2000–2004)", None, None, None, None, None))
with (REV / "task3_rate_ratios_2020_2023_vs_2000_2004.csv").open(encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        rows.append((None, f"Ages {row['age_group']} years", float(row["rate_ratio"]),
                     float(row["rr_lcl"]), float(row["rr_ucl"]), None))
rows.append(("Five-year overall survival, unadjusted hazard ratio", None, None, None, None, None))
rows.append((None, "Black vs White", 1.26, 0.88, 1.80, None))
rows.append((None, "Hispanic vs non-Hispanic", 1.05, 0.67, 1.65, None))
rows.append(("Adjusted OS hazard ratio by histology code (equivalence margin 0.80–1.25)", None, None, None, None, None))
rows.append((None, "8244/3 vs 8243/3", 1.063, 0.909, 1.243, "equivalent"))
rows.append((None, "8245/3 vs 8243/3", 1.134, 0.925, 1.392, "indeterminate"))

fig, ax = plt.subplots(figsize=(9.2, 8.2), dpi=600)
yticks, ylabels = [], []
y = 0
for section, label, est, lo, hi, note in rows:
    if section:
        y -= 1.0
        yticks.append(y)
        ylabels.append(("section", section))
        continue
    y -= 0.62
    color = "#1F5C8B" if note != "indeterminate" else "#888888"
    ax.plot([lo, hi], [y, y], color=color, linewidth=1.4, solid_capstyle="round")
    ax.plot([est], [y], marker="s", markersize=6.5, color=color)
    yticks.append(y)
    suffix = {"equivalent": "  (equivalence established)", "indeterminate": "  (underpowered)"}.get(note, "")
    ylabels.append(("item", f"{label}   {est:.2f} ({lo:.2f}–{hi:.2f}){suffix}"))

ax.axvline(1.0, color="black", linewidth=1.0)
ax.axvline(0.80, color="#C44E52", linestyle="--", linewidth=1.1)
ax.axvline(1.25, color="#C44E52", linestyle="--", linewidth=1.1)
ax.text(2.4, y - 0.45, "Solid line: null (1.0)\nDashed red: equivalence\nmargin (0.80–1.25)", fontsize=8.5,
        color="#555555", va="top", ha="left")

ax.set_yticks(yticks)
ax.set_yticklabels([t[1] for t in ylabels], fontsize=9)
for tick, (kind, _) in zip(ax.get_yticklabels(), ylabels):
    if kind == "section":
        tick.set_fontweight("bold")
        tick.set_fontsize(9.5)
ax.set_xscale("log")
ax.set_xlim(0.55, 22)
ax.set_xticks([0.67, 1, 2, 3, 5, 10, 17])
ax.set_xticklabels(["0.67", "1", "2", "3", "5", "10", "17"], fontsize=9.5)
ax.set_xlabel("Rate ratio or hazard ratio (log scale)", fontsize=10.5)
ax.set_ylim(y - 1.35, 1.0)
ax.tick_params(axis="y", length=0)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
fig.tight_layout()
out4 = REV / "FigureS4_subgroup_forest_plot.png"
fig.savefig(out4, dpi=600, facecolor="white")
plt.close(fig)
print("FigureS4:", out4)
