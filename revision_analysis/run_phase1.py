# -*- coding: utf-8 -*-
"""Revision-phase analyses (Tasks 1-4) for the appendiceal GCA SEER study.

Study: how nested ICD-O-3 case definitions affect reported incidence, stage
distribution, and survival estimates for appendiceal goblet cell
adenocarcinoma (GCA), SEER-17, 2000-2023.

Tasks:
  1. Table 1 characteristics of the strict (D1, ICD-O-3 8243/3) cohort by
     race and by Hispanic ethnicity (NHIA).
  2. Race/ethnicity-stratified 5-year overall survival: Kaplan-Meier
     estimates with log-log confidence intervals, k-group log-rank tests,
     and univariate Cox hazard ratios (Breslow ties).
  3. Age-specific crude incidence trends (2000-2023) and rate ratios
     (2020-2023 vs 2000-2004) with Poisson confidence intervals.
  4. Equivalence assessment of adjusted histology-code hazard ratios
     against a 0.80-1.25 margin, plus Schoenfeld power calculations.

Required private inputs (not distributed; see README):
  raw/gca_case_listing_seer17_2000_2023.csv
  raw/gca_incidence_age_hist_year_seer17_2000_2023.csv
  outputs/formal_next_phase/exclusive_histology_adjusted_cox_models.csv
      (produced by R/01_main_analysis.R)

Outputs are written to outputs/revision_analysis/ (excluded from version
control).

Dependencies: Python 3.x, numpy, pandas, scipy, matplotlib.
"""
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats, optimize

RAW = Path("raw")
OUT = Path("outputs") / "revision_analysis"
OUT.mkdir(parents=True, exist_ok=True)

CASE = RAW / "gca_case_listing_seer17_2000_2023.csv"
INC = RAW / "gca_incidence_age_hist_year_seer17_2000_2023.csv"
COX = Path("outputs") / "formal_next_phase" / "exclusive_histology_adjusted_cox_models.csv"

# ---------------- load D1 ----------------
df = pd.read_csv(CASE, dtype=str)
d1 = df[df["ICD-O-3 Hist/behav, malignant"].str.startswith("8243")].copy()
assert len(d1) == 2571, len(d1)

RACE = "Race recode (with detailed Asian and Native Hawaiian other PI)"
API_LBL = {"Chinese", "Filipino", "Hawaiian", "Asian Indian, Pakistani",
           "Other Asian American", "Korean", "Vietnamese", "Japanese",
           "Laotian", "Samoan", "Guamanian/Chamorro", "Other Pacific Islander"}
def race4(x):
    if x == "White": return "White"
    if x == "Black": return "Black"
    if x in API_LBL: return "Asian/Pacific Islander"
    return "Other/Unknown"   # AIAN + Unknown Race
d1["race4"] = d1[RACE].map(race4)
d1["hisp"] = np.where(d1["Origin recode NHIA (Hispanic, Non-Hisp)"].str.startswith("Spanish"),
                      "Hispanic", "Non-Hispanic")
d1["year"] = d1["Year of diagnosis"].astype(int)
d1["period"] = np.where(d1["year"] <= 2009, "2000-2009", "2010-2023")

# age midpoint for grouped median
def age_lo(a):
    a = a.replace(" years", "").replace("+", "-999")
    return int(a.split("-")[0])
d1["age_lo"] = d1["Age recode with <1 year olds and 90+"].map(age_lo)
AGE_ORDER = sorted(d1["Age recode with <1 year olds and 90+"].unique(), key=age_lo)

def grouped_median(sub):
    vc = sub["Age recode with <1 year olds and 90+"].value_counts()
    n = vc.sum(); half = n / 2.0; cum = 0
    for g in AGE_ORDER:
        c = vc.get(g, 0)
        if cum + c >= half:
            lo = age_lo(g)
            width = 5 if "85" not in g and "90+" not in g else 5
            prev = sum(vc.get(x, 0) for x in AGE_ORDER[:AGE_ORDER.index(g)])
            return lo + (half - prev) / c * width if c else np.nan
        cum += c
    return np.nan

def mar3(x):
    if x in ("Married (including common law)", "Unmarried or Domestic Partner"): return "Married/partnered"
    if x == "Single (never married)": return "Single (never married)"
    if x == "Unknown": return "Unknown"
    return "Divorced/Widowed/Separated"
d1["marital"] = d1["Marital status at diagnosis"].map(mar3)

STG = "Combined Summary Stage with Expanded Regional Codes (2004+)"
def stage3(x):
    if x == "Localized only": return "Localized"
    if x == "Distant site(s)/node(s) involved": return "Distant"
    if x == "Unknown/unstaged/unspecified/DCO": return "Unknown"
    return "Regional"
d1["stage3"] = d1[STG].map(stage3)

# ---------------- Task 1: Table 1 ----------------
def pct(n, N): return f"{n:,} ({100*n/N:.1f}%)"

def table1_block(sub, colname):
    N = len(sub); rows = {}
    rows["N"] = f"{N:,}"
    med = grouped_median(sub)
    rows["Median age, years*"] = f"{med:.1f}"
    nm = (sub["Sex"] == "Male").sum()
    rows["Male sex"] = pct(nm, N)
    for p in ["2000-2009", "2010-2023"]:
        n = (sub["period"] == p).sum()
        rows[f"Period {p}"] = pct(n, N)
    s04 = sub[sub["year"] >= 2004]; n04 = len(s04)
    for s in ["Localized", "Regional", "Distant", "Unknown"]:
        n = (s04["stage3"] == s).sum()
        rows[f"Stage (2004+): {s}"] = f"{n:,} ({100*n/n04:.1f}%)" if n04 else "NA"
    for m in ["Married/partnered", "Single (never married)", "Divorced/Widowed/Separated", "Unknown"]:
        n = (sub["marital"] == m).sum()
        rows[f"Marital: {m}"] = pct(n, N)
    return pd.Series(rows, name=colname)

groups = [("White",), ("Black",), ("Asian/Pacific Islander",), ("Other/Unknown",)]
cols = [table1_block(d1[d1["race4"] == g[0]], g[0]) for g in groups]
cols.append(table1_block(d1, "Total"))
t1_race = pd.concat(cols, axis=1)
t1_race.to_csv(OUT / "task1_table1_by_race.csv", encoding="utf-8-sig")

cols2 = [table1_block(d1[d1["hisp"] == "Hispanic"], "Hispanic"),
         table1_block(d1[d1["hisp"] == "Non-Hispanic"], "Non-Hispanic"),
         table1_block(d1, "Total")]
t1_hisp = pd.concat(cols2, axis=1)
t1_hisp.to_csv(OUT / "task1_table1_by_hispanic.csv", encoding="utf-8-sig")

with open(OUT / "task1_table1.md", "w", encoding="utf-8") as f:
    f.write("**Table 1a. Characteristics of patients with appendiceal goblet cell adenocarcinoma (ICD-O-3 8243/3), SEER-17, 2000-2023, by race.**\n\n")
    f.write(t1_race.to_markdown())
    f.write("\n\n**Table 1b. Same characteristics by Hispanic ethnicity (NHIA).**\n\n")
    f.write(t1_hisp.to_markdown())
    f.write("\n\n\\* Median age approximated by linear interpolation from the 5-year age recode. "
            "Stage percentages use cases diagnosed 2004 or later (Combined Summary Stage available 2004+; "
            "2000-2003 cases excluded from stage rows). 'Other/Unknown' combines American Indian/Alaska Native "
            "and Unknown race. 'Married/partnered' includes common-law marriage and unmarried/domestic partner.\n")
print("Task1 done")
print(t1_race)
print(t1_hisp)

# ---------------- Task 2: survival ----------------
d1["surv_m"] = pd.to_numeric(d1["Survival months"], errors="coerce")
d1["dead"] = (d1["Vital status recode (study cutoff used)"] == "Dead").astype(int)
sv = d1[d1["year"] <= 2018].dropna(subset=["surv_m"]).copy()
# 5-year OS: censor at 60
sv["t"] = sv["surv_m"].clip(upper=60)
sv["e"] = np.where((sv["dead"] == 1) & (sv["surv_m"] <= 60), 1, 0)

def km_at(sub, t_eval=60):
    t = sub["t"].values; e = sub["e"].values
    times = np.sort(np.unique(t[e == 1]))
    S, var_greenwood = 1.0, 0.0
    surv_at = 1.0; var_at = 0.0
    for tt in times:
        n_risk = (t >= tt).sum()
        d = ((t == tt) & (e == 1)).sum()
        if n_risk == 0: continue
        S *= (1 - d / n_risk)
        if n_risk - d > 0:
            var_greenwood += d / (n_risk * (n_risk - d))
        if tt <= t_eval:
            surv_at, var_at = S, var_greenwood
    se = surv_at * math.sqrt(var_at)
    # log-log CI
    if 0 < surv_at < 1 and se > 0:
        theta = math.log(-math.log(surv_at))
        se_theta = se / (surv_at * abs(math.log(surv_at)))
        lo = math.exp(-math.exp(theta + 1.96 * se_theta))
        hi = math.exp(-math.exp(theta - 1.96 * se_theta))
    else:
        lo, hi = np.nan, np.nan
    n_events = int(sub["e"].sum())
    return surv_at, lo, hi, n_events, len(sub)

def logrank(groups_dict):
    # k-group log-rank
    labels = list(groups_dict)
    all_t = np.concatenate([g["t"].values for g in groups_dict.values()])
    all_e = np.concatenate([g["e"].values for g in groups_dict.values()])
    times = np.sort(np.unique(all_t[all_e == 1]))
    k = len(labels)
    O = np.zeros(k); V = np.zeros((k, k))
    for tt in times:
        n = np.array([((g["t"].values) >= tt).sum() for g in groups_dict.values()], float)
        d = np.array([(((g["t"].values) == tt) & ((g["e"].values) == 1)).sum() for g in groups_dict.values()], float)
        N = n.sum(); D = d.sum()
        if N <= 1: continue
        p = n / N
        cov = D * (N - D) / (N - 1) * (np.diag(p) - np.outer(p, p))
        V += cov
        O += d - D * n / N
    idx = list(range(k - 1))
    Ov = O[idx]; Vv = V[np.ix_(idx, idx)]
    chi2 = float(Ov @ np.linalg.pinv(Vv) @ Ov)
    pval = 1 - stats.chi2.cdf(chi2, k - 1)
    return chi2, k - 1, pval

def cox_hr(sub, group_col, ref, cmp):
    """Univariate Cox, Breslow ties. Returns HR, CI, p, n, events."""
    dat = sub[sub[group_col].isin([ref, cmp])].copy()
    x = (dat[group_col] == cmp).astype(float).values
    t = dat["t"].values; e = dat["e"].values
    order = np.argsort(-t)
    t, e, x = t[order], e[order], x[order]
    def nll(beta):
        xb = beta * x
        # log partial likelihood (Breslow)
        ll = 0.0
        ev_t = np.unique(t[e == 1])
        mx = xb.max()
        for tt in ev_t:
            mask_d = (t == tt) & (e == 1)
            risk = t >= tt
            ll += xb[mask_d].sum() - mask_d.sum() * (mx + math.log(np.exp(xb[risk] - mx).sum()))
        return -ll
    def grad(beta):
        xb = beta * x
        g = 0.0
        ev_t = np.unique(t[e == 1])
        mx = xb.max()
        for tt in ev_t:
            mask_d = (t == tt) & (e == 1)
            risk = t >= tt
            w = np.exp(xb[risk] - mx)
            g += x[mask_d].sum() - mask_d.sum() * (w @ x[risk]) / w.sum()
        return -np.array([g])
    res = optimize.minimize(nll, 0.0, jac=grad, method="BFGS")
    beta = res.x[0]
    # observed information via numeric 2nd derivative of nll
    eps = 1e-5
    info = (nll(beta + eps) - 2 * nll(beta) + nll(beta - eps)) / eps**2
    se = math.sqrt(1.0 / info)
    hr = math.exp(beta)
    return hr, math.exp(beta - 1.96 * se), math.exp(beta + 1.96 * se), \
           2 * (1 - stats.norm.cdf(abs(beta / se))), len(dat), int(dat["e"].sum())

surv_rows = []
r3 = {g: sv[sv["race4"] == g] for g in ["White", "Black", "Asian/Pacific Islander"]}
chi2, dfree, pval = logrank(r3)
surv_rows.append({"comparison": "Race: White vs Black vs API", "test": "log-rank",
                  "stat": f"chi2({dfree})={chi2:.3f}", "p": pval})
for g in ["White", "Black", "Asian/Pacific Islander"]:
    s5, lo, hi, ne, n = km_at(r3[g])
    flag = "yes" if ne >= 10 else "NO - <10 events within 5y, do not report"
    surv_rows.append({"comparison": f"Race {g}: 5-year OS (KM)", "test": "KM",
                      "stat": f"{100*s5:.1f}% (95% CI {100*lo:.1f}-{100*hi:.1f})" if not np.isnan(lo) else f"{100*s5:.1f}%",
                      "p": np.nan, "n": n, "events": ne, "reportable": flag})
for cmp in ["Black", "Asian/Pacific Islander"]:
    hr, lo, hi, p, n, ne = cox_hr(sv[sv["race4"].isin(["White", cmp])], "race4", "White", cmp)
    ne_cmp = int(sv[(sv["race4"] == cmp)]["e"].sum())
    flag = "yes" if ne_cmp >= 10 else "NO - <10 events in comparison group, do not report"
    surv_rows.append({"comparison": f"Cox HR {cmp} vs White (unadjusted)", "test": "Cox",
                      "stat": f"HR {hr:.2f} (95% CI {lo:.2f}-{hi:.2f})", "p": p, "n": n, "events": ne,
                      "reportable": flag})
h2 = {"Hispanic": sv[sv["hisp"] == "Hispanic"], "Non-Hispanic": sv[sv["hisp"] == "Non-Hispanic"]}
chi2b, dfb, pb = logrank(h2)
surv_rows.append({"comparison": "Ethnicity: Hispanic vs Non-Hispanic", "test": "log-rank",
                  "stat": f"chi2({dfb})={chi2b:.3f}", "p": pb})
for g in ["Hispanic", "Non-Hispanic"]:
    s5, lo, hi, ne, n = km_at(h2[g])
    flag = "yes" if ne >= 10 else "NO - <10 events within 5y, do not report"
    surv_rows.append({"comparison": f"{g}: 5-year OS (KM)", "test": "KM",
                      "stat": f"{100*s5:.1f}% (95% CI {100*lo:.1f}-{100*hi:.1f})" if not np.isnan(lo) else f"{100*s5:.1f}%",
                      "p": np.nan, "n": n, "events": ne, "reportable": flag})
hr, lo, hi, p, n, ne = cox_hr(sv, "hisp", "Non-Hispanic", "Hispanic")
surv_rows.append({"comparison": "Cox HR Hispanic vs Non-Hispanic (unadjusted)", "test": "Cox",
                  "stat": f"HR {hr:.2f} (95% CI {lo:.2f}-{hi:.2f})", "p": p, "n": n, "events": ne})

t2 = pd.DataFrame(surv_rows)
t2.to_csv(OUT / "task2_survival_by_race_ethnicity.csv", index=False, encoding="utf-8-sig")
print("Task2 done"); print(t2.to_string())

# ---------------- Task 3: age-specific rate trends ----------------
inc = pd.read_csv(INC)
inc = inc[(inc["ICD-O-3 Hist/behav, malignant"].str.startswith("8243"))
          & (inc["Year of diagnosis"] != "2000-2023")].copy()
inc["year"] = inc["Year of diagnosis"].astype(int)
AGE_MAP = {
    "20-24 years": "20-44", "25-29 years": "20-44", "30-34 years": "20-44",
    "35-39 years": "20-44", "40-44 years": "20-44",
    "45-49 years": "45-54", "50-54 years": "45-54",
    "55-59 years": "55-64", "60-64 years": "55-64",
    "65-69 years": "65-74", "70-74 years": "65-74",
    "75-79 years": "75+", "80-84 years": "75+", "85-89 years": "75+", "90+ years": "75+",
}
inc["age_grp"] = inc["Age recode with <1 year olds and 90+"].map(AGE_MAP)
inc = inc.dropna(subset=["age_grp"])
agg = inc.groupby(["age_grp", "year"])[["Count", "Population"]].sum().reset_index()
agg["rate_per_100k"] = agg["Count"] / agg["Population"] * 1e5
wide = agg.pivot(index="year", columns="age_grp", values="rate_per_100k")
wide.to_csv(OUT / "task3_age_specific_rates_annual.csv", encoding="utf-8-sig")
agg.to_csv(OUT / "task3_age_specific_rates_long.csv", index=False, encoding="utf-8-sig")

# rate ratio 2000-2004 vs 2020-2023
rr_rows = []
for g in ["20-44", "45-54", "55-64", "65-74", "75+"]:
    sub = agg[agg["age_grp"] == g]
    e = sub[sub["year"].between(2000, 2004)]
    l = sub[sub["year"].between(2020, 2023)]
    r_e = e["Count"].sum() / e["Population"].sum() * 1e5
    r_l = l["Count"].sum() / l["Population"].sum() * 1e5
    # Poisson CI for RR
    c_e, c_l = e["Count"].sum(), l["Count"].sum()
    se = math.sqrt(1 / c_l + 1 / c_e)
    rr = r_l / r_e
    rr_rows.append({"age_group": g, "rate_2000_2004": round(r_e, 4), "rate_2020_2023": round(r_l, 4),
                    "count_2000_2004": int(c_e), "count_2020_2023": int(c_l),
                    "rate_ratio": round(rr, 2),
                    "rr_lcl": round(math.exp(math.log(rr) - 1.96 * se), 2),
                    "rr_ucl": round(math.exp(math.log(rr) + 1.96 * se), 2)})
t3rr = pd.DataFrame(rr_rows)
t3rr.to_csv(OUT / "task3_rate_ratios_2020_2023_vs_2000_2004.csv", index=False, encoding="utf-8-sig")
print("Task3 done"); print(t3rr.to_string())

# plot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 5))
colors = {"20-44": "#4C72B0", "45-54": "#55A868", "55-64": "#DD8452",
          "65-74": "#C44E52", "75+": "#8172B3"}
for g in ["20-44", "45-54", "55-64", "65-74", "75+"]:
    ax.plot(wide.index, wide[g], marker="o", ms=3, lw=1.6, color=colors[g], label=g)
ax.axvline(2019, color="gray", ls="--", lw=1)
ax.text(2019.2, ax.get_ylim()[1] * 0.92, "2019", color="gray", fontsize=9)
ax.set_xlabel("Year of diagnosis")
ax.set_ylabel("Crude incidence rate per 100,000 person-years")
ax.set_title("Age-specific incidence of appendiceal goblet cell adenocarcinoma (ICD-O-3 8243/3),\nSEER-17, 2000-2023")
ax.legend(title="Age group (years)", frameon=False)
ax.set_xlim(2000, 2023); ax.set_ylim(bottom=0)
fig.tight_layout()
fig.savefig(OUT / "task3_age_specific_rate_trends.png", bbox_inches="tight", dpi=200)
plt.close(fig)
print("plot saved")

# ---------------- Task 4: equivalence / power ----------------
cox = pd.read_csv(COX)
rows4 = []
for ep in ["OS", "CSS"]:
    for term, lbl in [("hist_group8244", "8244 vs 8243"), ("hist_group8245", "8245 vs 8243")]:
        r = cox[(cox["endpoint"] == ep) & (cox["term"] == term)].iloc[0]
        in_margin = (r["LCL_95"] >= 0.80) and (r["UCL_95"] <= 1.25)
        rows4.append({"endpoint": ep, "comparison": lbl, "HR": round(r["HR"], 3),
                      "LCL": round(r["LCL_95"], 3), "UCL": round(r["UCL_95"], 3),
                      "equivalence_0.80_1.25": "YES" if in_margin else "NO"})
t4 = pd.DataFrame(rows4)
# Schoenfeld power: E required = (z_a+z_b)^2 / (p1*p2*(lnHR)^2), p = fraction in smaller group
# Approximation: OS deaths from the full exclusive 2000-2018 cohort (adjusted model used
# the 2004-2018 subset, so true event counts are somewhat smaller and power slightly lower).
za = 1.95996  # alpha=0.05 two-sided
grp = {"8243": {"n": 1324, "os_deaths": 286 + 233 + 10},
       "8244": {"n": 511, "os_deaths": 228 + 51 + 7},
       "8245": {"n": 249, "os_deaths": 57 + 62 + 1}}
power_rows = []
for cmp in ["8244", "8245"]:
    E = grp["8243"]["os_deaths"] + grp[cmp]["os_deaths"]
    p = grp[cmp]["n"] / (grp["8243"]["n"] + grp[cmp]["n"])
    zpower = math.sqrt(E * p * (1 - p)) * abs(math.log(1.25)) - za
    power = stats.norm.cdf(zpower)
    E80 = (za + 0.84162) ** 2 / (p * (1 - p) * math.log(1.25) ** 2)
    E90 = (za + 1.28155) ** 2 / (p * (1 - p) * math.log(1.25) ** 2)
    power_rows.append({"comparison": f"{cmp} vs 8243 (OS)", "os_events_approx": E,
                       "smaller_group_fraction": round(p, 3),
                       "power_to_detect_HR1.25": round(power, 3),
                       "events_needed_80pct": math.ceil(E80),
                       "events_needed_90pct": math.ceil(E90)})
t4p = pd.DataFrame(power_rows)
t4.to_csv(OUT / "task4_equivalence_assessment.csv", index=False, encoding="utf-8-sig")
t4p.to_csv(OUT / "task4_power_schoenfeld.csv", index=False, encoding="utf-8-sig")
print("Task4 done"); print(t4.to_string()); print(t4p.to_string())
print("ALL DONE")
