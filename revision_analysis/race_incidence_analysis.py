# -*- coding: utf-8 -*-
"""Race/origin-specific age-standardized incidence of appendiceal GCA, SEER-17 2000-2023.

Study: appendiceal goblet cell adenocarcinoma (GCA) case-definition research
in SEER (revision phase).

Direct standardization to the 2000 US Standard Population; Fay-Feuer gamma
confidence intervals (same method as R/01_main_analysis.R).

Analyses:
  B1. D1 (strict, ICD-O-3 8243/3) age-standardized rate by race/origin,
      2000-2023.
  B2. D1 rates in 2000-2009 vs 2010-2023 and rate ratios with CIs.
  B3. Annual D1 age-standardized rates by race/origin (CSV; figure in
      race_incidence_figure.py).
  B4. D1/D2/D3 rates by race/origin, 2000-2023.

Required private input (not distributed; see README):
  raw/gca_incidence_age_race_nhia_hist_year_seer17_2000_2023.txt
      (SEER*Stat rate-session export stratified by age, race/origin
      [NHW, NHB, NHAIAN, NHAPI, Hispanic], histology, and year;
      run qc_race_incidence_export.py first to validate this export)

Outputs are written to outputs/revision_analysis/ (excluded from version
control).

Dependencies: Python 3.x, numpy, pandas, scipy.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import gamma as gamma_dist

RAW = Path("raw")
OUT = Path("outputs") / "revision_analysis"
OUT.mkdir(parents=True, exist_ok=True)

NEW = RAW / "gca_incidence_age_race_nhia_hist_year_seer17_2000_2023.txt"

STD_POP = {"00": 3794901, "01-04": 15191619, "05-09": 19919840, "10-14": 20056779,
           "15-19": 19819518, "20-24": 18257225, "25-29": 17722067, "30-34": 19511370,
           "35-39": 22179956, "40-44": 22479229, "45-49": 19805793, "50-54": 17224359,
           "55-59": 13307234, "60-64": 10654272, "65-69": 9409940, "70-74": 8725574,
           "75-79": 7414559, "80-84": 4900234, "85-89": 2678567, "90+": 1580606}
AGE_ORDER = list(STD_POP.keys())
STD_W = {k: v / sum(STD_POP.values()) for k, v in STD_POP.items()}

RACE_MAP = {
    "Non_Hispanic_White": "NHW",
    "Non_Hispanic_Black": "NHB",
    "Non_Hispanic_American_IndianAlaska_Native": "NHAIAN",
    "Non_Hispanic_Asian_or_Pacific_Islander": "NHAPI",
    "Hispanic_All_Races": "Hispanic",
}
RACE_ORDER = ["NHW", "NHB", "NHAIAN", "NHAPI", "Hispanic"]
DEFS = {"D1": [8243], "D2": [8243, 8245], "D3": [8243, 8244, 8245]}


def age_key(v):
    t = re.sub(r"(?i)years?|yrs?", "", str(v)).strip().lower()
    t = re.sub(r"[_\s]+", "-", t).strip("-")
    if t in {"0", "00", "<1"}:
        return "00"
    if t.startswith("90"):
        return "90+"
    nums = re.findall(r"\d+", t)
    if len(nums) >= 2:
        return f"{int(nums[0]):02d}-{int(nums[1]):02d}"
    return None


def fay_feuer(cases, pops, weights, alpha=0.05):
    """cases, pops, weights aligned arrays. Returns rate, se, lcl, ucl per 100k."""
    d = np.asarray(cases, float)
    n = np.asarray(pops, float)
    w = np.asarray(weights, float)
    ok = np.isfinite(d) & np.isfinite(n) & (n > 0)
    d, n, w = d[ok], n[ok], w[ok]
    r = float(np.sum(w * d / n))
    v = float(np.sum(w**2 * d / n**2))
    wm = float(np.max(w / n))
    if r <= 0 or v <= 0:
        return 0.0, 0.0, 0.0, -np.log(alpha) * wm * 1e5
    lcl = gamma_dist.ppf(alpha / 2, a=r**2 / v, scale=v / r)
    ucl = gamma_dist.ppf(1 - alpha / 2, a=(r + wm) ** 2 / (v + wm**2),
                         scale=(v + wm**2) / (r + wm))
    return r * 1e5, np.sqrt(v) * 1e5, lcl * 1e5, ucl * 1e5


print("loading new file ...", flush=True)
use = ["ICD-O-3 Hist/behav, malignant", "Year of diagnosis",
       "Age recode with <1 year olds and 90+",
       "Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)", "Count", "Pop"]
df = pd.read_csv(NEW, dtype=str, usecols=use)
df.columns = ["hist_label", "year_label", "age_label", "race_label", "count_raw", "pop_raw"]
df["hist"] = pd.to_numeric(df["hist_label"].str.extract(r"(?<!\d)(\d{4})")[0], errors="coerce")
df["year"] = pd.to_numeric(df["year_label"], errors="coerce")
df["age"] = df["age_label"].map(age_key)
df["cases"] = pd.to_numeric(df["count_raw"].replace({" ": np.nan}), errors="coerce").fillna(0).astype(int)
df["pop"] = pd.to_numeric(df["pop_raw"].replace({" ": np.nan}), errors="coerce")
df["race"] = df["race_label"].map(RACE_MAP)
df = df[df["year"].between(2000, 2023) & df["age"].isin(AGE_ORDER) & df["race"].notna()]
print("analysis rows:", len(df), flush=True)

# population is constant across histology within year x age x race -> dedupe once
pop_cells = (df.drop_duplicates(["year", "age", "race"])[["year", "age", "race", "pop"]])

def period_cells(codes, y0, y1):
    sub = df[df["hist"].isin(codes) & df["year"].between(y0, y1)]
    cases = sub.groupby(["race", "age"])["cases"].sum()
    pops = pop_cells[pop_cells["year"].between(y0, y1)].groupby(["race", "age"])["pop"].sum()
    out = pd.DataFrame({"cases": cases, "pop": pops}).fillna(0)
    return out

def rate_table(codes, y0, y1):
    cells = period_cells(codes, y0, y1)
    rows = []
    for race in RACE_ORDER:
        cases = [cells["cases"].get((race, a), 0) for a in AGE_ORDER]
        pops = [cells["pop"].get((race, a), 0) for a in AGE_ORDER]
        w = [STD_W[a] for a in AGE_ORDER]
        r, se, l, u = fay_feuer(cases, pops, w)
        rows.append({"race": race, "cases": int(sum(cases)), "person_years": int(sum(pops)),
                     "asr_per_100k": r, "se": se, "lcl95": l, "ucl95": u})
    return pd.DataFrame(rows)

# ---- sanity: pooled all-race D1 2000-2023 should reproduce 0.1208 ----
cells_all = period_cells(DEFS["D1"], 2000, 2023).groupby("age").sum()
r_all = fay_feuer([cells_all["cases"].get(a, 0) for a in AGE_ORDER],
                  [cells_all["pop"].get(a, 0) for a in AGE_ORDER],
                  [STD_W[a] for a in AGE_ORDER])
print(f"sanity pooled D1 2000-2023 (excl unknown race): cases={int(cells_all['cases'].sum())}, ASR={r_all[0]:.4f} (benchmark 0.1208, cases 2571 incl 11 unknown-race)")

# ---- B1: D1 2000-2023 by race ----
b1 = rate_table(DEFS["D1"], 2000, 2023)
b1.insert(0, "definition", "D1")
b1.to_csv(OUT / "b1_d1_asr_by_race_2000_2023.csv", index=False)
print("\nB1 D1 2000-2023 ASR by race:"); print(b1.to_string(index=False))

# ---- B2: two periods + RR ----
p1 = rate_table(DEFS["D1"], 2000, 2009).add_suffix("_p1").rename(columns={"race_p1": "race"})
p2 = rate_table(DEFS["D1"], 2010, 2023).add_suffix("_p2").rename(columns={"race_p2": "race"})
b2 = p1.merge(p2, on="race")
rr = b2["asr_per_100k_p2"] / b2["asr_per_100k_p1"]
se_log = np.sqrt((b2["se_p1"] / b2["asr_per_100k_p1"]) ** 2 + (b2["se_p2"] / b2["asr_per_100k_p2"]) ** 2)
b2["rr_p2_vs_p1"] = rr
b2["rr_lcl95"] = np.exp(np.log(rr) - 1.96 * se_log)
b2["rr_ucl95"] = np.exp(np.log(rr) + 1.96 * se_log)
b2.to_csv(OUT / "b2_d1_asr_by_race_periods_rr.csv", index=False)
print("\nB2 D1 periods 2000-2009 vs 2010-2023:"); print(b2.to_string(index=False))

# ---- B3: annual D1 ASR by race (csv; figure in separate script) ----
ann_rows = []
for race in RACE_ORDER:
    for yr in range(2000, 2024):
        sub = df[(df["hist"].isin(DEFS["D1"])) & (df["year"] == yr) & (df["race"] == race)]
        pc = pop_cells[(pop_cells["year"] == yr) & (pop_cells["race"] == race)]
        cases = [sub.loc[sub["age"] == a, "cases"].sum() for a in AGE_ORDER]
        pops = [pc.loc[pc["age"] == a, "pop"].sum() for a in AGE_ORDER]
        r, se, l, u = fay_feuer(cases, pops, [STD_W[a] for a in AGE_ORDER])
        ann_rows.append({"race": race, "year": yr, "cases": int(sum(cases)),
                         "asr_per_100k": r, "lcl95": l, "ucl95": u})
ann = pd.DataFrame(ann_rows)
ann.to_csv(OUT / "b3_d1_asr_annual_by_race.csv", index=False)
print("\nB3 annual rows:", len(ann))

# ---- B4: D2, D3 combined 2000-2023 ----
frames = []
for dname in ["D1", "D2", "D3"]:
    t = rate_table(DEFS[dname], 2000, 2023)
    t.insert(0, "definition", dname)
    frames.append(t)
b4 = pd.concat(frames, ignore_index=True)
b4.to_csv(OUT / "b4_d1d2d3_asr_by_race_2000_2023.csv", index=False)
print("\nB4 D1/D2/D3 2000-2023 by race:"); print(b4.to_string(index=False))
print("\nDONE", flush=True)
