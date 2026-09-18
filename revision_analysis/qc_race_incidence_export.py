# -*- coding: utf-8 -*-
"""QC of the race-stratified SEER*Stat export against the original benchmark export.

Study: appendiceal goblet cell adenocarcinoma (GCA) case-definition research
in SEER-17, 2000-2023 (revision phase).

Checks:
  A1. Per-histology case totals (single years 2000-2023) in the new
      race-stratified export vs the original benchmark export, and vs the
      new file's combined (all-years) row. The strict definition (ICD-O-3
      8243/3) must total 2571 cases.
  A2. Population consistency: within each year x age x race cell the
      population must be identical across histologies; race-group population
      sums must reproduce the benchmark file's year x age populations.
  A3. Race/origin distribution of D1 (8243) cases.
  A4. Per-year D1 totals: new export vs benchmark export.

Required private inputs (not distributed; see README):
  raw/gca_incidence_age_race_nhia_hist_year_seer17_2000_2023.txt
      (SEER*Stat rate-session export stratified by age, race/origin
      [NHW, NHB, NHAIAN, NHAPI, Hispanic], histology, and year)
  raw/gca_incidence_age_hist_year_seer17_2000_2023.csv
      (original benchmark incidence export)

Outputs are written to outputs/revision_analysis/ (excluded from version
control).

Dependencies: Python 3.x, numpy, pandas.
"""
import re
from pathlib import Path

import numpy as np
import pandas as pd

RAW = Path("raw")
OUT = Path("outputs") / "revision_analysis"
OUT.mkdir(parents=True, exist_ok=True)

NEW = RAW / "gca_incidence_age_race_nhia_hist_year_seer17_2000_2023.txt"
OLD = RAW / "gca_incidence_age_hist_year_seer17_2000_2023.csv"

def hist_code(v):
    m = re.search(r"(?<!\d)(\d{4})(?:/\d)?", str(v))
    return int(m.group(1)) if m else np.nan

def is_single_year(v):
    return bool(re.fullmatch(r"(?:19|20)\d{2}", str(v).strip()))

def to_num(s):
    return pd.to_numeric(s.replace({" ": np.nan, "": np.nan, "~": np.nan}), errors="coerce")

print("=== Loading OLD benchmark file ===", flush=True)
old = pd.read_csv(OLD, dtype=str)
old.columns = [c.strip() for c in old.columns]
print("old columns:", list(old.columns), "rows:", len(old))
old["hist"] = old["ICD-O-3 Hist/behav, malignant"].map(hist_code)
old["count"] = to_num(old["Count"]).fillna(0).astype(int)
old["single_year"] = old["Year of diagnosis"].map(is_single_year)
old_hist_label = old.drop_duplicates("hist")[["hist", "ICD-O-3 Hist/behav, malignant"]]
print("old hist codes n =", len(old_hist_label))

old_tot = old[old["single_year"]].groupby("hist")["count"].sum().rename("old_count_2000_2023")

# old per-year D1
old_yr = old[old["single_year"]].copy()
old_yr["year"] = old_yr["Year of diagnosis"].astype(int)
old_d1_year = old_yr[old_yr["hist"] == 8243].groupby("year")["count"].sum().rename("old_d1")

print("\n=== Loading NEW race file (usecols) ===", flush=True)
usecols = ["ICD-O-3 Hist/behav, malignant", "Year of diagnosis",
           "Age recode with <1 year olds and 90+",
           "Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)",
           "Count", "Pop"]
new = pd.read_csv(NEW, dtype=str, usecols=usecols)
new.columns = ["hist_label", "year_label", "age_label", "race_label", "count_raw", "pop_raw"]
print("new rows:", len(new), flush=True)
new["hist"] = new["hist_label"].map(hist_code)
new["count"] = to_num(new["count_raw"]).fillna(0).astype(int)
new["pop"] = to_num(new["pop_raw"])
new["single_year"] = new["year_label"].map(is_single_year)
print("new hist codes n =", new["hist"].nunique())
print("\nnew year labels:", sorted(new["year_label"].unique()))
print("\nnew age labels:", sorted(new["age_label"].unique()))
print("\nnew race labels:", sorted(new["race_label"].unique()))

# ---- A1: per-histology totals, single years 2000-2023 ----
new_tot = new[new["single_year"]].groupby("hist")["count"].sum().rename("new_count_2000_2023")
# cross-check: combined row in new file (all races incl unknown)
combined = new[~new["single_year"]].groupby("hist")["count"].sum().rename("new_combined_row")
cmp1 = pd.concat([old_hist_label.set_index("hist")["ICD-O-3 Hist/behav, malignant"].rename("label"),
                  old_tot, new_tot, combined], axis=1).fillna(0).astype({"old_count_2000_2023": int, "new_count_2000_2023": int, "new_combined_row": int})
cmp1["match_old_vs_new"] = cmp1["old_count_2000_2023"] == cmp1["new_count_2000_2023"]
cmp1["match_new_vs_combinedrow"] = cmp1["new_count_2000_2023"] == cmp1["new_combined_row"]
print("\n=== A1 histology totals comparison (nonzero rows) ===")
print(cmp1[(cmp1["old_count_2000_2023"] > 0) | (cmp1["new_count_2000_2023"] > 0)].to_string())
cmp1.to_csv(OUT / "qc_a1_histology_totals.csv")
all_match = bool(cmp1["match_old_vs_new"].all())
n_codes = len(cmp1)
d1_total = int(cmp1.loc[8243, "new_count_2000_2023"]) if 8243 in cmp1.index else None
print(f"\nA1 RESULT: {n_codes} histology codes, all match old = {all_match}; D1(8243) total = {d1_total} (expected 2571)")
mism = cmp1[~cmp1["match_old_vs_new"]]
print(f"A1 mismatched codes: {len(mism)}")
if len(mism):
    print(mism.to_string())

# ---- A2: population consistency within year x age x race across histologies ----
print("\n=== A2 population consistency check ===", flush=True)
sub = new[new["single_year"]].copy()
g = sub.groupby(["year_label", "age_label", "race_label"])["pop"].nunique()
bad = g[g > 1]
print(f"year x age x race cells checked: {len(g)}; cells with >1 distinct Pop across histologies: {len(bad)}")
if len(bad):
    print(bad.head(20))
# also verify new-file pop sums (5 race groups, excl unknown) reproduce old-file population per year x age
def age_key(v):
    t = re.sub(r"(?i)years?|yrs?", "", str(v)).strip().lower()
    t = re.sub(r"[_\s]+", "-", t).strip("-")
    if t in {"0", "00", "<1"}: return "00"
    if t.startswith("90"): return "90+"
    if t == "unknown": return "Unknown"
    nums = re.findall(r"\d+", t)
    if len(nums) >= 2: return f"{int(nums[0]):02d}-{int(nums[1]):02d}"
    return t
new["age"] = new["age_label"].map(age_key)
new["year"] = pd.to_numeric(new["year_label"], errors="coerce").fillna(-1).astype(int)
old["age"] = old["Age recode with <1 year olds and 90+"].map(age_key)
old["year"] = pd.to_numeric(old["Year of diagnosis"], errors="coerce").fillna(-1).astype(int)
old["pop"] = to_num(old["Population"])

# note: pop within a race x year x age cell is identical across hist -> take max per race first
new_pop_cell = (new[(new["single_year"]) & (new["age"] != "Unknown") & (new["pop"].notna())]
                .groupby(["year", "age", "race_label"])["pop"].max()
                .groupby(["year", "age"]).sum().rename("new_pop_sumraces"))
old_pop_cell = (old[(old["single_year"]) & (old["age"] != "Unknown")]
                .groupby(["year", "age"])["pop"].max().rename("old_pop"))
popcmp = pd.concat([old_pop_cell, new_pop_cell], axis=1)
popcmp["diff"] = popcmp["new_pop_sumraces"] - popcmp["old_pop"]
n_bad = int((popcmp["diff"].fillna(1) != 0).sum())
print(f"year x age cells pop compare (new summed over races incl unknown vs old): {len(popcmp)} cells, mismatched: {n_bad}")
print(popcmp[popcmp["diff"].fillna(1) != 0].head(15).to_string())
popcmp.to_csv(OUT / "qc_a2_population_cells.csv")

# ---- A3: race label counts distribution (D1) ----
print("\n=== A3 race distribution (D1=8243, 2000-2023 single years) ===")
d1 = new[(new["single_year"]) & (new["hist"] == 8243)]
race_dist = d1.groupby("race_label")["count"].sum().rename("d1_cases")
race_dist_all = new[new["single_year"]].groupby("race_label")["count"].sum().rename("all_hist_cases")
a3 = pd.concat([race_dist, race_dist_all], axis=1).fillna(0).astype(int)
a3.to_csv(OUT / "qc_a3_race_distribution.csv")
print(a3.to_string())

# ---- A4: per-year D1 totals new vs old ----
print("\n=== A4 per-year D1 totals ===")
new_d1_year = d1.groupby("year")["count"].sum().rename("new_d1")
a4 = pd.concat([old_d1_year, new_d1_year], axis=1).fillna(0).astype(int)
a4["match"] = a4["old_d1"] == a4["new_d1"]
a4.to_csv(OUT / "qc_a4_d1_by_year.csv")
print(a4.to_string())
print(f"A4 RESULT: all years match = {bool(a4['match'].all())}; total new D1 = {a4['new_d1'].sum()}")

qc_pass = all_match and (d1_total == 2571) and (len(bad) == 0) and bool(a4["match"].all())
print(f"\n======= QC OVERALL: {'PASS' if qc_pass else 'FAIL'} =======")
with open(OUT / "qc_overall.txt", "w", encoding="utf-8") as f:
    f.write(f"QC_PASS={qc_pass}\nA1_all_match={all_match} n_codes={n_codes} D1_total={d1_total}\n"
            f"A2_bad_cells={len(bad)} A2_pop_mismatch_cells={n_bad}\nA4_all_years_match={bool(a4['match'].all())}\n")
