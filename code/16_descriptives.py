"""Weighted descriptive statistics (Table 1), weighted correlations (Table S1) and the list
of response-style items (Table S6).

Writes descriptives.csv, correlations.csv and response_style_items.csv to output/results.
Statistics use the final student weight; statistics involving science are averaged over
the 10 plausible values. Correlations use pairwise deletion.
"""
from itertools import combinations
import numpy as np
import pandas as pd
import pyreadstat
from common import DERIVED, RESULTS, pisa_file

d = pd.read_parquet(DERIVED / "analysis_sample.parquet")
d = d[d["lca_in"] == 1].copy()
d = d.merge(pd.read_parquet(DERIVED / "response_style.parquet"), on="CNTSTUID", how="left")
d["grade_below"] = (d["GRADE"] < 0).astype(float)
d["eval_ai"] = d["class_eval_ai"] - 1
d["private"] = (d["PRIVATESCH"] == 2).astype(float)
w = d["W_FSTUWT"].to_numpy()

VARS = ["ST438Q01DA", "ST438Q02DA", "ST438Q03DA", "ST438Q04DA", "ESCS", "female", "AGE", "grade_below",
        "eval_ai", "ai_school_avail", "ai_home_avail", "private", "CURIO", "PERSEV", "COGABIL", "ENPROBS",
        "ERS", "ERS_ict", "ERS_neu", "ARS_bal", "SCIE"]
CORR_VARS = [v for v in VARS if v not in ("AGE", "ERS_ict", "ERS_neu")]


def col(v, m=1):
    return d[f"PV{m}SCIE"].to_numpy(float) if v == "SCIE" else d[v].to_numpy(float)


def wstats(x):
    ok = ~np.isnan(x)
    m = np.average(x[ok], weights=w[ok])
    return m, np.sqrt(np.average((x[ok] - m) ** 2, weights=w[ok])), 100 * (1 - ok.mean())


rows = []
for v in VARS:
    st = np.array([wstats(col(v, m)) for m in range(1, 11)]).mean(0) if v == "SCIE" else wstats(col(v))
    rows.append({"variable": v, "mean": st[0], "sd": st[1], "missing_pct": st[2]})
desc = pd.DataFrame(rows)
desc.to_csv(RESULTS / "descriptives.csv", index=False)
print(desc.round(3).to_string(index=False))


def wcorr(x, y):
    ok = ~np.isnan(x) & ~np.isnan(y)
    x, y, ww = x[ok], y[ok], w[ok]
    mx, my = np.average(x, weights=ww), np.average(y, weights=ww)
    return np.average((x - mx) * (y - my), weights=ww) / np.sqrt(
        np.average((x - mx) ** 2, weights=ww) * np.average((y - my) ** 2, weights=ww))


k = len(CORR_VARS)
corr = np.eye(k)
for i, j in combinations(range(k), 2):
    vi, vj = CORR_VARS[i], CORR_VARS[j]
    if "SCIE" in (vi, vj):
        r = np.mean([wcorr(col(vi, m), col(vj, m)) for m in range(1, 11)])
    else:
        r = wcorr(col(vi), col(vj))
    corr[i, j] = corr[j, i] = r
corr = pd.DataFrame(corr, index=CORR_VARS, columns=CORR_VARS)
corr.to_csv(RESULTS / "correlations.csv")
print("\nweighted correlations:\n", corr.round(2).to_string())

# Items of the response-style indices with wording, keying and item sets
GRIDS = {"ST034": "Sense of belonging", "ST265": "Feeling safe", "ST131": "Epistemic beliefs about science",
         "ST324": "Attitudes towards school", "ST420": "Fixed beliefs about ability", "ST493": "Fixed beliefs about status",
         "IC179": "Views on ICT rules", "IC172": "School ICT conditions", "ST499": "Environmental efficacy (community)",
         "ST500": "Environmental efficacy (personal)"}
NEGATIVE = ["ST034Q01TA", "ST034Q04TA", "ST034Q06TA", "ST324Q10JA", "ST324Q11JA",
            "ST420Q01DA", "ST420Q02DA", "ST493Q01DA", "ST131Q17DA"]
NEUTRAL = [f"IC179Q0{i}JA" for i in range(1, 7)] + ["ST131Q05DA"]
_, meta = pyreadstat.read_sav(str(pisa_file("STU_FILE")), metadataonly=True)
items = []
for c in meta.column_names:
    if c[:5] in GRIDS and c[5] == "Q":
        sets = ["All"] + (["ICT"] if c[:5] in ("IC179", "IC172") else []) + \
               ([] if c[:5] in ("ST131", "ST324", "ST499", "ST500") else ["Excl. related"])
        items.append({"item": c, "scale": GRIDS[c[:5]],
                      "wording": meta.column_names_to_labels[c].split(" - ")[-1].replace("<", "").replace(">", ""),
                      "keying": "negative" if c in NEGATIVE else "neutral" if c in NEUTRAL else "positive",
                      "reverse_coded": c[:5] in ("ST034", "ST265"), "sets": ", ".join(sets)})
pd.DataFrame(items).to_csv(RESULTS / "response_style_items.csv", index=False)
print(f"\nresponse-style items: {len(items)}")
