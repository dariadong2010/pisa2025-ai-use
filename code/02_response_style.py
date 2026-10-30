"""Response-style indices from 43 four-point agreement items in 10 scales that do not
concern AI use or the four disposition indices.

ERS          share of answers in an end category (1 or 4), all items
ERS_pos      ERS on the 27 positively worded items
ERS_neg      ERS on the 9 negatively worded items
ERS_neu      ERS on the 7 content-neutral items
ERS_norel    ERS without the scales whose content may relate to the dispositions
             (ST131 epistemic beliefs, ST324 attitudes towards school, ST499/ST500
             environmental efficacy), 26 items
ERS_ict      ERS on the 13 ICT items (IC179 views on ICT rules, IC172 school ICT conditions)
EXT_agree    share of "strongly agree" answers, all items
EXT_disagree share of "strongly disagree" answers, all items
ARS_pos / ARS_neg   share of agreement with positively / negatively worded items
ARS_bal      mean of ARS_pos and ARS_neg (balanced acquiescence)

Writes output/derived/response_style.parquet.
"""
import pandas as pd
import pyreadstat
from common import DERIVED, RESULTS, pisa_file

GRIDS = ["ST034", "ST265", "ST131", "ST324", "ST420", "ST493", "IC179", "IC172", "ST499", "ST500"]
REVERSED = ("ST034", "ST265")          # coded 1 = strongly agree in the public-use file
NEGATIVE = ["ST034Q01TA", "ST034Q04TA", "ST034Q06TA", "ST324Q10JA", "ST324Q11JA",
            "ST420Q01DA", "ST420Q02DA", "ST493Q01DA", "ST131Q17DA"]
NEUTRAL = [f"IC179Q0{i}JA" for i in range(1, 7)] + ["ST131Q05DA"]   # enter only the ERS indices
RELATED = ("ST131", "ST324", "ST499", "ST500")
ICT = ("IC179", "IC172")

stu_file = str(pisa_file("STU_FILE"))
_, meta = pyreadstat.read_sav(stu_file, metadataonly=True)
items = [c for c in meta.column_names if c[:5] in GRIDS and c[5] == "Q"]
d, _ = pyreadstat.read_sav(stu_file, usecols=["CNT", "CNTSTUID"] + items)
d = d[d["CNT"] == "QCI"].copy()

x = d[items].where((d[items] >= 1) & (d[items] <= 4))
for c in items:
    if c.startswith(REVERSED):
        x[c] = 5 - x[c]                # 1 = strongly disagree ... 4 = strongly agree
positive = [c for c in items if c not in NEGATIVE + NEUTRAL]

out = pd.DataFrame({"CNTSTUID": d["CNTSTUID"].values})
out["rs_n_items"] = x.notna().sum(axis=1).values
out["ERS"] = x.isin([1, 4]).sum(axis=1).values / out["rs_n_items"].where(out["rs_n_items"] >= 20)
agree = x.ge(3).astype(float).where(x.notna())
out["ARS_pos"] = agree[positive].mean(axis=1).values
out["ARS_neg"] = agree[NEGATIVE].mean(axis=1).values
out["ARS_bal"] = (out["ARS_pos"] + out["ARS_neg"]) / 2

sets = {
    "ERS_pos": positive, "ERS_neg": NEGATIVE, "ERS_neu": NEUTRAL,
    "ERS_norel": [c for c in items if not c.startswith(RELATED)],
    "ERS_ict": [c for c in items if c.startswith(ICT)],
}
for name, cols in sets.items():
    n = x[cols].notna().sum(axis=1)
    out[name] = (x[cols].isin([1, 4]).sum(axis=1) / n.where(n >= max(3, len(cols) // 2))).values
n_all = x.notna().sum(axis=1).where(lambda s: s >= 20)
out["EXT_agree"] = (x.eq(4).sum(axis=1) / n_all).values
out["EXT_disagree"] = (x.eq(1).sum(axis=1) / n_all).values
out.to_parquet(DERIVED / "response_style.parquet", index=False)

# Cronbach's alpha of the extreme-answer indicators (1 = end category), complete cases
ext = x.isin([1, 4]).astype(float).where(x.notna())


def alpha(cols):
    s = ext[cols].dropna()
    k = len(cols)
    return k / (k - 1) * (1 - s.var().sum() / s.sum(axis=1).var())


rel = pd.DataFrame([{"index": "ERS", "n_items": len(items), "alpha": alpha(items)}] +
                   [{"index": k, "n_items": len(v), "alpha": alpha(v)} for k, v in sets.items()])
rel.to_csv(RESULTS / "response_style_alpha.csv", index=False)
print(rel.round(2).to_string(index=False))

print(f"items {len(items)} (positive {len(positive)}, negative {len(NEGATIVE)}, neutral {len(NEUTRAL)})")
print({k: len(v) for k, v in sets.items()})
print(out.drop(columns="CNTSTUID").describe().round(3).to_string())
print("correlations:\n", out.drop(columns=["CNTSTUID", "rs_n_items"]).corr().round(2).to_string())
