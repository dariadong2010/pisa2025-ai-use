"""Analysis file for the B-S-J-Z (China) sample of PISA 2025.

Reads the student, school, questionnaire-timing and LDW item files and writes
output/derived/analysis_sample.parquet (one row per student, CNT == "QCI").
"""
import pandas as pd
import pyreadstat
from common import DERIVED, pisa_file

USE = [f"ST438Q0{i}DA" for i in range(1, 5)]
CURIO_ITEMS = ["ST301Q01JA", "ST301Q04JA", "ST301Q06JA", "ST301Q08JA", "ST301Q10JA"]
PERSEV_ITEMS = ["ST307Q01JA", "ST307Q02JA", "ST307Q03JA", "ST307Q07JA", "ST307Q09JA", "ST307Q10JA"]
PV = [f"PV{i}{d}" for d in ("SCIE", "CMPS", "READ", "MATH") for i in range(1, 11)]
RWT = [f"W_FSTURWT{i}" for i in range(1, 81)]
STU_VARS = (
    ["CNT", "CNTSCHID", "CNTSTUID", "STRATUM", "BOOKID", "W_FSTUWT", "MALE", "GRADE", "AGE", "ESCS",
     "CURIO", "PERSEV", "COGABIL", "ENPROBS", "AIUSESCH",
     "ST436Q16DA", "IC170Q10DA", "IC171Q10DA", "ST331Q03JA"]
    + USE + CURIO_ITEMS + PERSEV_ITEMS + PV + RWT
)
SCH_VARS = ["CNT", "CNTSCHID", "W_NRASCHBWT", "PRIVATESCH", "SC001Q01TA"]

stu, _ = pyreadstat.read_sav(str(pisa_file("STU_FILE")), usecols=STU_VARS)
stu = stu[stu["CNT"] == "QCI"].copy()
sch, _ = pyreadstat.read_sav(str(pisa_file("SCH_FILE")), usecols=SCH_VARS)
sch = sch[sch["CNT"] == "QCI"].drop(columns="CNT")
df = stu.merge(sch, on="CNTSCHID", how="left")

# Time spent on the ST438 screen (the timing file has to be read as latin1)
tt, _ = pyreadstat.read_sav(str(pisa_file("TT_FILE")), usecols=["CNT", "CNTSTUID", "ST438_TT"],
                            encoding="latin1")
tt = tt[tt["CNT"] == "QCI"].drop(columns="CNT")
df = df.merge(tt, on="CNTSTUID", how="left")

# Students who took the learning-in-the-digital-world (LDW) domain: at least one scored item
ldw, _ = pyreadstat.read_sav(str(pisa_file("LDW_FILE")))
ldw = ldw[ldw["CNT"] == "QCI"]
scored = [c for c in ldw.columns if c.endswith("S") and c[:2] in ("P1", "P2", "P3", "P4")]
taker_ids = ldw.loc[ldw[scored].notna().any(axis=1), "CNTSTUID"]
df["ldw_taker"] = df["CNTSTUID"].isin(taker_ids).astype(int)

# Latent class indicators keep their 1-5 coding
for v in USE:
    df[v] = df[v].where(df[v].between(1, 5))
n_use = df[USE].notna().sum(axis=1)
df["lca_in"] = (n_use > 0).astype(int)

# Response patterns used to flag possibly careless respondents
same4 = df[USE].nunique(axis=1).eq(1) & n_use.eq(4)
df["use_all_never"] = (same4 & df["ST438Q01DA"].eq(1)).astype(int)
df["use_same_nonnever"] = (same4 & df["ST438Q01DA"].gt(1)).astype(int)
for stem, cols in (("curio", CURIO_ITEMS), ("persev", PERSEV_ITEMS)):
    complete = df[cols].notna().all(axis=1)
    df[f"{stem}_straight"] = (complete & df[cols].nunique(axis=1).eq(1)).astype(float).where(complete)

# Correlates
df["female"] = 1 - df["MALE"]
df["class_eval_ai"] = df["ST436Q16DA"]
for src, stem in (("IC170Q10DA", "ai_school"), ("IC171Q10DA", "ai_home")):
    df[f"{stem}_avail"] = (df[src] != 6).astype(float).where(df[src].notna())
df["q_effort"] = df["ST331Q03JA"]

df.to_parquet(DERIVED / "analysis_sample.parquet", index=False)

print(f"students {len(df)}, schools {df['CNTSCHID'].nunique()}, in latent class sample {df['lca_in'].sum()}")
print(f"LDW takers {df['ldw_taker'].sum()}; by booklet:")
print(df.groupby("BOOKID")["ldw_taker"].agg(["size", "sum"]).to_string())
print(f"all never {df['use_all_never'].mean():.3f}, same answer other than never {df['use_same_nonnever'].mean():.3f}")
print(f"identical answers ST301 {df['curio_straight'].mean():.3f}, ST307 {df['persev_straight'].mean():.3f}")
chk = ["ESCS", "female", "GRADE", "AGE", "class_eval_ai", "ai_school_avail", "ai_home_avail",
       "PRIVATESCH", "SC001Q01TA", "CURIO", "PERSEV", "COGABIL", "ENPROBS", "q_effort", "ST438_TT"]
print(df[chk].describe().T[["count", "mean", "min", "max"]].round(2).to_string())
