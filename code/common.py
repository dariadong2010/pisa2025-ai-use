"""Paths shared by the Python scripts. All scripts are run from the repository root."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def read_config(file=ROOT / "config.txt"):
    cfg = {}
    for line in file.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#") and "=" in line:
            key, value = line.split("=", 1)
            cfg[key.strip()] = value.strip()
    return cfg


CFG = read_config()
PISA_DIR = (ROOT / CFG["PISA_DIR"]).resolve()
OUTPUT = (ROOT / CFG.get("OUTPUT_DIR", "output")).resolve()
DERIVED = OUTPUT / "derived"
RESULTS = OUTPUT / "results"
for p in (DERIVED, RESULTS):
    p.mkdir(parents=True, exist_ok=True)


def pisa_file(key):
    """Full path of one of the PISA files named in config.txt."""
    return (PISA_DIR / CFG[key]).resolve()


CLASSES = ["non-user", "occasional", "learning-focused", "moderate all-purpose",
           "frequent all-purpose", "daily all-purpose"]
