"""Run the complete analysis in order and write a log for every step.

Usage (from the repository root):
    python run_all.py                 # all steps
    python run_all.py --from 05       # start at step 05
    python run_all.py --only 15 16
Rscript is taken from the RSCRIPT environment variable or from the PATH.
"""
import argparse
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / "code"))
from common import OUTPUT  # noqa: E402

STEPS = [
    ("01", "code/01_prepare_data.py"),
    ("02", "code/02_response_style.py"),
    ("03", "code/03_enumerate_classes.R"),
    ("04", "code/04_careless_responding.R"),
    ("05", "code/05_correlates_three_step.R"),
    ("06", "code/06_correlates_imputation.R"),
    ("07", "code/07_outcomes_bch.R"),
    ("08", "code/08_outcomes_composite.R"),
    ("09", "code/09_response_style_bch.R"),
    ("10", "code/10_unweighted_multilevel.R"),
    ("11", "code/11_within_school.R"),
    ("12", "code/12_school_composition.R"),
    ("13", "code/13_adjacent_solutions.R"),
    ("14", "code/14_continuous_and_ers_correlates.R"),
    ("15", "code/15_local_dependence.R"),
    ("16", "code/16_descriptives.py"),
]


def command(script):
    if script.endswith(".R"):
        rscript = os.environ.get("RSCRIPT") or shutil.which("Rscript")
        if not rscript:
            sys.exit("Rscript not found; set the RSCRIPT environment variable")
        return [rscript, script]
    return [sys.executable, script]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--from", dest="start", default=None)
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()
    ids = [s for s, _ in STEPS]
    todo = STEPS
    if args.only:
        todo = [s for s in STEPS if s[0] in args.only]
    elif args.start:
        todo = STEPS[ids.index(args.start):]
    logs = OUTPUT / "logs"
    logs.mkdir(parents=True, exist_ok=True)
    for step, scripts in todo:
        for script in (scripts if isinstance(scripts, list) else [scripts]):
            log = logs / (Path(script).stem + ".log")
            t0 = time.time()
            print(f"[{step}] {script} ...", flush=True)
            with open(log, "w", encoding="utf-8") as fh:
                rc = subprocess.run(command(script), cwd=ROOT, stdout=fh, stderr=subprocess.STDOUT).returncode
            print(f"      finished in {time.time() - t0:.0f} s (log: {log})", flush=True)
            if rc != 0:
                sys.exit(f"{script} failed with exit code {rc}; see {log}")


if __name__ == "__main__":
    main()
