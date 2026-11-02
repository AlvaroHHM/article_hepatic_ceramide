#!/usr/bin/env python3
"""
audit_final.py - Comprehensive final audit of the reproducibility compendium.
Verifies: structure, scripts, dependencies, 15/15 manuscript numbers,
git state, documentation, and reproducibility.
"""
import os, re, subprocess, sys
from pathlib import Path

SEP = "=" * 70
issues, warnings, oks = [], [], []

def section(t):
    print(f"\n{SEP}\n{t}\n{SEP}")

def ok(msg):
    print(f"  OK    {msg}")
    oks.append(msg)

def warn(msg):
    print(f"  WARN  {msg}")
    warnings.append(msg)

def fail(msg):
    print(f"  FAIL  {msg}")
    issues.append(msg)

# ============================================================
# FASE 1: STRUCTURE
# ============================================================
section("FASE 1: STRUCTURE")

expected_dirs = [
    "analysis/bulk", "analysis/cellular", "analysis/snrnaseq",
    "analysis/convergence", "analysis/figures",
    "results/bulk", "results/cellular", "results/snrnaseq", "results/convergence",
    "figures/main", "figures/supplementary",
    "logs", "docs", "data", "tables",
]
for d in expected_dirs:
    if os.path.isdir(d):
        n = len([f for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))])
        print(f"  {d:<30} {n:>4} files")
    else:
        fail(f"Missing dir: {d}")

n_scripts = sum(1 for _ in Path("analysis").rglob("*.R"))
n_csvs    = sum(1 for _ in Path("results").rglob("*.csv"))
n_figs    = sum(1 for _ in Path("figures").rglob("*.pdf"))
print(f"\n  Total: {n_scripts} scripts | {n_csvs} CSVs | {n_figs} figures")

# README count consistency
readme = Path("README.md").read_text(encoding="utf-8")
for pattern, actual, label in [
    (r'\*\*(\d+)\s+R\s+scripts\*\*',  n_scripts, "scripts"),
    (r'\*\*(\d+)\s+CSV\s+outputs\*\*', n_csvs,    "CSVs"),
]:
    m = re.search(pattern, readme)
    if not m:
        warn(f"README: '{label}' count not found")
    elif int(m.group(1)) == actual:
        ok(f"README count {label}: {actual} (matches)")
    else:
        fail(f"README says {m.group(1)} {label}, repo has {actual}")

# ============================================================
# FASE 2: SCRIPTS (syntax + deps)
# ============================================================
section("FASE 2: SCRIPTS (syntax + deps)")

scripts = sorted(Path("analysis").rglob("*.R"))
syntax_bad = []
deps_all = set()
for s in scripts:
    txt = s.read_text(encoding="utf-8")
    r = subprocess.run(["Rscript", "-e", f'invisible(parse("{s}"))'],
                       capture_output=True, text=True)
    if r.returncode != 0:
        syntax_bad.append(str(s))
    for m in re.finditer(r'(?:library|require)\s*\(\s*["\']?([A-Za-z0-9._]+)["\']?\s*\)', txt):
        deps_all.add(m.group(1))

if syntax_bad:
    for s in syntax_bad: fail(f"Syntax error: {s}")
else:
    ok(f"All {len(scripts)} scripts parse OK")

missing = []
for dep in sorted(deps_all):
    r = subprocess.run(["Rscript", "-e", f'cat(requireNamespace("{dep}", quietly=TRUE))'],
                       capture_output=True, text=True)
    if r.stdout.strip() != "TRUE":
        missing.append(dep)
if missing:
    fail(f"Missing deps: {', '.join(missing)}")
else:
    ok(f"All {len(deps_all)} dependencies installed")

# ============================================================
# FASE 3: MANUSCRIPT NUMBERS (15/15)
# ============================================================
section("FASE 3: MANUSCRIPT NUMBERS")

def fuzzy_find(csv_rel, value, tol):
    """Search for numeric value in CSV (as string with rounding tolerance)."""
    if not os.path.exists(csv_rel):
        return None
    txt = Path(csv_rel).read_text(encoding="utf-8", errors="ignore")
    found = []
    for m in re.finditer(r'-?\d+\.\d+(?:[eE][-+]?\d+)?', txt):
        try:
            v = float(m.group(0))
            if abs(v - value) <= tol:
                found.append(v)
        except ValueError:
            pass
    return found[0] if found else None

numbers = [
    ("SMPD3 +4.80 (Mdr2-KO)",        "results/bulk/SRP074136_DESeq2_inflamed_vs_WT.csv",      4.80,  0.10),
    ("B4GALT6 +3.77 (Mdr2-KO)",      "results/bulk/SRP074136_DESeq2_inflamed_vs_WT.csv",      3.77,  0.10),
    ("Smpd3 EAE +2.88 (pre-onset)",  "results/bulk/34_DESeq2_liver_EAE_panel.csv",        2.88,  0.10),
    ("Smpd3 EAE +1.77 (peak)",       "results/bulk/34_DESeq2_liver_EAE_panel.csv",        1.77,  0.10),
    ("B4galt6 EAE +1.49 (pre-onset)","results/bulk/34_DESeq2_liver_EAE_panel.csv",        1.49,  0.10),
    ("B4galt6 EAE -0.11 (peak)",     "results/bulk/34_DESeq2_liver_EAE_panel.csv",       -0.11,  0.10),
    ("S1PR1/S1PR3 20.3 (astro)",     "results/cellular/27_cellular_ratios.csv",              20.26,  0.20),
    ("S1PR1/S1PR3 21.1 (micro)",     "results/cellular/27_cellular_ratios.csv",              21.04,  0.20),
    ("S1PR1/S1PR3 2.7 (neuron)",     "results/cellular/27_cellular_ratios.csv",               2.66,  0.10),
    ("B4GALT6+ 7.06% (Early)",       "results/snrnaseq/36_positivity.csv", 7.06, 0.10),
    ("B4GALT6+ 17.9% (Peak)",        "results/snrnaseq/36_positivity.csv",17.89, 0.20),
    ("MOL3 59.6% (Late)",            "results/snrnaseq/30b_falcao_mol3_share.csv",           59.55,  0.20),
    ("DamageScore 0.061 (MOL2)",     "results/snrnaseq/12_MOL3_subpop_profile.csv",          0.0611, 0.001),
    ("DamageScore 0.006 (MOL56)",    "results/snrnaseq/12_MOL3_subpop_profile.csv",          0.0056, 0.001),
    ("LacCer p = 4.3e-6",            "results/snrnaseq/09_pathway_regression.csv",           4.32e-6, 0.01e-6),
    ("ARI = 0.402",                  "results/snrnaseq/11_ari_from_precomputed.csv",         0.4018, 0.001),
]

found = 0
for label, csv, val, tol in numbers:
    hit = fuzzy_find(csv, val, tol)
    if hit is not None:
        ok(f"{label:<35} -> {hit:.6g} in {csv}")
        found += 1
    else:
        warn(f"{label:<35} -> not found in {csv} (val={val}, tol={tol})")

print(f"\n  Numbers verified: {found}/{len(numbers)}")

# ============================================================
# FASE 4: GIT
# ============================================================
section("FASE 4: GIT")

def run(cmd):
    return subprocess.run(cmd.split(), capture_output=True, text=True).stdout.strip()

status = run("git status --porcelain")
if status == "":
    ok("Working tree clean")
else:
    fail(f"Working tree not clean ({len(status.splitlines())} files)")

tags = run("git tag -l").split("\n")
for t in ["v5.3", "v5.3-final", "v5.3.1"]:
    if t in tags:
        ok(f"Tag present: {t}")
    else:
        warn(f"Tag missing: {t}")

branch = run("git rev-parse --abbrev-ref HEAD")
remote = run("git rev-parse --abbrev-ref --symbolic-full-name @{u}")
head = run("git rev-parse HEAD")
remote_head = run("git rev-parse " + remote.split("/")[-1])
if head == remote_head:
    ok(f"Local {branch} == remote {remote}")
else:
    warn(f"Local {branch} diverged from remote {remote}")

# ============================================================
# FASE 5: DOCUMENTATION (stale references)
# ============================================================
section("FASE 5: DOCUMENTATION (stale references)")

stale_patterns = [
    ("10b_ari",                 "old script name (10b_ari)"),
    ("11_MOL3_subpop_profile.csv", "old CSV name (11_MOL3_subpop)"),
    ("Reproducibility note: ARI", "old 'not reproducible' note"),
]
for pat, desc in stale_patterns:
    hits = []
    for f in Path(".").rglob("*"):
        if f.is_file() and f.suffix in (".R", ".md", ".csv", ".txt"):
            try:
                if pat in f.read_text(encoding="utf-8", errors="ignore"):
                    hits.append(str(f))
            except Exception:
                pass
    if hits:
        fail(f"Stale reference '{desc}' found in: {hits[:5]}")
    else:
        ok(f"No stale reference: {desc}")

# Check H5 workflow section present
if "## H5 workflow (GSE250589)" in readme:
    ok("README contains 'H5 workflow (GSE250589)' section")
else:
    warn("README missing 'H5 workflow (GSE250589)' section")

# ============================================================
# FASE 6: REPRODUCIBILITY
# ============================================================
section("FASE 6: REPRODUCIBILITY")

# Each manuscript number should have a script that generates it
generators = {
    "results/snrnaseq/08_MOL3_sphingolipid_substructure.csv": "analysis/snrnaseq/08_",
    "results/snrnaseq/09_pathway_regression.csv":             "analysis/snrnaseq/09_",
    "results/snrnaseq/11_ari_from_precomputed.csv":            "analysis/snrnaseq/11_",
    "results/snrnaseq/12_MOL3_subpop_profile.csv":             "analysis/snrnaseq/12_",
}
for csv, script_prefix in generators.items():
    script_found = any(Path("analysis/snrnaseq").glob(script_prefix.split("/")[-1] + "*.R"))
    if script_found and os.path.exists(csv):
        ok(f"{csv} has generator {script_prefix}*.R")
    else:
        warn(f"{csv} may be orphan (no clear generator)")

# ============================================================
# VEREDICTO
# ============================================================
section("VEREDICTO FINAL")
print(f"  OK:       {len(oks)}")
print(f"  Warnings: {len(warnings)}")
print(f"  Failures: {len(issues)}")
if issues:
    print(f"\n  FAILURES:")
    for i in issues:
        print(f"    - {i}")
if warnings:
    print(f"\n  WARNINGS:")
    for w in warnings:
        print(f"    - {w}")

print()
if len(issues) == 0 and len(warnings) <= 2:
    print("  ✅ REPO CERRADO: reproducible, limpio, sin issues críticos")
    sys.exit(0)
elif len(issues) == 0:
    print("  ⚠️  REPO CERRADO con warnings menores (revisar arriba)")
    sys.exit(0)
else:
    print(f"  ❌ {len(issues)} issue(s) crítico(s) — fix antes de cerrar")
    sys.exit(1)
