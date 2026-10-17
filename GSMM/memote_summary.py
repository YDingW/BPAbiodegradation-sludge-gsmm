#!/usr/bin/env python3
"""Summary of the MEMOTE reports of the two model sets (CarveMe output and curated models).

Input   the MEMOTE report of every model as JSON, one folder per set (<OUT>/carveme_m9/*.json and
        <OUT>/curated/*.json, written by `memote run --solver-timeout 0 --filename <model>.json <model>.xml`),
        and, if present, <OUT>/energy_cycles.tsv (written by energy_cycles.py).
Output  memote_metrics.tsv       one row per model and set with the 20 MEMOTE metrics below (+ energy cycles)
        memote_summary.tsv       median, mean, minimum and maximum of every metric per set
        memote_distribution.svg / .png   distribution of every metric over the 90 models of each set
Fractions are those reported by MEMOTE (e.g. the fraction of reactions that are blocked); for mass balance, charge balance
and metabolite connectivity the fraction of reactions/metabolites that pass is shown. Stoichiometric consistency is
pass/fail (the number of models that pass is printed and is the mean in memote_summary.tsv); it and the two metrics that
are 0 in every model (metabolites without formula / charge) are in the tables but not in the figure.

Environment  Python 3.10, memote 0.17.0 (to produce the reports), pandas 2.3.3, numpy 2.2.6, matplotlib 3.10.8.

Usage   python memote_summary.py
"""
import glob
import json
import os
import sys

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

# --- Parameters ----------------------------------------------------------------------------------
SETS = {"carveme_m9": "CarveMe", "curated": "curated"}   # folder of each model set: display name
OUT = "results/memote"                                   # where the reports are and the tables are written
# display name: (memote test, how the value is taken from the test: the metric, 1 - metric, pass/fail, or the number of
# items in data; axis label)
METRICS = {
    "Total reactions": ("test_reactions_presence", "count", "n"),
    "Total metabolites": ("test_metabolites_presence", "count", "n"),
    "Total genes": ("test_genes_presence", "count", "n"),
    "Medium components": ("test_find_medium_metabolites", "count", "n (default medium)"),
    "Purely metabolic reactions": ("test_find_pure_metabolic_reactions", "metric", "fraction of reactions"),
    "Transport reactions": ("test_find_transport_reactions", "metric", "fraction of reactions"),
    "Metabolic coverage": ("test_metabolic_coverage", "metric", "reactions per gene"),
    "Stoichiometric consistency": ("test_stoichiometric_consistency", "passed", "models passing"),
    "Mass balance": ("test_reaction_mass_balance", "pass_fraction", "fraction of reactions balanced"),
    "Charge balance": ("test_reaction_charge_balance", "pass_fraction", "fraction of reactions balanced"),
    "Metabolites without formula": ("test_metabolites_formula_presence", "metric", "fraction of metabolites"),
    "Metabolites without charge": ("test_metabolites_charge_presence", "metric", "fraction of metabolites"),
    "Metabolite connectivity": ("test_find_disconnected", "pass_fraction", "fraction of metabolites connected"),
    "Universally blocked reactions": ("test_blocked_reactions", "metric", "fraction of reactions"),
    "Orphan metabolites": ("test_find_orphans", "metric", "fraction of metabolites"),
    "Dead-end metabolites": ("test_find_deadends", "metric", "fraction of metabolites"),
    "Not consumed in complete medium": ("test_find_metabolites_not_consumed_with_open_bounds", "metric", "fraction of metabolites"),
    "Not produced in complete medium": ("test_find_metabolites_not_produced_with_open_bounds", "metric", "fraction of metabolites"),
    "Unbounded flux in default medium": ("test_find_reactions_unbounded_flux_default_condition", "metric", "fraction of reactions"),
    "Stoichiometrically balanced cycles": ("test_find_stoichiometrically_balanced_cycles", "metric", "fraction of reactions"),
}
EGC = "Energy-generating cycles"          # from energy_cycles.tsv, not from MEMOTE (see scripts/energy_cycles.py)
NOT_PLOTTED = {"Stoichiometric consistency", "Metabolites without formula", "Metabolites without charge"}  # pass/fail; 0 in every model


def value(test, kind):
    if test is None:
        return float("nan")
    if kind == "count":
        return len(test["data"]) if isinstance(test.get("data"), list) else float("nan")
    if kind == "passed":
        return 1.0 if test.get("result") == "passed" else 0.0
    m = test.get("metric")
    if not isinstance(m, (int, float)):
        return float("nan")
    return 1.0 - m if kind == "pass_fraction" else float(m)


rows = []
for folder, label in SETS.items():
    files = sorted(glob.glob(f"{OUT}/{folder}/*.json"))
    if not files:
        sys.exit(f"STOP: no reports in {OUT}/{folder}")
    for f in files:
        tests = json.load(open(f))["tests"]
        row = {"model_set": label, "model": os.path.basename(f)[:-5]}
        row.update({name: value(tests.get(test), kind) for name, (test, kind, _) in METRICS.items()})
        rows.append(row)
metrics = pd.DataFrame(rows)
columns = list(METRICS)
egc_file = f"{OUT}/energy_cycles.tsv"
if os.path.exists(egc_file):
    egc = pd.read_csv(egc_file, sep="\t")
    egc["model_set"] = egc["model_set"].map(SETS)
    metrics = metrics.merge(egc[["model_set", "model", "couples_with_cycle"]].rename(columns={"couples_with_cycle": EGC}),
                            on=["model_set", "model"], how="left")
    columns.append(EGC)
else:
    print(f"note: {egc_file} not found, the energy-cycle panel is left out")
metrics.to_csv(f"{OUT}/memote_metrics.tsv", sep="\t", index=False, float_format="%.6g")
summary = metrics.groupby("model_set", sort=False)[columns].agg(["median", "mean", "min", "max"])
summary.columns = [f"{m} ({s})" for m, s in summary.columns]
summary.T.to_csv(f"{OUT}/memote_summary.tsv", sep="\t", float_format="%.6g")
print(f"{len(metrics)} reports; medians:\n" + summary.T.filter(like="median", axis=0).to_string())
passing = metrics.groupby("model_set", sort=False)["Stoichiometric consistency"].sum().astype(int)
print("stoichiometric consistency, models passing: " + ", ".join(f"{k} {v}/90" for k, v in passing.items()))

# memote_distribution.svg / .png: one panel per metric, one violin per model set with the individual models as points and the median as a bar
colours = {"CarveMe": "#eb6834", "curated": "#2a78d6"}
panels = [c for c in columns if c not in NOT_PLOTTED]
fig, axes = plt.subplots(3, 6, figsize=(17, 9))
rng = np.random.default_rng(0)
for ax, name in zip(axes.flat, panels):
    labels = list(SETS.values())
    for i, label in enumerate(labels):
        v = metrics.loc[metrics["model_set"] == label, name].dropna().values
        if len(v) == 0:
            continue
        if np.ptp(v) > 0:
            parts = ax.violinplot(v, positions=[i], widths=0.8, showextrema=False)
            for body in parts["bodies"]:
                body.set_facecolor(colours[label]); body.set_edgecolor("none"); body.set_alpha(0.45)
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, s=5, color="black", alpha=0.35, linewidths=0)
        ax.plot([i - 0.2, i + 0.2], [np.median(v)] * 2, color="black", linewidth=1.5)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_title(name, fontsize=9)
    ax.set_ylabel(METRICS[name][2] if name in METRICS else "energy couples with a cycle (of 16)", fontsize=7)
    ax.tick_params(axis="y", labelsize=7)
    ax.spines[["top", "right"]].set_visible(False)
for ax in list(axes.flat)[len(panels):]:                  # unused panels
    ax.axis("off")
fig.suptitle("MEMOTE metrics of the 90 models before (CarveMe) and after curation; points = models, bar = median", fontsize=10)
fig.tight_layout()
fig.savefig(f"{OUT}/memote_distribution.svg")
fig.savefig(f"{OUT}/memote_distribution.png", dpi=300)
