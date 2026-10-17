#!/usr/bin/env python3
"""Energy-generating cycles (EGCs) of every model, tested with the bounds the model actually has.

An EGC lets a model charge an energy carrier (ATP, NADH, a quinol, ...) without taking anything up. The test follows
Fritzemeier et al. 2017 (PLoS Comput Biol 13:e1005494), as in MEMOTE: every boundary reaction (exchange, sink, demand) is
closed, a dissipation reaction for the carrier is added (e.g. atp + h2o -> adp + pi + h) and its flux is maximised; any
flux > 0 means an EGC exists. Unlike MEMOTE's version, the bounds of the internal reactions are kept as they are in the
model file (MEMOTE resets them before the test, which re-opens reactions that were closed during curation).

Input   the models (*.xml), one sub-folder per set of models, for example
        models/carveme_m9/bin.1.xml and models/curated/bin.1.xml.
Output  energy_cycles.tsv: one row per model, with its set, the maximal dissipation flux of every energy couple
        (NaN if the model does not contain the couple) and the number of couples with an EGC.

Environment  Python 3.10, cobra 0.31.1, pandas 2.3.3, and a solver (CPLEX 22.2 or the free GLPK through
        swiglpk 5.0.13).

Usage   python energy_cycles.py [--solver cplex|glpk]
"""
import argparse
import glob
import os
import sys
import warnings

import cobra
import pandas as pd
from cobra import Reaction

warnings.filterwarnings("ignore")

# --- Parameters ----------------------------------------------------------------------------------
SETS = ["carveme_m9", "curated"]          # one sub-folder per set of models
# energy couples of MEMOTE (support/consistency.py), as BiGG metabolite IDs: charged form, discharged form, dissipation type
COUPLES = {
    "ATP": ("atp", "adp", "ntp"), "CTP": ("ctp", "cdp", "ntp"), "GTP": ("gtp", "gdp", "ntp"),
    "UTP": ("utp", "udp", "ntp"), "ITP": ("itp", "idp", "ntp"),
    "NADH": ("nadh", "nad", "nad"), "NADPH": ("nadph", "nadp", "nad"),
    "FADH2": ("fadh2", "fad", "redox"), "FMNH2": ("fmnh2", "fmn", "redox"), "Q8H2": ("q8h2", "q8", "redox"),
    "MQL8": ("mql8", "mqn8", "redox"), "MQL7": ("mql7", "mqn7", "redox"), "MQL6": ("mql6", "mqn6", "redox"),
    "DMMQL8": ("2dmmql8", "2dmmq8", "redox"),
    "ACCOA": ("accoa", "coa", "accoa"), "GLU": ("glu__L", "akg", "glu"),
}
TOL = 1e-6


def met(model, bigg):
    mid = f"{bigg}_c"
    return model.metabolites.get_by_id(mid) if mid in model.metabolites else None


def dissipation(model, charged, discharged, kind):
    a, b = met(model, charged), met(model, discharged)
    if a is None or b is None:
        return None
    h, h2o = met(model, "h"), met(model, "h2o")
    r = Reaction("dissipation")
    if kind == "ntp":          # atp + h2o -> adp + pi + h
        r.add_metabolites({a: -1, h2o: -1, b: 1, met(model, "pi"): 1, h: 1})
    elif kind == "nad":        # nadh -> nad + h
        r.add_metabolites({a: -1, b: 1, h: 1})
    elif kind == "redox":      # q8h2 -> q8 + 2 h
        r.add_metabolites({a: -1, b: 1, h: 2})
    elif kind == "accoa":      # accoa + h2o -> coa + ac + h
        r.add_metabolites({a: -1, h2o: -1, b: 1, met(model, "ac"): 1, h: 1})
    elif kind == "glu":        # glu__L + h2o -> akg + nh4 + 2 h
        r.add_metabolites({a: -1, h2o: -1, b: 1, met(model, "nh4"): 1, h: 2})
    r.bounds = (0, 1000)
    return r


p = argparse.ArgumentParser()
p.add_argument("--solver", default="cplex")
p.add_argument("--models", default="models")
p.add_argument("--out", default="results/memote/energy_cycles.tsv")
a = p.parse_args()
cobra.Configuration().solver = a.solver

rows = []
for s in SETS:
    files = sorted(glob.glob(f"{a.models}/{s}/*.xml"))
    if not files:
        sys.exit(f"STOP: no models in {a.models}/{s}")
    for f in files:
        model = cobra.io.read_sbml_model(f)
        for r in model.boundary:                       # nothing in, nothing out
            r.bounds = (0, 0)
        if model.slim_optimize(error_value=0.0) > TOL:  # self-check: with closed boundaries there can be no growth
            sys.exit(f"STOP: {f} grows with all boundary reactions closed")
        row = {"model_set": s, "model": os.path.basename(f)[:-4]}
        n = 0
        for name, (charged, discharged, kind) in COUPLES.items():
            r = dissipation(model, charged, discharged, kind)
            if r is None:
                row[name] = float("nan")
                continue
            with model:
                model.add_reactions([r])
                model.objective = r
                v = model.slim_optimize(error_value=float("nan"))
            row[name] = v
            n += v > TOL
        row["couples_with_cycle"] = n
        rows.append(row)
        print(f"{s}/{row['model']}: {n} energy couples with a cycle", flush=True)

os.makedirs(os.path.dirname(a.out), exist_ok=True)
df = pd.DataFrame(rows)
df.to_csv(a.out, sep="\t", index=False, float_format="%.6g")
print("\nmodels (of 90) in which the energy couple can be charged with every boundary reaction closed:")
print((df.set_index(["model_set", "model"])[list(COUPLES)] > TOL).groupby("model_set").sum().T.to_string())
