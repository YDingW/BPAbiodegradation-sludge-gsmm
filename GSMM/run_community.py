#!/usr/bin/env python3
"""Community simulation of one time point: the members and their abundances, growth rates pinned to the
values derived from 16S qPCR and MAG abundances, and the fluxes chosen by expression-weighted pFBA.

1. Members: the MAGs marked as members in results/eflux/<time point>_mags.tsv, with their metaWRAP abundance
   of that time point (normalised to 1 over the members).
2. Models: the curated model of each member, with every constrained reaction given its E-Flux bounds from
   results/eflux/<time point>_reactions.tsv. The models are merged into a MICOM community.
3. Medium: the column of the time point in MEDIUM_FILE (components with a maximal uptake of 0 are left out).
   The BPA uptake of the community is confined to the fraction range of its maximal uptake given for the time
   point.
4. Growth: (i) the maximal community growth rate is recorded for reference; (ii) every member is required to
   grow at its target rate mu_i (results/growth_rates.tsv): maximise t with mu_i >= t * target_i for every
   member with a target; if t < 1 the members whose constraint is binding cannot reach their target, they
   are pinned at the rate they reached, and the LP is repeated for the others; (iii) with all members held at
   their target (or pinned) rates, minimise sum_j price_j |v_j| over all reactions, where price = (1 / c)^PRICE_POWER
   for a reaction with expression coefficient c and PRICE_OTHER for every other reaction.
5. Checks: no member below its pinned rate, carbon entering the community >= carbon leaving it, and no flux
   through boundary reactions outside the shared medium.

Input   the curated models (*.xml), the MAG abundance table, the growth rates (growth_rates.py), the members
        and E-Flux bounds of the time point (eflux_bounds.py) and medium.tsv (one row per exchange reaction,
        one column per time point: the maximal uptake of the community in mmol/gDW/h).
Output  growth.tsv (member, abundance, target, pinned, achieved), fluxes.tsv (member x reaction, mmol/gDW/h of
        the member; the "medium" row holds the community exchange fluxes), medium.tsv (what was applied) and
        summary.tsv.

Environment  Python 3.10, micom 0.39.0, cobra 0.31.1, optlang 1.9.0, pandas 2.3.3, and IBM ILOG CPLEX 22.2 with
        its Python interface (the full version; the Community Edition is too small for these models).

Usage   python run_community.py --time-point 24h
        python run_community.py --time-point 24h --out community/24h
"""
import argparse
import csv
import logging
import os
import shutil
import sys
import tempfile
import time

import cobra
import pandas as pd
from micom import Community
from micom.solution import CommunitySolution
from micom.util import _apply_min_growth
from optlang.symbolics import Zero

# --- Parameters --------------------------------------------------------------------------------
# Per time point: the column of the abundance table, the range of the community BPA uptake as a fraction of
# its maximum, and the MAGs that are members although their transcript detection is below the cutoff of
# eflux_bounds.py.
TIME_POINTS = {
    "0h":  {"abundance_column": "0H",  "bpa_uptake_fraction": [0.01, 1.0], "forced_members": ["bin.47", "bin.104"]},
    "24h": {"abundance_column": "24H", "bpa_uptake_fraction": [0.7, 1.0],  "forced_members": []},
    "26h": {"abundance_column": "26H", "bpa_uptake_fraction": [0.7, 1.0],  "forced_members": []},
    "30h": {"abundance_column": "30H", "bpa_uptake_fraction": [0.7, 1.0],  "forced_members": []},
}
MEDIUM_FILE = "medium.tsv"             # shared medium: one row per exchange reaction, one column per time point
BPA_EXCHANGE = "EX_bpa_m"              # BPA exchange of the community
GLUCOSE_EXCHANGE = "EX_glc__D_m"       # glucose exchange of the community (reported in summary.tsv)
PRICE_POWER = 2                        # price of a reaction with expression coefficient c is (1 / c) ** 2
PRICE_OTHER = 1600                     # price of a reaction without an expression coefficient
SOLVER = "cplex"
THREADS = 4
TIME_LIMIT_S = 21600                   # per LP
GROWTH_RTOL = GROWTH_ATOL = 1e-6       # tolerance when checking that a member reached its target rate
CARBON_CLOSURE_MIN = 0.9               # minimal ratio of the carbon leaving the community to the carbon entering it
# -----------------------------------------------------------------------------------------------

CPLEX_BASE = {"advance": 0, "numerical": 1, "scale": 1, "barrier_convergetol": 1e-9, "feasibility": 1e-6, "optimality": 1e-6}
CPLEX_LP = {"dual_simplex": dict(CPLEX_BASE, lpmethod=2, presolve=True),          # used for the growth and pFBA problems
            "barrier": dict(CPLEX_BASE, lpmethod=4, scale=0, solutiontype=1, presolve=False)}   # reference maximum only


def log(msg):
    print(time.strftime("%H:%M:%S ") + msg, flush=True)


def read_tsv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def configure_solver(com, method, threads, time_limit):
    """CPLEX parameters: numerical emphasis, cold start, chosen LP method; the simplex runs with presolve."""
    if "cplex" not in type(com.solver).__module__:
        return                                          # other solvers (tests) keep their defaults
    p = com.solver.problem.parameters
    cfg = CPLEX_LP[method]
    p.advance.set(cfg["advance"]); p.emphasis.numerical.set(cfg["numerical"]); p.read.scale.set(cfg["scale"])
    p.lpmethod.set(cfg["lpmethod"]); p.qpmethod.set(cfg["lpmethod"]); p.solutiontype.set(cfg.get("solutiontype", 0))
    p.barrier.convergetol.set(cfg["barrier_convergetol"]); p.barrier.algorithm.set(0); p.barrier.startalg.set(1)
    p.simplex.tolerances.feasibility.set(cfg["feasibility"]); p.simplex.tolerances.optimality.set(cfg["optimality"])
    p.threads.set(threads); p.timelimit.set(time_limit)
    com.solver.configuration.presolve = cfg["presolve"]


def solve(com, label):
    t0 = time.time()
    com.solver.optimize()
    status = str(com.solver.status)
    value = float(com.solver.objective.value) if status == "optimal" else None
    log(f"   {label}: {status}" + (f", objective {value:.6g}" if value is not None else "") + f" ({(time.time() - t0) / 60:.1f} min)")
    if status != "optimal":
        sys.exit(f"STOP: {label} did not reach an optimal solution (status {status})")
    return value


def close_free_sources(com):
    """MICOM leaves the boundary reactions of metabolites whose id contains 'bm' outside the shared medium (it takes
    them for biomass reactions). Such a reaction must not act as a free source: allow it only to export."""
    exchange_ids = {r.id for r in com.exchanges}
    free = [r for r in com.reactions if r.boundary and r.id not in exchange_ids]
    closed = []
    for r in free:
        for met, coef in r.metabolites.items():
            if coef < 0 and r.lower_bound < 0:
                r.lower_bound = 0.0; closed.append(r.id)
            elif coef > 0 and r.upper_bound > 0:
                r.upper_bound = 0.0; closed.append(r.id)
    return free, closed


def produced_by_free_sources(free, fluxes):
    """Largest net production (stoichiometry x flux) through the boundary reactions outside the medium."""
    worst = 0.0
    for r in free:
        v = fluxes.get(r.id)
        if v is None:
            continue
        for met, coef in r.metabolites.items():
            worst = max(worst, coef * v)
    return worst


def carbon_closure(com, fluxes):
    """Carbon entering / carbon leaving the community through the medium (>= 1 when nothing is created)."""
    c_in = c_out = 0.0
    for r in com.exchanges:
        carbons = max((m.elements.get("C", 0) if m.formula else 0) for m in r.metabolites)
        v = fluxes.get(r.id, 0.0)
        if carbons <= 0 or v != v:
            continue
        if v < 0:
            c_in += -v * carbons
        elif v > 0:
            c_out += v * carbons
    return (float("inf") if c_out <= 1e-12 else c_in / c_out), c_in, c_out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--time-point", required=True)
    ap.add_argument("--models", default="models/curated")
    ap.add_argument("--abundance", default="data/mag_abundance.csv")
    ap.add_argument("--growth-rates", default="results/growth_rates.tsv")
    ap.add_argument("--eflux", default="results/eflux", help="folder with <time point>_mags.tsv and <time point>_reactions.tsv")
    ap.add_argument("--out", default=None, help="default results/community/<time point>")
    ap.add_argument("--solver", default=None, help="override SOLVER (tests only)")
    a = ap.parse_args()
    tp = a.time_point
    if tp not in TIME_POINTS:
        sys.exit(f"STOP: time point {tp} is not in TIME_POINTS")
    tcfg = TIME_POINTS[tp]
    solver = a.solver or SOLVER
    out = a.out or f"results/community/{tp}"
    os.makedirs(out, exist_ok=True)
    logging.getLogger("cobra").setLevel(logging.ERROR)
    logging.getLogger("micom").setLevel(logging.ERROR)
    t_start = time.time()
    log(f"time point {tp}, solver {solver}, output {out}")

    # 1 members and abundances
    mags = [r for r in read_tsv(f"{a.eflux}/{tp}_mags.tsv") if r["member"] != "no"]
    members = [r["mag"] for r in mags]
    forced = [r["mag"] for r in mags if r["member"] == "forced"]
    if sorted(forced) != sorted(tcfg["forced_members"]):
        sys.exit(f"STOP: forced members in {a.eflux}/{tp}_mags.tsv are {forced}, TIME_POINTS says {tcfg['forced_members']}")
    abundance = {r[list(r)[0]]: float(r[tcfg["abundance_column"]]) for r in csv.DictReader(open(a.abundance, newline=""))}
    missing = [m for m in members if m not in abundance]
    if missing:
        sys.exit(f"STOP: no abundance for {missing[:5]}")
    total = sum(abundance[m] for m in members)
    share = {m: abundance[m] / total for m in members}
    targets = {r["mag"]: float(r["growth_rate"]) for r in read_tsv(a.growth_rates) if r["time_point"] == tp}
    missing = [m for m in members if m not in targets]
    if missing:
        sys.exit(f"STOP: no growth rate for {missing[:5]} at {tp}")
    log(f"1  {len(members)} members ({len(forced)} forced); growth rates for all of them")

    # 2 models with their E-Flux bounds, merged into the community
    bounds = {}
    for r in read_tsv(f"{a.eflux}/{tp}_reactions.tsv"):
        bounds.setdefault(r["mag"], {})[r["reaction"]] = (float(r["eflux_lower_bound"]), float(r["eflux_upper_bound"]))
    tmp = tempfile.mkdtemp(prefix="community_")
    rows = []
    for m in members:
        model = cobra.io.read_sbml_model(f"{a.models}/{m}.xml")
        for rid, (lb, ub) in bounds.get(m, {}).items():
            model.reactions.get_by_id(rid).bounds = (lb, ub)
        cobra.io.write_sbml_model(model, f"{tmp}/{m}.xml")
        rows.append({"id": m.replace(".", "_"), "genus": m, "species": m, "abundance": share[m],
                     "file": f"{tmp}/{m}.xml", "sample_id": tp})
    taxonomy = pd.DataFrame(rows)
    log(f"2  models prepared ({sum(len(b) for b in bounds.values())} E-Flux bounds applied); merging")
    com = Community(taxonomy, solver=solver, progress=False)
    shutil.rmtree(tmp)
    free, closed = close_free_sources(com)
    threads, time_limit = THREADS, TIME_LIMIT_S
    configure_solver(com, "dual_simplex", threads, time_limit)
    log(f"   community: {len(com.reactions)} reactions, {len(com.metabolites)} metabolites; "
        f"{len(free)} boundary reactions outside the medium, {len(closed)} closed as sources ({(time.time() - t_start) / 60:.1f} min)")

    # 3 medium
    medium = {r["exchange"]: float(r[tp]) for r in read_tsv(MEDIUM_FILE) if float(r[tp]) > 0}
    exchange_ids = {r.id for r in com.exchanges}
    missing = [k for k in medium if k not in exchange_ids]
    if missing:
        sys.exit(f"STOP: medium components without an exchange reaction in the community: {missing}")
    com.medium = medium
    bpa = com.reactions.get_by_id(BPA_EXCHANGE)
    lo, hi = (f * medium[bpa.id] for f in tcfg["bpa_uptake_fraction"])
    bpa.bounds = (-hi, -lo)
    with open(f"{out}/medium.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["exchange", "max_uptake"])
        w.writerows(sorted(medium.items()))
    log(f"3  medium: {len(medium)} components; BPA uptake between {lo:.6g} and {hi:.6g}")

    # 4 growth
    taxa = list(com.taxa)
    target = {sp: targets[sp.replace("_", ".")] for sp in taxa}
    need = [sp for sp in taxa if target[sp] > 1e-12]
    rtol, atol = GROWTH_RTOL, GROWTH_ATOL
    required = sum(float(com.abundances[sp]) * target[sp] for sp in taxa)
    with com:
        _apply_min_growth(com, {sp: 0.0 for sp in taxa}, atol, rtol)
        configure_solver(com, "barrier", threads, time_limit)
        com.objective = com.scale * com.variables.community_objective
        mu_max = solve(com, "4i maximal community growth") / com.scale
        configure_solver(com, "dual_simplex", threads, time_limit)
    log(f"   required community growth {required:.6f} = {required / mu_max:.3f} of the maximum {mu_max:.6f}")
    pinned, rounds = {}, 0
    while True:
        rounds += 1
        open_members = [sp for sp in need if sp not in pinned]
        if not open_members:
            break
        with com:
            t = com.problem.Variable("t", lb=0.0, ub=1.0)
            cons = {sp: com.problem.Constraint(com.constraints["objective_" + sp].expression - target[sp] * t, lb=0.0, name=f"t_{sp}")
                    for sp in open_members}
            com.add_cons_vars([t] + list(cons.values()))
            if pinned:
                _apply_min_growth(com, {sp: pinned.get(sp, 0.0) for sp in taxa}, atol, rtol)
            com.objective = com.problem.Objective(t, direction="max")
            t_star = solve(com, f"4ii round {rounds}, max t for {len(open_members)} members")
            if t_star >= 1.0 - 1e-6:
                break
            duals = {sp: abs(float(c.dual)) for sp, c in cons.items()}
            primals = {sp: abs(float(c.primal)) for sp, c in cons.items()}
            reached = {sp: float(com.constraints["objective_" + sp].primal) for sp in open_members}
        tight = [sp for sp in open_members if duals[sp] > 1e-9 or primals[sp] <= 1e-9 * max(1.0, target[sp])]
        pick = [sp for sp in tight if duals[sp] > 1e-9] or tight[:1]
        if not pick:
            sys.exit(f"STOP: t = {t_star:.6f} < 1 but no binding growth constraint was found")
        for sp in pick:
            pinned[sp] = max(0.0, reached[sp])
            log(f"   {sp} cannot reach {target[sp]:.6f}, pinned at {pinned[sp]:.6f} ({pinned[sp] / target[sp] * 100:.1f} %)")
    final = {sp: pinned.get(sp, target[sp]) for sp in taxa}
    prices, coef = {}, {}
    for r in read_tsv(f"{a.eflux}/{tp}_reactions.tsv"):
        coef[f"{r['reaction']}__{r['mag'].replace('.', '_')}"] = float(r["coefficient"])
    k, other = float(PRICE_POWER), float(PRICE_OTHER)
    for r in com.reactions:
        price = (1.0 / coef[r.id]) ** k if r.id in coef else other
        prices[r.forward_variable] = price
        prices[r.reverse_variable] = price
    missing = [rid for rid in coef if rid not in com.reactions]
    if missing:
        sys.exit(f"STOP: {len(missing)} reactions of the coefficient table are not in the community, e.g. {missing[:3]}")
    _apply_min_growth(com, final, atol, rtol)
    com.objective = Zero
    com.objective_direction = "min"
    com.objective.set_linear_coefficients(prices)
    total_cost = solve(com, f"4iii weighted pFBA ({len(coef)} priced reactions, others {other:g})")
    sol = CommunitySolution(com)

    # 5 checks and outputs
    rates = {sp: float(sol.members.at[sp, "growth_rate"]) for sp in taxa}
    low = [sp for sp in need if rates[sp] < (1 - rtol) * final[sp] - atol - 1e-9]
    if low:
        sys.exit(f"STOP: {len(low)} members ended below their pinned growth rate, e.g. {low[:3]}")
    fluxes = {r.id: float(r.flux) for r in com.reactions}
    worst = produced_by_free_sources(free, fluxes)
    if worst > 1e-6:
        sys.exit(f"STOP: a boundary reaction outside the medium produces metabolite ({worst:.3g})")
    closure, c_in, c_out = carbon_closure(com, fluxes)
    if closure < CARBON_CLOSURE_MIN:
        sys.exit(f"STOP: carbon closure {closure:.3f}: more carbon leaves the community than enters it")
    growth = float(sol.growth_rate)
    with open(f"{out}/growth.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["member", "abundance", "target_growth_rate", "pinned_growth_rate", "growth_rate", "pinned"])
        for sp in taxa:
            w.writerow([sp.replace("_", "."), repr(float(com.abundances[sp])), repr(target[sp]), repr(final[sp]), repr(rates[sp]),
                        "yes" if sp in pinned else "no"])
    fl = sol.fluxes.copy()
    fl.index = [i.replace("_", ".") if i != "medium" else i for i in fl.index]
    fl.index.name = "member"
    fl.to_csv(f"{out}/fluxes.tsv", sep="\t", float_format="%.12g")
    uptake = {m: -fluxes[m] for m in medium}
    summary = {"time_point": tp, "members": len(taxa), "forced_members": len(forced),
               "max_community_growth": repr(mu_max), "required_community_growth": repr(required), "community_growth": repr(growth),
               "achieved_fraction": repr(growth / required if required > 0 else float("nan")), "pinning_rounds": rounds,
               "pinned_members": len(pinned), "pfba_cost": repr(total_cost), "priced_reactions": len(coef),
               "bpa_uptake": repr(uptake[bpa.id]), "glucose_uptake": repr(uptake.get(GLUCOSE_EXCHANGE, 0.0)),
               "carbon_in": repr(c_in), "carbon_out": repr(c_out), "carbon_closure": repr(closure), "free_source_max": repr(worst),
               "solver": solver, "threads": threads, "minutes": round((time.time() - t_start) / 60, 1)}
    with open(f"{out}/summary.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(summary), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerow(summary)
    log(f"5  community growth {growth:.6f} ({growth / required * 100:.2f} % of the required); {len(pinned)} pinned; "
        f"BPA uptake {uptake[bpa.id]:.6g}; carbon closure {closure:.3f}; done in {summary['minutes']} min -> {out}")


if __name__ == "__main__":
    main()
