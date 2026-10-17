#!/usr/bin/env python3
"""Metabolite exchange tables from the community simulations of the four time points.

For every member and metabolite the exchange reaction of the member (EX_<metabolite>_e) gives the specific flux
(mmol per gDW of the member per h; negative = uptake, positive = secretion). The weighted flux is the specific
flux times the member's abundance share (mmol per gDW of community per h) and the carbon flux the weighted flux
times the carbon atoms of the metabolite. The weighted fluxes of all members sum to the community's exchange
with the medium, which is checked.

The SCR members (SCR below) are the 10 MAGs of the families Sphingomonadaceae, Comamonadaceae and
Rhizobiaceae, the three families identified as keystone BPA degraders in the paper.

Tables written to <out>:
  members.tsv              every member at every time point: abundance, growth rate, metabolites secreted / taken up
  exchange_members.tsv     every member x metabolite with a non-zero specific flux
  exchange_metabolites.tsv every metabolite at every time point: total secretion and uptake over the members, the
                           cross-fed amount min(secretion, uptake) and its carbon, the net exchange with the medium,
                           numbers of secreting / taking members, and the net flux of the SCR members
  cross_feeding_pairs.tsv  donor -> receiver pairs for organic metabolites and siderophores. A metabolite secreted by
                           several members and taken up by several others goes through a shared pool, so a pair gets
                           the cross-fed amount in proportion to the donor's share of the secretion and the receiver's
                           share of the uptake; surplus secretion goes to "medium" and uptake not covered by secretion
                           comes from "medium".
  scr_patterns.tsv         per metabolite, the net weighted flux of the SCR members at the four time points and its
                           pattern: A = uptake at 0 h and secretion at the three later points, B = none at 0 h and
                           secretion at the three later points; also the gross pattern (any SCR member secreting /
                           taking up): A'' = no SCR member secretes at 0 h but one takes up, and SCR members secrete
                           at the three later points; B'' = no SCR exchange at 0 h, secretion at the three later points.

Input   one sub-folder per time point with the fluxes.tsv and growth.tsv of the community simulation, and a
        table of metabolite classes, names and carbon numbers.
Output  the five tables above.

Environment  Python 3.10, pandas 2.3.3.

Usage   python exchange_tables.py --runs community --out community/exchange
"""
import argparse
import csv
import os
import sys

import pandas as pd

# --- Parameters --------------------------------------------------------------------------------
TIME_POINTS = ["0h", "24h", "26h", "30h"]
# The 10 MAGs of the three keystone families (SCR: Sphingomonadaceae, Comamonadaceae, Rhizobiaceae)
SCR_MEMBERS = ["bin.47", "bin.88", "bin.97", "bin.104", "bin.105", "bin.110", "bin.81", "bin.91", "bin.36", "bin.78"]
# -----------------------------------------------------------------------------------------------

EPS = 1e-6            # a weighted flux with |w| <= EPS counts as zero
FLUX_MIN = 1e-9       # a member x metabolite row needs |specific flux| >= FLUX_MIN
EXCHANGE_CLASSES = ["amino_acid", "peptide", "organic_acid", "vitamin", "sugar", "nucleoside", "other", "bpa_pathway", "siderophore"]
CLASS_ORDER = EXCHANGE_CLASSES + ["inorganic"]


def proportional_pairs(secretion, uptake):
    """{donor: amount}, {receiver: amount} -> [(donor, receiver, amount)] through a shared pool."""
    s_total, u_total = sum(secretion.values()), sum(uptake.values())
    rows = []
    if s_total <= 0 and u_total <= 0:
        return rows
    transfer = min(s_total, u_total)
    if s_total > 0 and u_total > 0:
        for donor, s in secretion.items():
            for receiver, u in uptake.items():
                rows.append((donor, receiver, transfer * (s / s_total) * (u / u_total)))
    if s_total - transfer > FLUX_MIN:              # surplus secretion leaves to the medium
        for donor, s in secretion.items():
            rows.append((donor, "medium", (s_total - transfer) * s / s_total))
    if u_total - transfer > FLUX_MIN:              # uptake not covered by secretion comes from the medium
        for receiver, u in uptake.items():
            rows.append(("medium", receiver, (u_total - transfer) * u / u_total))
    return [r for r in rows if r[2] > 0]


def pattern(values):
    """Net weighted SCR flux at the four time points -> (pattern, direction at 0 h, later secretions, later uptakes)."""
    v0, late = values[0], values[1:]
    d0 = "uptake" if v0 < -EPS else ("secretion" if v0 > EPS else "none")
    n_sec, n_upt = sum(v > EPS for v in late), sum(v < -EPS for v in late)
    if n_sec == 3 and d0 == "uptake":
        return "A", d0, n_sec, n_upt
    if n_sec == 3 and d0 == "none":
        return "B", d0, n_sec, n_upt
    return "-", d0, n_sec, n_upt


def gross_pattern(gross):
    """[(gross secretion, gross uptake)] at the four time points -> A'', B'' or -."""
    s0, u0 = gross[0]
    late = all(s > EPS for s, u in gross[1:])
    if late and s0 <= EPS and u0 > EPS:
        return "A''"
    if late and s0 <= EPS and u0 <= EPS:
        return "B''"
    return "-"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--runs", default="results/community", help="folder with one sub-folder per time point (fluxes.tsv, growth.tsv)")
    ap.add_argument("--out", default=None, help="default <runs>/exchange")
    ap.add_argument("--classes", default="data/metabolite_classes.tsv")
    a = ap.parse_args()
    tps = TIME_POINTS
    scr = set(SCR_MEMBERS)
    out = a.out or f"{a.runs}/exchange"
    os.makedirs(out, exist_ok=True)
    classes = pd.read_csv(a.classes, sep="\t", index_col="metabolite")
    cls_of = lambda m: classes.at[m, "class"] if m in classes.index else "unclassified"
    name_of = lambda m: classes.at[m, "name"] if m in classes.index else m
    carbons_of = lambda m: float(classes.at[m, "carbons"]) if m in classes.index else 0.0
    order = {c: i for i, c in enumerate(CLASS_ORDER)}

    members_rows, member_rows, met_rows, pair_rows = [], [], [], []
    scr_net, scr_gross, present = {}, {}, {tp: set() for tp in tps}
    for tp in tps:
        d = f"{a.runs}/{tp}"
        if not os.path.exists(f"{d}/fluxes.tsv"):
            sys.exit(f"STOP: {d}/fluxes.tsv not found")
        fluxes = pd.read_csv(f"{d}/fluxes.tsv", sep="\t", index_col=0)
        growth = pd.read_csv(f"{d}/growth.tsv", sep="\t", index_col="member")
        members = [m for m in fluxes.index if m != "medium"]
        if sorted(members) != sorted(growth.index):
            sys.exit(f"STOP: {tp}: the members of fluxes.tsv and growth.tsv differ")
        share = growth["abundance"]
        if abs(share.sum() - 1) > 1e-9:
            sys.exit(f"STOP: {tp}: abundance shares do not sum to 1")
        ex = [c for c in fluxes.columns if c.startswith("EX_") and c.endswith("_e")]
        V = fluxes.loc[members, ex].fillna(0.0)
        V.columns = [c[3:-2] for c in ex]
        W = V.mul(share, axis=0)
        medium = {c[3:-2]: float(fluxes.at["medium", c]) for c in fluxes.columns if c.startswith("EX_") and c.endswith("_m")}
        present[tp] = set(members)
        for m in members:
            members_rows.append({"time_point": tp, "member": m, "scr": "yes" if m in scr else "no", "abundance": share[m],
                                 "growth_rate": growth.at[m, "growth_rate"], "metabolites_secreted": int((V.loc[m] > 1e-6).sum()),
                                 "metabolites_taken_up": int((V.loc[m] < -1e-6).sum())})
        for met in V.columns:
            v, w = V[met], W[met]
            total, med = float(w.sum()), medium.get(met)
            if med is not None and abs(total - med) > 1e-6:
                sys.exit(f"STOP: {tp} {met}: the weighted member fluxes sum to {total:.6g} but the medium exchange is {med:.6g}")
            cls, carbons = cls_of(met), carbons_of(met)
            for m in members:
                if abs(v[m]) < FLUX_MIN:
                    continue
                member_rows.append({"time_point": tp, "member": m, "scr": "yes" if m in scr else "no", "abundance": share[m], "metabolite": met,
                                    "name": name_of(met), "class": cls, "carbons": carbons, "direction": "secretion" if v[m] > 0 else "uptake",
                                    "specific_flux": v[m], "weighted_flux": w[m], "carbon_flux": w[m] * carbons})
            w_scr = w[[m for m in members if m in scr]]
            if cls != "inorganic":
                scr_net.setdefault(met, {})[tp] = float(w_scr.sum())
                scr_gross.setdefault(met, {})[tp] = (float(w_scr.clip(lower=0).sum()), float((-w_scr).clip(lower=0).sum()))
            if (w.abs() < 1e-12).all() and abs(med or 0.0) < 1e-12:
                continue
            sec, upt = float(w[w > 0].sum()), float(-w[w < 0].sum())
            met_rows.append({"time_point": tp, "metabolite": met, "name": name_of(met), "class": cls, "carbons": carbons,
                             "secretion": sec, "uptake": upt, "cross_fed": min(sec, upt), "cross_fed_carbon": min(sec, upt) * carbons,
                             "medium_net": med if med is not None else float("nan"), "secreting_members": int((w > 1e-9).sum()),
                             "taking_members": int((w < -1e-9).sum()), "scr_net": float(w_scr.sum())})
            if cls in EXCHANGE_CLASSES:
                secretion = {m: float(w[m]) for m in members if w[m] > 1e-9}
                uptake = {m: -float(w[m]) for m in members if w[m] < -1e-9}
                for donor, receiver, amount in proportional_pairs(secretion, uptake):
                    pair_rows.append({"time_point": tp, "metabolite": met, "name": name_of(met), "class": cls, "carbons": carbons,
                                      "donor": donor, "receiver": receiver, "weighted_flux": amount, "carbon_flux": amount * carbons,
                                      "donors": len(secretion), "receivers": len(uptake), "cross_fed": min(sum(secretion.values()), sum(uptake.values()))})
    scr_rows = []
    for met, d in scr_net.items():
        values = [d.get(tp, 0.0) for tp in tps]
        if all(abs(x) <= EPS for x in values):
            continue
        p, d0, n_sec, n_upt = pattern(values)
        gross = [scr_gross[met].get(tp, (0.0, 0.0)) for tp in tps]
        scr_rows.append({"metabolite": met, "name": name_of(met), "class": cls_of(met), "carbons": carbons_of(met),
                         **{f"scr_net_{tp}": d.get(tp, float("nan")) for tp in tps}, "pattern": p, "direction_0h": d0,
                         "later_secretions": n_sec, "later_uptakes": n_upt,
                         **{f"scr_secretion_{tp}": gross[i][0] for i, tp in enumerate(tps)}, **{f"scr_uptake_{tp}": gross[i][1] for i, tp in enumerate(tps)},
                         "gross_pattern": gross_pattern(gross)})
    tables = {"members.tsv": pd.DataFrame(members_rows), "exchange_members.tsv": pd.DataFrame(member_rows),
              "exchange_metabolites.tsv": pd.DataFrame(met_rows), "cross_feeding_pairs.tsv": pd.DataFrame(pair_rows),
              "scr_patterns.tsv": pd.DataFrame(scr_rows)}
    for name, T in tables.items():
        if len(T) and "class" in T.columns:
            T["_c"] = T["class"].map(lambda c: order.get(c, 99))
            if "time_point" in T.columns:
                T["_t"] = T["time_point"].map(tps.index)
                T = T.sort_values(["_t", "_c", "metabolite"] + (["weighted_flux"] if "weighted_flux" in T.columns else []),
                                  ascending=[True, True, True] + ([False] if "weighted_flux" in T.columns else []))
            else:
                T = T.sort_values(["pattern", "_c", "metabolite"])
            T = T.drop(columns=[c for c in ("_c", "_t") if c in T.columns])
        T.to_csv(f"{out}/{name}", sep="\t", index=False, float_format="%.15g")
        print(f"{name:<26} {len(T):>6} rows")
    print(f"Done: tables written to {out}")


if __name__ == "__main__":
    main()
