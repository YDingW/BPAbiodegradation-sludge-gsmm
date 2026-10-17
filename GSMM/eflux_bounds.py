#!/usr/bin/env python3
"""Community members, expression coefficients and E-Flux bounds for every simulated time point, from the
metatranscriptome.

Expression of a gene at a time point = square root of the mean of its three replicates.

1. Members. A MAG is a community member at a time point when transcripts (mean > 0) are found for at least
   DETECTION_CUTOFF % of its protein-coding genes; bin.47 and bin.104 fall below the cutoff at 0 h (1.5 % and
   3.1 %) but are high at 24, 26 and 30 h (38-66 %), so they are kept as members at 0 h as well (see FORCED).
2. Reference. G is the Q_REFERENCE-th percentile of the expression of all detected protein-coding genes of all
   MAGs over the four time points; one G for all time points.
3. Coefficient of a reaction = min(1, score / G), where the score is its GPR evaluated on gene expression
   (maximum over OR, minimum over AND, genes without transcript data ignored). It is not allowed below the
   floor of the MAG, the Q_FLOOR-th percentile of min(1, expression / G) over the MAG's detected
   protein-coding genes. Reactions with a GPR but no data and metabolic reactions without a GPR get the
   floor. Exempt, with bounds unchanged: the biomass reaction, spontaneous reactions, and exchange and
   transport reactions without a GPR.
4. E-Flux bounds: every non-zero bound is multiplied by the coefficient (the lower bound if negative, the
   upper bound if positive), so an irreversible reaction keeps its direction.

Input   the metatranscriptome table (one row per gene, one column per sample), the curated models (*.xml) and
        the protein FASTA files the models were built from (*.faa, one per MAG).
Output  per time point, <out>/<time point>_mags.tsv (every MAG: detection, membership, floor) and
        <out>/<time point>_reactions.tsv (every constrained reaction of every member: class, source,
        coefficient, model bounds and E-Flux bounds); <out>/summary.tsv has G and the member counts.

Environment  Python 3.10, cobra 0.31.1, pandas 2.3.3, numpy 2.2.6.

Usage   python eflux_bounds.py
        python eflux_bounds.py --metat metaT.csv --models curated --proteomes proteomes --out eflux
"""
import argparse
import csv
import logging
import os
import re
import sys
import time

import cobra
import numpy as np
import pandas as pd

# --- Parameters --------------------------------------------------------------------------------
TIME_POINTS = ["0h", "24h", "26h", "30h"]          # simulated time points (hours after the BPA spike)
REPLICATES = ["a", "b", "c"]                       # the three metatranscriptome replicates
SAMPLE_COLUMN = "12C-{replicate}-{time_point}"     # column of replicate `replicate` at time point `time_point`
DETECTION_CUTOFF = 5.0                             # % of protein-coding genes with transcripts to be a member
Q_REFERENCE = 99                                   # percentile of the reference expression G
Q_FLOOR = 10                                       # percentile of the floor of each MAG
FORCED = {"0h": ["bin.47", "bin.104"], "24h": [], "26h": [], "30h": []}
# FORCED: members regardless of the cutoff (see point 1 above).
SAMPLES = {tp: [SAMPLE_COLUMN.format(replicate=r, time_point=tp) for r in REPLICATES] for tp in TIME_POINTS}
# -----------------------------------------------------------------------------------------------

# source of a coefficient: expression = from the GPR score; below_floor = score below the floor;
# no_data = GPR without transcript data; no_gpr = metabolic reaction without GPR


def gpr_score(rule, expression):
    """Expression of a GPR: max over OR, min over AND; genes without data are ignored. None when nothing has data."""
    if not rule.strip() or rule.strip().lower() == "spontaneous":
        return None
    tokens = re.findall(r"\(|\)|\band\b|\bor\b|[^\s()]+", rule, flags=re.I)
    pos = 0

    def or_expr():
        nonlocal pos
        values = [and_expr()]
        while pos < len(tokens) and tokens[pos].lower() == "or":
            pos += 1
            values.append(and_expr())
        values = [v for v in values if v is not None]
        return max(values) if values else None

    def and_expr():
        nonlocal pos
        values = [atom()]
        while pos < len(tokens) and tokens[pos].lower() == "and":
            pos += 1
            values.append(atom())
        values = [v for v in values if v is not None]
        return min(values) if values else None

    def atom():
        nonlocal pos
        if pos >= len(tokens):
            return None
        token = tokens[pos]
        pos += 1
        if token == "(":
            value = or_expr()
            if pos < len(tokens) and tokens[pos] == ")":
                pos += 1
            return value
        return expression.get(token)

    return or_expr()


def classify(rxn):
    """Reaction class and whether it has a GPR."""
    rule = rxn.gene_reaction_rule.strip()
    if rxn.objective_coefficient != 0 or rxn.id == "Growth" or "biomass" in rxn.id.lower():
        return "biomass", bool(rule)
    if rule.lower() == "spontaneous":
        return "spontaneous", False
    if rxn.id.startswith("EX_"):
        return "exchange", bool(rule)
    if len({m.compartment for m in rxn.metabolites}) > 1:
        return "transport", bool(rule)
    return "metabolic", bool(rule)


def protein_ids(faa):
    """Gene IDs of a protein FASTA file, one per header: the first word, after a '|' if there is one."""
    ids = []
    with open(faa) as f:
        for line in f:
            if line.startswith(">"):
                word = line[1:].split()[0]
                ids.append(word.split("|", 1)[1] if "|" in word else word)
    return ids


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--metat", default="data/metaT.csv", help="metatranscriptome table")
    ap.add_argument("--models", default="models/curated", help="folder with the curated models (*.xml)")
    ap.add_argument("--proteomes", default="models/proteomes", help="folder with the protein FASTA files (*.faa)")
    ap.add_argument("--out", default="results/eflux")
    a = ap.parse_args()
    logging.getLogger("cobra").setLevel(logging.ERROR)
    os.makedirs(a.out, exist_ok=True)
    t0 = time.time()

    mags = sorted((f[:-4] for f in os.listdir(a.models) if f.endswith(".xml")), key=lambda m: int(m.split(".")[1]))
    missing = [m for m in mags if not os.path.exists(f"{a.proteomes}/{m}.faa")]
    if missing:
        sys.exit(f"STOP: no protein file for {missing[:5]}")
    proteins = {m: protein_ids(f"{a.proteomes}/{m}.faa") for m in mags}
    protein_sets = {m: set(ids) for m, ids in proteins.items()}
    forced = {tp: FORCED.get(tp, []) for tp in TIME_POINTS}
    unknown = sorted({m for tp in TIME_POINTS for m in forced[tp]} - set(mags))
    if unknown:
        sys.exit(f"STOP: forced member {unknown} is not in {a.models}")

    mt = pd.read_csv(a.metat)
    columns = [c for tp in TIME_POINTS for c in SAMPLES[tp] if c not in mt.columns]
    if columns:
        sys.exit(f"STOP: {a.metat} has no column {columns}")
    mt["mag"] = mt["bin_gene"].astype(str).str.split("|").str[0]
    mt["gene"] = mt["bin_gene"].astype(str).str.split("|").str[1].astype(str)
    is_protein = np.array([g in protein_sets.get(m, ()) for m, g in zip(mt["mag"], mt["gene"])])
    in_models = mt["mag"].isin(mags).to_numpy()
    if in_models.sum() == 0 or is_protein[in_models].mean() < 0.95:
        sys.exit("STOP: the gene IDs of the metatranscriptome do not match the protein files")
    expression = {tp: mt[SAMPLES[tp]].apply(pd.to_numeric, errors="coerce").fillna(0.0).mean(axis=1).to_numpy()
                  for tp in TIME_POINTS}
    detected = [np.sqrt(x[in_models & is_protein & (x > 0)]) for x in expression.values()]
    G = float(np.percentile(np.concatenate(detected), Q_REFERENCE))
    print(f"reference G = {G:.4f} ({sum(len(d) for d in detected)} detected protein-coding genes over {len(TIME_POINTS)} time points)")

    mag_rows = {tp: [] for tp in TIME_POINTS}
    rxn_rows = {tp: [] for tp in TIME_POINTS}
    for i, m in enumerate(mags, 1):
        rows = (mt["mag"] == m).to_numpy()
        genes = mt.loc[rows, "gene"].tolist()
        model = None
        for tp in TIME_POINTS:
            x = expression[tp]
            n_detected = int((x[rows] > 0).sum())
            pct = n_detected / len(proteins[m]) * 100 if proteins[m] else 0.0
            member = "forced" if m in forced[tp] else ("yes" if pct >= DETECTION_CUTOFF else "no")
            record = {"mag": m, "protein_genes": len(proteins[m]), "detected_genes": n_detected, "detected_pct": repr(pct),
                      "member": member, "floor": "", "reactions_with_data": "", "reactions_no_data": "", "reactions_no_gpr": ""}
            if member != "no":
                if model is None:
                    model = cobra.io.read_sbml_model(f"{a.models}/{m}.xml")
                sqrt_x = np.sqrt(x)
                expr = dict(zip(genes, sqrt_x[rows].tolist()))
                scores = {r.id: s for r in model.reactions for s in [gpr_score(r.gene_reaction_rule, expr)] if s is not None}
                if not scores or max(scores.values()) <= 0:
                    sys.exit(f"STOP: {m} is a member at {tp} but none of its reactions has a GPR score")
                floor = float(np.percentile(np.minimum(1.0, sqrt_x[rows & is_protein & (x > 0)] / G), Q_FLOOR))
                counts = {"expression": 0, "below_floor": 0, "no_data": 0, "no_gpr": 0}
                for r in model.reactions:
                    cls, has_gpr = classify(r)
                    if cls in ("biomass", "spontaneous") or (cls in ("exchange", "transport") and not has_gpr):
                        continue
                    if not has_gpr:
                        coef, source = floor, "no_gpr"
                    elif r.id not in scores:
                        coef, source = floor, "no_data"
                    else:
                        coef = min(1.0, scores[r.id] / G)
                        if coef <= floor:
                            coef, source = floor, "below_floor"
                        else:
                            source = "expression"
                    counts[source] += 1
                    lb, ub = r.lower_bound, r.upper_bound
                    rxn_rows[tp].append([m, r.id, cls, "yes" if has_gpr else "no", source, repr(coef), repr(lb), repr(ub),
                                         repr(float(lb * coef) if lb < 0 else lb), repr(float(ub * coef) if ub > 0 else ub)])
                record.update(floor=repr(floor), reactions_with_data=counts["expression"] + counts["below_floor"],
                              reactions_no_data=counts["no_data"], reactions_no_gpr=counts["no_gpr"])
            mag_rows[tp].append(record)
        print(f"[{i:>3}/{len(mags)}] {m:<8} member at {sum(r[-1]['member'] != 'no' for r in mag_rows.values())} time points   {time.time() - t0:6.0f} s", flush=True)

    summary = []
    for tp in TIME_POINTS:
        with open(f"{a.out}/{tp}_mags.tsv", "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(mag_rows[tp][0]), delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(mag_rows[tp])
        with open(f"{a.out}/{tp}_reactions.tsv", "w", newline="") as f:
            w = csv.writer(f, delimiter="\t", lineterminator="\n")
            w.writerow(["mag", "reaction", "class", "gpr", "source", "coefficient", "lower_bound", "upper_bound",
                        "eflux_lower_bound", "eflux_upper_bound"])
            w.writerows(rxn_rows[tp])
        members = sum(r["member"] != "no" for r in mag_rows[tp])
        summary.append([tp, members, len(forced[tp]), len(rxn_rows[tp]), repr(G)])
        print(f"{tp:>6}: {members} members ({len(forced[tp])} forced), {len(rxn_rows[tp])} constrained reactions")
    with open(f"{a.out}/summary.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["time_point", "members", "forced_members", "constrained_reactions", "reference_G"])
        w.writerows(summary)
    print(f"Done: tables written to {a.out}")


if __name__ == "__main__":
    main()
