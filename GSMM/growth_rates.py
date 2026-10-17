#!/usr/bin/env python3
"""Growth rate of every MAG at each simulated time point, from 16S rRNA gene qPCR and MAG abundances.

Over an interval of h hours the whole community grows at the rate

    r = mean over the three replicates of ln(N_end / N_start) / h        (N = 16S rRNA gene copies)

and the relative abundance a_i of MAG i (metaWRAP) changes with the difference between its own growth
rate and r, so

    mu_i = ln(a_i,end / a_i,start) / h + r

A negative mu_i is set to 0: a MAG whose abundance falls faster than the community grows is simulated as
not growing.

Input   16S rRNA gene copies per time point, three replicates (columns time_h, replicate_a, replicate_b,
        replicate_c) and the metaWRAP relative abundance of every MAG (one row per MAG, one column per
        time point: 0H, 24H, 26H, 30H).
Output  one row per simulated time point and MAG: the interval used, the two abundances, the community
        rate r, and mu_i before and after the negative values are set to 0.

Environment  Python 3.10, standard library only.

Usage   python growth_rates.py
        python growth_rates.py --qpcr 16s_qpcr.tsv --abundance mag_abundance.csv --out growth_rates.tsv
"""
import argparse
import csv
import math
import os
import sys

# --- Parameters --------------------------------------------------------------------------------
# Interval (hours) over which the growth rate of each simulated time point is computed.
GROWTH_INTERVALS = {"0h": (0, 24), "24h": (24, 30), "26h": (24, 30), "30h": (24, 30)}
# -----------------------------------------------------------------------------------------------


def read_qpcr(path):
    """{hours: [copies of replicate a, b, c]}"""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    copies = {int(r["time_h"]): [float(r[k]) for k in ("replicate_a", "replicate_b", "replicate_c")] for r in rows}
    if len(copies) != len(rows) or any(v <= 0 for vals in copies.values() for v in vals):
        sys.exit("STOP: the qPCR table must have one row per time and only positive counts")
    return copies


def read_abundance(path):
    """{MAG: {column: value as written in the file}}"""
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    key = list(rows[0])[0]                           # the first column holds the MAG names
    abundance = {r[key]: r for r in rows}
    if len(abundance) != len(rows):
        sys.exit("STOP: a MAG appears twice in the abundance table")
    return abundance


def community_rate(copies, start, end):
    hours = end - start
    return sum(math.log(n2 / n1) / hours for n1, n2 in zip(copies[start], copies[end])) / len(copies[start])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--qpcr", default="data/16s_qpcr.tsv", help="16S rRNA gene copies per time, three replicates")
    ap.add_argument("--abundance", default="data/mag_abundance.csv", help="metaWRAP abundance of every MAG, columns 0H, 24H, 26H, 30H")
    ap.add_argument("--out", default="results/growth_rates.tsv")
    a = ap.parse_args()
    intervals = GROWTH_INTERVALS

    copies = read_qpcr(a.qpcr)
    abundance = read_abundance(a.abundance)
    times = {t for interval in intervals.values() for t in interval}
    if times - set(copies):
        sys.exit(f"STOP: the qPCR table has no row for {sorted(times - set(copies))} h")
    columns = {f"{t}H" for t in times} - set(next(iter(abundance.values())))
    if columns:
        sys.exit(f"STOP: the abundance table has no column {sorted(columns)}")
    os.makedirs(os.path.dirname(a.out) or ".", exist_ok=True)
    with open(a.out, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["time_point", "mag", "interval", "abundance_start", "abundance_end",
                    "community_rate", "growth_rate_uncapped", "growth_rate"])
        for tp, (start, end) in intervals.items():
            r = community_rate(copies, start, end)
            growing = 0
            for mag in sorted(abundance, key=lambda m: int(m.split(".")[1])):
                a1, a2 = abundance[mag][f"{start}H"], abundance[mag][f"{end}H"]
                if float(a1) <= 0 or float(a2) <= 0:
                    sys.exit(f"STOP: {mag} has a non-positive abundance at {start} or {end} h")
                mu = math.log(float(a2) / float(a1)) / (end - start) + r
                growing += mu > 0
                w.writerow([tp, mag, f"{start}-{end}h", a1, a2, f"{r:.6g}", f"{mu:.6g}", f"{max(0.0, mu):.6g}"])
            print(f"{tp:>6}: interval {start}-{end} h, community rate {r:.5f} /h, "
                  f"{growing} of {len(abundance)} MAGs with a positive growth rate")
    print(f"Done: {len(intervals) * len(abundance)} rows written to {a.out}")


if __name__ == "__main__":
    main()
