"""Sankey diagram of the MAGs across abundance cluster -> family -> expression cluster.
Input:  results/clusters_abundance_vs_expression.csv (written by mfuzz_clustering.R) and data/mag_taxonomy.csv
        (MAG, family: the three families highlighted in the paper, every other MAG is "Others").
Output: results/sankey_plot.png, 600 dpi. Ribbon widths are numbers of MAGs, colours are the families.
Environment: Python 3.13.16 with pandas 2.3.3 and matplotlib 3.11.2.
Run from this folder:  python sankey_plot.py
"""
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import PathPatch, Rectangle
from matplotlib.path import Path

matplotlib.rcParams.update({"font.family": ["Arial", "Liberation Sans", "DejaVu Sans"], "font.weight": "bold"})

clusters = pd.read_csv("results/clusters_abundance_vs_expression.csv").dropna(subset=["expression_cluster"])   # MAGs in both clusterings
df = clusters.merge(pd.read_csv("data/mag_taxonomy.csv"), on="MAG")
df["abundance_cluster"] = df["abundance_cluster"].astype(int)
df["expression_cluster"] = df["expression_cluster"].astype(int)

families = ["Sphingomonadaceae", "Rhizobiaceae", "Comamonadaceae", "Others"]
colour = {"abundance": {1: "#A6CEE3", 2: "#B2DF8A", 3: "#CAB2D6"}, "expression": {1: "#FDBF6F", 2: "#FB9A99", 3: "#D9D9D9"},
          "family": {"Sphingomonadaceae": "#E64B35", "Rhizobiaceae": "#4DBBD5", "Comamonadaceae": "#00A087", "Others": "#BDBDBD"}}
BAR_W, GAP = 0.12, 2.5
X = {"abundance": 1.10, "family": 2.00, "expression": 2.90}


def stack(sizes):
    """Top and bottom of one bar per category, stacked from the top with gaps; y unit = one MAG."""
    pos, y = {}, sizes.sum() + GAP * (len(sizes) - 1)
    for k, n in sizes.items():
        pos[k] = (y - n, y)
        y -= n + GAP
    return pos


sizes = {"abundance": df["abundance_cluster"].value_counts().reindex([1, 2, 3], fill_value=0),
         "family": df["family"].value_counts().reindex(families, fill_value=0),
         "expression": df["expression_cluster"].value_counts().reindex([1, 2, 3], fill_value=0)}
height = max(s.sum() + GAP * (len(s) - 1) for s in sizes.values())
pos = {}
for col, s in sizes.items():                                   # centre the three columns vertically
    shift = (height - s.sum() - GAP * (len(s) - 1)) / 2
    pos[col] = {k: (b + shift, t + shift) for k, (b, t) in stack(s).items()}

fig, ax = plt.subplots(figsize=(7.5, 8))


def ribbon(x0, y0, x1, y1, n, c, z):
    """Band of n MAGs from the top y0 of a bar at x0 to the top y1 of a bar at x1."""
    cx = (x0 + x1) / 2
    verts = [(x0, y0), (cx, y0), (cx, y1), (x1, y1), (x1, y1 - n), (cx, y1 - n), (cx, y0 - n), (x0, y0 - n), (x0, y0)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), facecolor=c, edgecolor="none", alpha=0.55, lw=0, zorder=z))


top = {col: {k: t for k, (b, t) in p.items()} for col, p in pos.items()}       # next free y at the top of every bar
top["family_right"] = dict(top["family"])
order = {f: i for i, f in enumerate(families)}
left = df.groupby(["abundance_cluster", "family"]).size().reset_index(name="n")
for _, r in left.sort_values(["abundance_cluster", "family"], key=lambda s: s.map(order) if s.name == "family" else s).iterrows():
    a, f, n = r["abundance_cluster"], r["family"], r["n"]
    ribbon(X["abundance"] + BAR_W / 2, top["abundance"][a], X["family"] - BAR_W / 2, top["family"][f], n, colour["family"][f], 3 if f == "Others" else 4)
    top["abundance"][a] -= n; top["family"][f] -= n
right = df.groupby(["family", "expression_cluster"]).size().reset_index(name="n")
for _, r in right.sort_values(["family", "expression_cluster"], key=lambda s: s.map(order) if s.name == "family" else s).iterrows():
    f, e, n = r["family"], r["expression_cluster"], r["n"]
    ribbon(X["family"] + BAR_W / 2, top["family_right"][f], X["expression"] - BAR_W / 2, top["expression"][e], n, colour["family"][f], 3 if f == "Others" else 4)
    top["family_right"][f] -= n; top["expression"][e] -= n

for col, label in (("abundance", "AC"), ("family", None), ("expression", "EC")):
    for k, (b, t) in pos[col].items():
        ax.add_patch(Rectangle((X[col] - BAR_W / 2, b), BAR_W, t - b, facecolor=colour[col][k], edgecolor="black", lw=0.8, zorder=5))
        if label:
            ax.text(X[col], (b + t) / 2, f"{label}{k}", ha="center", va="center", rotation=90, fontsize=11, fontweight="bold", zorder=6)
ax.set_xlim(X["abundance"] - 0.6, X["expression"] + 0.6)
ax.set_ylim(-3, height + 10)
ax.axis("off")
plt.tight_layout()
plt.savefig("results/sankey_plot.png", dpi=600, bbox_inches="tight", facecolor="white")
