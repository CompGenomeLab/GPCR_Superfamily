"""Consensus similarity of STE2 and STE3 to every other family on the ML-tree motif alignment
(representatives_all_Oct08_ginsi_ep0.123_treein; logo JSONs in GPCRevolution
public/superfamily_logo_mappings_treein, top-5 mapping).

For a family pair and a threshold t, the compared columns are those where both families have a
consensus: the most frequent residue (or similarity group) reaches a share >= t of all sequences and
the column is >= 50% occupied in both. The score is the fraction of compared columns whose consensus
is the same. Writes a table and a supplementary figure (one heatmap per query family).

    python ste_consensus_similarity.py                       # main alignment (--ep 0.123)
    python ste_consensus_similarity.py --mappings superfamily_logo_mappings_lower_gap_treein         --out ste_consensus_similarity_lower_gap               # plain G-INS-i, same ML guide tree
"""
import argparse
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ap = argparse.ArgumentParser()
ap.add_argument("--mappings", default="superfamily_logo_mappings_treein",
                help="mapping folder under GPCRevolution/public")
ap.add_argument("--out", default="ste_consensus_similarity", help="output name (tsv/pdf/png)")
args = ap.parse_args()
sys.path.insert(0, str(HERE.parent))
from gpcr_paths import stage  # noqa: E402
D = stage("10_motif_analysis", "consensus_mapping", "logo_conf0.33_top5") / args.mappings
OUT = stage("10_motif_analysis") / args.out
G = {a: g for g in ["ED", "WYHF", "RK", "QN", "VI", "ML", "ST"] for a in g}
# families in the order the Results text discusses them
FAM = ["cAMP", "classC", "classA", "Mth", "classB2", "classB1", "FSLB", "classF", "GP143", "TM116",
       "GP157", "Olfactory", "Vomeronasal1", "classT", "Nematode", "Vomeronasal2", "GPR1", "STE3", "STE2"]
NICE = {"classA": "Class A", "classB1": "Class B1", "classB2": "Class B2", "classC": "Class C",
        "classF": "Class F", "classT": "Taste-2", "FSLB": "FSL", "Vomeronasal1": "Vomeronasal 1",
        "Vomeronasal2": "Vomeronasal 2"}
THRESHOLDS = [0, 0.4, 0.6, 0.8]
QUERIES = ["STE2", "STE3"]

P = {k: json.load(open(D / f"{k}.json"))["positions"] for k in FAM}
n = len(P["classA"])


def consensus(p, grp):
    rc, t = p.get("residueCounts") or {}, p.get("totalSequences") or 0
    if not t or sum(rc.values()) / t < 0.5:
        return None, 0
    agg = {}
    for a, v in rc.items():
        k = G.get(a, a) if grp else a
        agg[k] = agg.get(k, 0) + v
    top = max(agg, key=agg.get)
    return top, agg[top] / t


C = {grp: {k: [consensus(P[k][c], grp) for c in range(n)] for k in FAM} for grp in (False, True)}
res = {}  # (query, family, grp, t) -> (matches, compared)
for q in QUERIES:
    for k in FAM:
        if k == q:
            continue
        for grp in (False, True):
            for t in THRESHOLDS:
                cols = [c for c in range(n) if C[grp][q][c][0] and C[grp][k][c][0]
                        and C[grp][q][c][1] >= t and C[grp][k][c][1] >= t]
                res[q, k, grp, t] = (sum(C[grp][q][c][0] == C[grp][k][c][0] for c in cols), len(cols))

with open(f"{OUT}.tsv", "w") as f:
    f.write("query\tfamily\tcomparison\tthreshold\tmatches\tcompared\tfraction\n")
    for (q, k, grp, t), (m, nc) in res.items():
        f.write(f"{q}\t{k}\t{'similarity' if grp else 'identity'}\t{t}\t{m}\t{nc}\t{m / nc if nc else 0:.3f}\n")

cols = [(grp, t) for grp in (False, True) for t in THRESHOLDS]
xlab = [("any" if t == 0 else f"≥{t:.0%}") for _, t in cols]
fig, axes = plt.subplots(1, 2, figsize=(11, 7.2), constrained_layout=True)
for ax, q, letter in zip(axes, QUERIES, "ab"):
    rows = [k for k in FAM if k != q]
    M = [[(res[q, k, g, t][0] / res[q, k, g, t][1]) if res[q, k, g, t][1] else float("nan") for g, t in cols]
         for k in rows]
    im = ax.imshow(M, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    for i, k in enumerate(rows):
        for j, (g, t) in enumerate(cols):
            m, nc = res[q, k, g, t]
            if nc:
                ax.text(j, i, f"{m / nc:.0%}\n{m}/{nc}", ha="center", va="center", fontsize=6.5,
                        color="white" if m / nc > 0.55 else "#222222")
            else:
                ax.text(j, i, "–", ha="center", va="center", fontsize=7, color="#888888")
    ax.set_xticks(range(len(cols)), xlab, fontsize=8)
    ax.set_yticks(range(len(rows)), [NICE.get(k, k) for k in rows], fontsize=8.5)
    ax.axvline(len(THRESHOLDS) - 0.5, color="white", lw=3)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(f"{letter}   {q}", loc="left", fontsize=11, fontweight="bold")
    for x, name in [(1.5, "Identical residue"), (5.5, "Same similarity group")]:
        ax.text(x, -0.9, name, ha="center", va="bottom", fontsize=8.5)
    ax.set_xlabel("Consensus share threshold", fontsize=8.5)
cb = fig.colorbar(im, ax=axes, shrink=0.5, pad=0.01)
cb.set_label("Fraction of compared columns with the same consensus", fontsize=8.5)
cb.ax.tick_params(labelsize=7.5)
for ext in ("pdf", "png"):
    fig.savefig(f"{OUT}.{ext}", dpi=300)
print(f"wrote {OUT}.tsv/.pdf/.png")
