#!/usr/bin/env python3
"""
discover_groupings.py - find every multi-family split in the ML trees without naming
clades in advance, and report how often each is recovered and what sCF it carries.

    python discover_groupings.py <dir_of_cf_trees> [<dir2> ...] -o grouping_discovery \
           --labels L1 L2 ... --divergent STE2 STE3 classC --mode excluded|included

SPLIT OF A BRANCH
    Every internal branch cuts the unrooted tree into two sides. Each family is assigned
    to the side that holds MORE THAN 50% of its leaves in that tree (same rule as
    gather_scf.py), so a branch becomes a split of the FAMILIES into two
    groups, written with BOTH sides, the group with fewer families first:
    "GP143+TM116 | cAMP+classA+...". sCF scores the branch, i.e. both sides at once.
    Naming does not depend on how many leaves each side has, so the same split gets the
    same name in every dataset. Nested branches that give the same split (e.g. differing
    only by a few stray leaves) are collapsed; per tree the highest-sCF one is kept (biased
    upward, but equally in every dataset).

--mode excluded
    The --divergent families (default STE2, STE3, classC; Vomeronasal2 is merged into
    classC in the family map) are left out when the split is named, so a split is the same
    whether or not those long-branch families sit inside it. The column
    "divergent_placement" records, per tree, which divergent families fell on side A
    ("in: ...") - this is where they attach.
--mode included
    Every family counts, divergent ones included. Splits are named as they are.

WHAT IS REPORTED
    Splits whose named group has >= --min-families families. A split is listed if, in at
    least one dataset, it is recovered in >= --min-trees trees with a median sCF above
    --threshold. For every listed split all datasets are shown, including those where it
    is rare or weak. Rows are ordered by the best median sCF among the datasets where the
    split passes (most supported first).

Self-contained: standard library only. matplotlib optional, for the figure.
sCF noise floor is ~33%.
"""
import os, sys, glob, argparse, statistics as st
from collections import defaultdict, Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from gather_scf import parse_newick, accession, scf_of, tree_index, ACC2FAM   # noqa
from gather_scf import MAJORITY as MAJ   # noqa

TAB, NL = "\t", "\n"
FAM_ORDER = ["classA", "Nematode", "Olfactory", "classT", "Vomeronasal1", "GPR1",
             "cAMP", "GP157", "GP143", "TM116", "classB1", "classB2", "Mth", "classF", "FSLB",
             "STE2", "STE3", "classC"]
DISPLAY = {"classA": "Class A", "classB1": "Class B1", "classB2": "Class B2",
           "classC": "Class C", "classF": "Class F", "classT": "Class T"}


def branch_family_sets(path):
    """-> (branches, ftot) where branches is a list of
       (scf, Counter_of_family_counts_on_smaller_side, n_tips)"""
    tips, internals = parse_newick(open(path, encoding='utf-8').read())
    tip_fam, unknown = {}, set()
    for t in tips:
        f = ACC2FAM.get(accession(t))
        if f is None:
            unknown.add(accession(t))
        tip_fam[t] = f
    if unknown:
        raise SystemExit("ERROR " + path + ": " + str(len(unknown)) +
                         " tips not in family map, e.g. " +
                         str(sorted(x for x in unknown if x)[:5]))
    ftot = Counter(tip_fam.values())
    N, allt = len(tips), set(tips)
    out = []
    for label, clade in internals:
        v = scf_of(label)
        if v is None:
            continue
        small = clade if len(clade) * 2 <= N else (allt - clade)
        if len(small) < 2 or len(small) > N - 2:
            continue
        out.append((v, Counter(tip_fam[t] for t in small), len(small)))
    return out, ftot


def fam_key(fams):
    return "+".join(sorted(fams, key=str.lower))


def sides(key):
    a, b = key.split(" | ")
    return a.split("+"), b.split("+")


def splits(path, min_fam, divergent, mode):
    """-> {split name: (best sCF, divergent placement on that branch)} for one tree"""
    branches, ftot = branch_family_sets(path)
    ignore = divergent if mode == "excluded" else frozenset()
    universe = frozenset(f for f in ftot if f not in ignore)
    best = {}
    for v, counts, nt in branches:
        # counts are for one side; family goes to the side holding > 50% of its leaves
        side = frozenset(f for f in universe if counts.get(f, 0) > MAJ * ftot[f])
        other = frozenset(f for f in universe
                          if ftot[f] - counts.get(f, 0) > MAJ * ftot[f])
        if not side or not other:
            continue                      # branch inside one group of families
        # name by the group with fewer families; tie -> alphabetical
        a = min((side, other), key=lambda s: (len(s), fam_key(s).lower()))
        if len(a) < min_fam:
            continue
        key = fam_key(a) + " | " + fam_key(side if a == other else other)
        place = ""
        if mode == "excluded":
            named_is_side = (a == side)
            inside = sorted(f for f in divergent if f in ftot and
                            ((counts.get(f, 0) if named_is_side else ftot[f] - counts.get(f, 0))
                             > MAJ * ftot[f]))
            place = ("in: " + "+".join(inside)) if inside else "none in"
        if key not in best or v > best[key][0]:
            best[key] = (v, place)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dirs', nargs='+')
    ap.add_argument('-o', '--out', default='grouping_discovery')
    ap.add_argument('--labels', nargs='*', default=None)
    ap.add_argument('--mode', choices=['excluded', 'included'], default='excluded')
    ap.add_argument('--divergent', nargs='*', default=["STE2", "STE3", "classC"])
    ap.add_argument('--min-families', type=int, default=2)
    ap.add_argument('--min-trees', type=int, default=10,
                    help='recovered in at least this many trees in some dataset')
    ap.add_argument('--plot-min-share', type=float, default=0.5,
                    help='figure only: show splits recovered in at least this share of the '
                         'trees in at least one dataset (tables keep every split)')
    ap.add_argument('--threshold', type=float, default=33.0,
                    help='median sCF must be above this in that dataset')
    args = ap.parse_args()
    divergent = frozenset(args.divergent)
    excl = args.mode == "excluded"

    res = defaultdict(lambda: defaultdict(list))      # ds -> split -> [sCF per tree]
    plc = defaultdict(lambda: defaultdict(list))      # ds -> split -> [placement per tree]
    dsnames, ntrees = [], {}
    for d in args.dirs:
        files = sorted(glob.glob(os.path.join(d, '*.cf.tree')), key=tree_index)
        if not files:
            print("  !! no *.cf.tree in " + d, file=sys.stderr)
            continue
        ds = os.path.basename(os.path.normpath(d))
        dsnames.append(ds)
        ntrees[ds] = len(files)
        print(ds + ": " + str(len(files)) + " trees")
        for p in files:
            for key, (v, place) in splits(p, args.min_families, divergent, args.mode).items():
                res[ds][key].append(v)
                plc[ds][key].append(place)
    LAB = dict(zip(dsnames, args.labels)) if args.labels and len(args.labels) == len(dsnames) \
        else {d: d for d in dsnames}

    def passes(ds, key):
        v = res[ds].get(key, [])
        return len(v) >= args.min_trees and st.median(v) > args.threshold

    def support(key):
        """best median sCF over the datasets where the split passes; ties -> more trees"""
        return max((st.median(res[ds][key]), len(res[ds][key]))
                   for ds in dsnames if passes(ds, key))

    def placement(ds, key):
        c = Counter(plc[ds].get(key, []))
        if not c or not excl:
            return ""
        return "; ".join("%s (%d trees)" % (p, n) for p, n in c.most_common())

    allkeys = set(k for ds in dsnames for k in res[ds])
    # most supported first
    keep = sorted((k for k in allkeys if any(passes(ds, k) for ds in dsnames)),
                  key=lambda k: (-support(k)[0], -support(k)[1], k.lower()))

    # ------------------------------------------------------------ long TSV
    with open(args.out + '_long.tsv', 'w', encoding='utf-8', newline='') as fh:
        fh.write(TAB.join(["side_A", "side_B", "n_families_A", "dataset", "label",
                           "recovered", "n_trees",
                           "median_sCF", "q1_sCF", "q3_sCF", "mean_sCF",
                           "pct_trees_above_threshold", "passes", "divergent_placement"]) + NL)
        for k in keep:
            for ds in dsnames:
                v = res[ds].get(k, [])
                if v:
                    q = st.quantiles(v, n=4) if len(v) > 3 else [min(v), st.median(v), max(v)]
                    above = 100.0 * sum(1 for x in v if x > args.threshold) / len(v)
                    fh.write(TAB.join(k.split(" | ") + [str(len(sides(k)[0])), ds, LAB[ds],
                                       str(len(v)),
                                       str(ntrees[ds]), "%.2f" % st.median(v), "%.2f" % q[0],
                                       "%.2f" % q[2], "%.2f" % st.mean(v), "%.0f" % above,
                                       "yes" if passes(ds, k) else "no",
                                       placement(ds, k)]) + NL)
                else:
                    fh.write(TAB.join(k.split(" | ") + [str(len(sides(k)[0])), ds, LAB[ds], "0",
                                       str(ntrees[ds])] + [""] * 5 + ["no", ""]) + NL)

    # ------------------------------------------------------------ wide TSV
    with open(args.out + '_wide.tsv', 'w', encoding='utf-8', newline='') as fh:
        hdr = ["side_A", "side_B", "n_families_A", "best_median_sCF"]
        for ds in dsnames:
            hdr += [LAB[ds] + " recovered", LAB[ds] + " median_sCF"]
        fh.write(TAB.join(hdr) + NL)
        for k in keep:
            row = k.split(" | ") + [str(len(sides(k)[0])), "%.2f" % support(k)[0]]
            for ds in dsnames:
                v = res[ds].get(k, [])
                row += [str(len(v)), "%.2f" % st.median(v) if v else ""]
            fh.write(TAB.join(row) + NL)
    print("\nwrote " + args.out + "_long.tsv and " + args.out + "_wide.tsv")

    how = ("divergent families %s left out of the split names" % ",".join(sorted(divergent))
           if excl else "all families counted, divergent ones included")
    print("\n%d splits (named group of >= %d families; %s) recovered in >= %d trees with "
          "median sCF > %.0f in at least one dataset,\nsorted by their best median sCF (most "
          "supported first). Each split is printed as side A, with side B below it.\n"
          "Cells: '<trees> tr, <median sCF>%%' (* = passes)\n"
          % (len(keep), args.min_families, how, args.min_trees, args.threshold))
    for j, ds in enumerate(dsnames):
        print("  [%d] %s" % (j + 1, LAB[ds]))
    print("\n%-60s" % "split" + "".join("%17s" % ("[%d]" % (j + 1)) for j in range(len(dsnames))))
    for k in keep:
        cells = []
        for ds in dsnames:
            v = res[ds].get(k, [])
            cells.append("%17s" % (("%d tr, %.1f%%%s" % (len(v), st.median(v), "*" if passes(ds, k) else ""))
                                   if v else "-"))
        a, b = k.split(" | ")
        print("%-60s" % a[:60] + "".join(cells))
        print("    | " + b)
    print("\nNOTE  a split is an exact grouping of families. Nested splits (e.g. the same")
    print("      lineage with and without one more family) are listed separately; read them")
    print("      together. Per tree the best-sCF branch of a split is used, which is biased")
    print("      upward equally in every dataset.")
    if excl:
        print("      Where the divergent families attach is in the divergent_placement column")
        print("      of the _long.tsv.")

    # ---------------------------------------------------------------- figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap, Normalize
    except ImportError:
        print("\n(matplotlib unavailable - figure skipped)")
        return
    # figure only: splits recovered in >= plot_min_share of the trees in some dataset
    keep = [k for k in keep
            if any(len(res[ds].get(k, [])) >= args.plot_min_share * ntrees[ds] for ds in dsnames)]
    if not keep:
        return
    INK, INK2, GRID = "#0b0b0b", "#52514e", "#dcdcd8"
    # sequential blue ramp (steps 150 -> 700), one hue light -> dark
    cmap = LinearSegmentedColormap.from_list(
        "scf", ["#b7d3f6", "#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281", "#0d366b"])
    vals = [st.median(res[ds][k]) for k in keep for ds in dsnames if res[ds].get(k)]
    norm = Normalize(vmin=args.threshold, vmax=max(max(vals), args.threshold + 10))

    nr, nc = len(keep), len(dsnames)
    present = set(f for k in keep for sd in sides(k) for f in sd)
    fams = [f for f in FAM_ORDER if f in present] + sorted(present - set(FAM_ORDER))
    nf = len(fams)
    fig, (axm, ax) = plt.subplots(
        1, 2, sharey=True, figsize=(0.32 * nf + 2.3 * nc + 3.0, 2.4 + 0.36 * nr),
        gridspec_kw=dict(width_ratios=[0.32 * nf, 2.3 * nc], wspace=0.04))
    fig.patch.set_facecolor("white")

    def membership(axm):
        """grid of families: filled = side A (fewer families), hollow = side B"""
        for i, k in enumerate(keep):
            y = nr - 1 - i
            a, b = sides(k)
            for x, f in enumerate(fams):
                if f in a:
                    axm.scatter(x, y, s=46, color=INK, zorder=3)
                elif f in b:
                    axm.scatter(x, y, s=46, facecolors="white", edgecolors=INK2,
                                linewidths=1.0, zorder=3)
            xa = [x for x, f in enumerate(fams) if f in a]
            axm.plot([min(xa), max(xa)], [y, y], color=INK, lw=1.2, zorder=2)
        axm.set_xticks(range(nf))
        axm.set_xticklabels([DISPLAY.get(f, f) for f in fams], rotation=90, fontsize=8.5,
                            color=INK)
        axm.xaxis.tick_top()
        axm.set_xlim(-0.6, nf - 0.4)
        axm.set_ylim(-0.7, nr - 0.3)
        axm.set_yticks(range(nr))
        axm.set_yticklabels([])
        for sp in axm.spines.values():
            sp.set_visible(False)
        axm.tick_params(length=0)
        axm.yaxis.grid(True, color=GRID, lw=0.6)
        axm.set_axisbelow(True)
        axm.legend(handles=[plt.scatter([], [], s=46, color=INK, label="side A"),
                            plt.scatter([], [], s=46, facecolors="white", edgecolors=INK2,
                                        linewidths=1.0, label="side B")],
                   title="Split: which side each family is on", title_fontsize=8.5,
                   loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=2, frameon=False,
                   fontsize=8.5, labelcolor=INK)

    membership(axm)
    for i, k in enumerate(keep):
        y = nr - 1 - i
        for j, ds in enumerate(dsnames):
            v = res[ds].get(k, [])
            if not v:
                ax.text(j, y, "-", ha="center", va="center", fontsize=8, color=INK2)
                continue
            med = st.median(v)
            size = 20 + 380 * len(v) / ntrees[ds]
            if med > args.threshold:
                ax.scatter(j, y, s=size, color=cmap(norm(med)), edgecolors="white",
                           linewidths=1.5, zorder=3)
            else:
                ax.scatter(j, y, s=size, facecolors="none", edgecolors="#8a8a84",
                           linewidths=1.2, zorder=3)
            ax.text(j + 0.30, y, "%d trees, sCF %.0f%%" % (len(v), med), ha="left", va="center",
                    fontsize=7.5, color=INK2)
    ax.tick_params(axis="y", labelleft=False)
    ax.set_xticks(range(nc))
    ax.set_xticklabels([LAB[ds].replace(", ", "\n") for ds in dsnames], fontsize=9, color=INK)
    ax.xaxis.tick_top()
    ax.set_xlim(-0.6, nc - 0.1)
    ax.set_ylim(-0.7, nr - 0.3)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)
    ax.yaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=ax, fraction=0.025, pad=0.02)
    cb.set_label("median sCF (%)", fontsize=9, color=INK)
    cb.outline.set_visible(False)
    cb.ax.tick_params(labelsize=8, colors=INK2, length=0)

    # legend: dot area = how many of the trees recover the split
    ntr = max(ntrees.values())
    handles = [plt.scatter([], [], s=20 + 380 * n / 100.0, color="#5598e7",
                           edgecolors="white", linewidths=1.5,
                           label="%d%% of trees" % n) for n in (10, 50, 100)]
    handles.append(plt.scatter([], [], s=20 + 380 * 0.5, facecolors="none",
                               edgecolors="#8a8a84", linewidths=1.2,
                               label="median sCF at or below %.0f%%" % args.threshold))
    ax.legend(handles=handles, title="Dot area = share of the %d ML trees recovering the split"
              % ntr, title_fontsize=8.5, loc="upper center", bbox_to_anchor=(0.5, -0.01),
              ncol=4, frameon=False, fontsize=8.5, labelcolor=INK, borderpad=0.2,
              handletextpad=0.4, columnspacing=1.8, scatterpoints=1)

    dname = ", ".join(DISPLAY.get(f, f) for f in sorted(divergent))
    title = ("Family splits across 100 ML trees per dataset - divergent families (%s) excluded"
             % dname if excl else
             "Family splits across 100 ML trees per dataset - all families included")
    fig.suptitle(title, fontsize=12.5, color=INK, x=0.02, ha="left")
    note = ("Split names ignore %s, so a split is the same whether or not they sit inside it."
            % dname if excl else "Split names count every family, divergent ones included.")
    fig.text(0.02, 0.005, "Each row is one split (branch): side A | side B. %s\nShown: splits "
             "recovered in at least %.0f%% of the trees in at least one dataset, ordered by best "
             "median sCF, most supported first. Labels: trees recovering the split, and median "
             "sCF." % (note, 100 * args.plot_min_share), fontsize=8, color=INK2, ha="left")
    fig.tight_layout(rect=[0, 0.09, 1, 0.97])
    for ext in ("svg", "png"):
        fig.savefig(args.out + "." + ext, dpi=200, facecolor="white", bbox_inches="tight")
    print("\nwrote " + args.out + ".svg / .png")

    # ------------------------------------------------ box plots of the same splits
    boxplot_figure(keep, fams, res, dsnames, ntrees, LAB, args.threshold,
                   title.replace("Family splits", "sCF of family splits"), note,
                   args.out + "_boxplot")


def split_label(key, fams):
    """side A as text, one family per line; families on neither side are named below it"""
    a, b = sides(key)
    order = {f: i for i, f in enumerate(FAM_ORDER)}
    lab = "\n+ ".join(DISPLAY.get(f, f) for f in sorted(a, key=lambda f: (order.get(f, 99), f)))
    split = [DISPLAY.get(f, f) for f in fams if f not in a and f not in b]
    if split:
        lab += "\n(split: " + ", ".join(split) + ")"
    return lab


def boxplot_figure(keep, fams, res, dsnames, ntrees, LAB, threshold, title, note, out):
    """one group of boxes per split, named as text, one box per
    dataset built from the trees that recover the split; tree counts under the boxes"""
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import random
    INK, INK2, GRID = "#0b0b0b", "#52514e", "#dcdcd8"
    # validated categorical slots 1-5 (fixed order, one per dataset)
    COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
    n = len(dsnames)
    if n > len(COLORS):
        print("\n(more than %d datasets - box plot skipped)" % len(COLORS))
        return
    w = 0.8 / n
    offs = [(j - (n - 1) / 2) * w for j in range(n)]
    fig, ax = plt.subplots(figsize=(max(11.0, 0.55 * n * len(keep)), 7.0))
    fig.patch.set_facecolor("white")
    pos, series, cols = [], [], []
    for i, k in enumerate(keep):
        for j, ds in enumerate(dsnames):
            v = res[ds].get(k, [])
            if v:
                pos.append(i + offs[j])
                series.append(v)
                cols.append(COLORS[j])
    bp = ax.boxplot(series, positions=pos, widths=w * 0.8, patch_artist=True, showfliers=False,
                    medianprops=dict(color="white", lw=2.0),
                    whiskerprops=dict(color=INK2, lw=1.2),
                    capprops=dict(color=INK2, lw=1.2), boxprops=dict(lw=0))
    for patch, c in zip(bp['boxes'], cols):
        patch.set_facecolor(c)
        patch.set_edgecolor("white")
        patch.set_linewidth(2)
    random.seed(7)
    ymin = min(min(min(s) for s in series), threshold)
    for i, k in enumerate(keep):
        for j, ds in enumerate(dsnames):
            v = res[ds].get(k, [])
            ax.scatter([i + offs[j] + random.uniform(-w * 0.25, w * 0.25) for _ in v], v,
                       s=8, color=COLORS[j], alpha=0.35, linewidths=0, zorder=3)
            ax.text(i + offs[j], ymin - 4, "%d" % len(v), ha="center", fontsize=7.5,
                    color=INK2 if v else "#b3261e", fontweight="normal" if v else "bold")
    ax.axhline(threshold, ls=(0, (5, 4)), lw=1.4, color="#8a8a84", zorder=0)
    ax.text(-0.52, threshold - 1, "noise floor %.0f%%" % threshold, va="top", ha="left",
            fontsize=9.5, color=INK2)
    ax.set_xticks(range(len(keep)))
    labels = [split_label(k, fams) for k in keep]
    ax.set_xticklabels(labels, fontsize=9.5, color=INK)
    nlines = max(lab.count("\n") + 1 for lab in labels)
    ax.set_xlim(-0.55, len(keep) - 0.45)
    ax.set_ylim(ymin - 6, None)
    ax.set_ylabel("sCF of the split in each tree that recovers it (%)", fontsize=11, color=INK)
    ax.yaxis.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(GRID)
    ax.tick_params(axis="y", colors=INK2, length=0, labelsize=10)
    ax.tick_params(axis="x", length=0)
    ax.legend(handles=[mpatches.Patch(color=COLORS[j], label="%s  (n=%d trees)"
                                      % (LAB[ds], ntrees[ds]))
                       for j, ds in enumerate(dsnames)],
              loc="upper center", bbox_to_anchor=(0.5, -0.14 - 0.05 * nlines), ncol=min(n, 3),
              frameon=False, fontsize=10, labelcolor=INK)
    fig.suptitle(title, fontsize=13.5, color=INK, x=0.045, ha="left", y=0.985)
    fig.text(0.045, 0.935, "Each group is one split, named by its smaller side (those families "
             "together, against all the others). Each point is one ML tree that recovers the "
             "split.\nSmall numbers give how many trees recover it; a red 0 means that dataset "
             "never does. " + note, fontsize=9.5, color=INK2, ha="left", va="top")
    fig.tight_layout(rect=[0, 0.02, 1, 0.87])
    for ext in ("svg", "png"):
        fig.savefig(out + "." + ext, dpi=200, facecolor="white", bbox_inches="tight")
    print("wrote " + out + ".svg / .png")

if __name__ == '__main__':
    main()
