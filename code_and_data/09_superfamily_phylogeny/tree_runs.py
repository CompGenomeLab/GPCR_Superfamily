"""Summary of the eight maximum likelihood runs (Table S6a), read from their RAxML-NG logs.

Per run: sequences, alignment columns after kpic-gappy trimming, the model chosen by ModelFinder, the
log-likelihood of the best tree, the number of ML searches and of bootstrap trees (autoMRE).

Reads  09_superfamily_phylogeny/runs/<set>/<set>_<alignment>_kpicgappy.raxml.log
       08_superfamily_alignment/alignments/<set>_<alignment>_kpicgappy.fasta
Writes 09_superfamily_phylogeny/tree_support_tables/S_tree_runs.tsv

Usage:
    python tree_runs.py
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402

SETS = [("All 19 families", "all"), ("All but STE2", "no_ste2"),
        ("Core (no STE2, STE3, Class C, V2)", "no_ste2ste3classCV2"), ("Core + 18 literature", "and_literature")]
ALIGNMENTS = [("G-INS-i --ep 0.123 (main)", "ginsi_ep0.123"), ("G-INS-i default --ep", "ginsi")]
HEADER = ["Representative set", "Alignment", "Sequences", "Alignment columns (kpic-gappy)",
          "Model (ModelFinder, BIC)", "Best-tree log-likelihood", "ML searches", "Bootstrap trees (autoMRE)"]


def main():
    rows = []
    for label, s in SETS:
        for alabel, a in ALIGNMENTS:
            name = f"representatives_{s}_{a}_kpicgappy"
            log = stage("09_superfamily_phylogeny", "runs", f"representatives_{s}", f"{name}.raxml.log").read_text()
            fasta = stage("08_superfamily_alignment", "alignments", f"{name}.fasta").read_text()
            nseq = sum(1 for line in fasta.splitlines() if line.startswith(">"))
            sites = re.search(r"Alignment sites / patterns: (\d+)", log).group(1)
            model = re.search(r"--model (\S+)", log).group(1)
            ll = re.search(r"Final LogLikelihood: (\S+)", log).group(1)
            bs = re.findall(r"Bootstrap tree #(\d+)", log)
            conv = re.search(r"Bootstopping test converged after (\d+) trees", log)
            nbs = conv.group(1) if conv else (max(map(int, bs)) if bs else "")
            rows.append([label, alabel, nseq, sites, model, f"{float(ll):.1f}", 100, nbs])
    out = stage("09_superfamily_phylogeny", "tree_support_tables", "S_tree_runs.tsv")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("\t".join(HEADER) + "\n")
        fh.writelines("\t".join(map(str, r)) + "\n" for r in rows)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
