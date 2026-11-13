"""Supplementary Table S1: how each family was searched and how many sequences survived each stage.

Stage sources follow the latest file per stage:
  clade      step2_deeptmhmm/accessions_by_family_prefilter_Sep16.tsv
  7TM        step3 topology_msa_gt0.1/<family>_genes_7TMfilter_*.fasta (Sep22, incl. shard reruns),
             with the step1 nested-model exclusions applied (Class A - GPR1, Vomeronasal 1; Class C - Vomeronasal 2)
  topology   step3 topology_filter/accessions_by_family_7TM_topology0.8_Sep24.tsv
  final      step5 final_families/accessions_by_family_final_Sep24.tsv
Search setup is read from GPCR_Superfamily/code_and_data/homolog_search_pipeline (HMM and seed files)
and step1 identify_gene_clade.py (search used, tree variant). Sequences retrieved = top-ranked hits kept for the gene tree (Zenodo headers file).
The hit limits in today's hmm_search.sh were raised after some searches ran (cAMP: 10,000 retrieved), so they are not reported.

    python build_S1_search_and_filtering.py
"""
import ast, csv, gzip, re
from pathlib import Path

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from gpcr_paths import stage  # noqa: E402
CODE = Path(__file__).resolve().parent
HERE = stage("supplementary_tables")
HERE.mkdir(parents=True, exist_ok=True)
SEARCH = CODE.parent / "01_homology_search"
ZEN = stage("01_homology_search")

FAMILIES = ["classA", "Olfactory", "Vomeronasal1", "classT", "Nematode", "cAMP", "classB1", "classB2", "Mth",
            "classF", "FSLB", "GP143", "TM116", "GP157", "GPR1", "STE3", "classC", "Vomeronasal2", "STE2"]
CONTROLS = ["Galpha", "Gbeta", "Ggamma", "TBP", "RBP1"]
NICE = {"classA": "Class A", "classB1": "Class B1", "classB2": "Class B2", "classC": "Class C", "classF": "Class F",
        "classT": "Taste-2", "FSLB": "FSL", "Vomeronasal1": "Vomeronasal 1", "Vomeronasal2": "Vomeronasal 2",
        "Galpha": "G alpha", "Gbeta": "G beta", "Ggamma": "G gamma", "RBP1": "RPB1"}
ENTRY = {"GPR1": "GPR1_YEAST", "FSLB": "FSLB_DICDI", "classB1": "classB"}  # clade-extraction entry names
EXCLUDE = {"classA": ["GPR1", "Vomeronasal1"], "classC": ["Vomeronasal2"]}
TMFILE = "_genes_7TMfilter_db_FAMSA_topology_converted_gt0.1.fasta"


def columns(path):
    rows = list(csv.reader(open(path, encoding="utf-8"), delimiter="\t"))
    cols = {c: set() for c in rows[0]}
    for r in rows[1:]:
        for c, v in zip(rows[0], r):
            if v.strip():
                cols[c].add(v.strip())
    return cols


def acc(header):
    h = header[1:].split()[0].split("/")[0]
    return h.split("|")[1] if "|" in h else h


def fasta_ids(path):
    with open(path, encoding="utf-8") as f:
        return {acc(line) for line in f if line.startswith(">")}


def gene_boundaries():
    tree = ast.parse((CODE.parent / "02_clade_extraction" / "identify_gene_clade.py").read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.Assign) and getattr(node.targets[0], "id", "") == "gene_boundaries":
            return ast.literal_eval(node.value)


def query(search):
    d = SEARCH / "input_sequences" / search
    hmm = list(d.glob("*.hmm")) if d.exists() else []
    if hmm:
        t = hmm[0].read_text(encoding="utf-8", errors="ignore")
        name = re.search(r"^NAME\s+(\S+)", t, re.M).group(1)
        a = re.search(r"^ACC\s+(\S+)", t, re.M)
        return "HMMER profile", (f"{name} ({a.group(1).split('.')[0]})" if a else f"KEGG {name}")
    fa = list(d.glob("*.fasta")) if d.exists() else []
    if fa:
        h = open(fa[0], encoding="utf-8").readline()
        return "Iterative MMseqs2 profile", f"{h.split('|')[1]} ({h.split('|')[2].split()[0]})"
    # profile not staged in the repository: read it from the HMMER output on Zenodo
    with gzip.open(ZEN / search / f"{search}_hmmsearch_E1000.txt.gz", "rt") as f:
        head = f.read(4000)
    name = re.search(r"^Query:\s+(\S+)", head, re.M).group(1)
    a = re.search(r"^Accession:\s+(\S+)", head, re.M)
    return "HMMER profile", f"{name} ({a.group(1).split('.')[0]})" if a else name


def hits_retrieved(search):
    d = ZEN / search
    f = next(iter(sorted(d.glob("*headers*.txt.gz"))), None)
    if f is None:
        return ""
    with gzip.open(f, "rt") as h:
        return sum(1 for line in h if line.strip())


def main():
    gb = gene_boundaries()
    clade = columns(stage("02_clade_extraction", "clade_membership.tsv"))
    tm = columns(stage("03_topology_filter", "core_membership.tsv"))
    topo = columns(stage("03_topology_filter", "topology_membership.tsv"))
    final = columns(stage("04_overlap_assignment", "final_membership.tsv"))

    used_by = {}
    for k in FAMILIES + CONTROLS:
        e = gb[ENTRY.get(k, k)]
        used_by.setdefault(e.get("search_name", ENTRY.get(k, k)), []).append(NICE.get(k, k))
    rows = []
    for k in FAMILIES + CONTROLS:
        e = gb[ENTRY.get(k, k)]
        search = e.get("search_name", ENTRY.get(k, k))
        method, q = query(search)
        hr = e.get("hmmregion")
        tree = {1: "HMM-matched region", 0: "full-length"}.get(hr, "full-length")
        fam = k in FAMILIES
        rows.append({
            "Family": NICE.get(k, k),
            "Type": "GPCR family" if fam else "control",
            "Search method": method,
            "Query": q,
            "Search shared with": ", ".join(x for x in used_by[search] if x != NICE.get(k, k)),
            "Sequences retrieved": hits_retrieved(search),
            "Gene tree built from": tree,
            # clades added to the family; a third element -1 marks a clade subtracted from it
            "Extracted clades": sum(1 for b in e["boundary"] if len(b) < 3 or b[2] != -1),
            "Clade members": len(clade.get(k, ())) if fam else "",
            "Single 7TM core": len(tm[k]) if fam else "",
            "Topology score >= 0.8": len(topo.get(k, ())) if fam else "",
            "Final (after overlap assignment)": len(final.get(k, ())) if fam else "",
        })
    union = lambda d: len(set().union(*[d[k] for k in FAMILIES if k in d]))
    rows.append({"Family": "Total (unique sequences)", "Clade members": union(clade), "Single 7TM core": union(tm),
                 "Topology score >= 0.8": union(topo), "Final (after overlap assignment)": union(final)})

    out = HERE / "Table_S1_search_and_filtering.tsv"
    with open(out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print("wrote", out)


if __name__ == "__main__":
    main()
