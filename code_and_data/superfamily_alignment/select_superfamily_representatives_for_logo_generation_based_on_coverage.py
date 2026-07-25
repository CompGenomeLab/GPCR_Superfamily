
import os
import pandas as pd
from collections import Counter
from typing import List, Tuple, Optional, Dict

# ---- MSA I/O ----
def read_fasta_alignment(path):
    names, seqs = [], []
    name, chunks = None, []
    def flush():
        nonlocal name, chunks
        if name is not None:
            seq = "".join(chunks).strip().upper()
            if seq:
                names.append(name); seqs.append(seq)
        name, chunks = None, []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                flush(); name = line[1:].strip()
            else:
                chunks.append(line)
        flush()
    if not seqs: raise ValueError("No sequences found")
    if len({len(s) for s in seqs}) != 1: raise ValueError("Sequences must be same length")
    return names, seqs

# ---- Coverage (signed, independent of seq length) ----
def compute_signed_coverage(seqs, gap_chars="-.", weight_power=1.0, normalize_to_unit=False):
    n, L = len(seqs), len(seqs[0])
    gap = set(gap_chars)
    w = []
    for j in range(L):
        c = sum(seqs[i][j] not in gap for i in range(n))
        w.append((c / n) ** weight_power)
    scores = []
    for i in range(n):
        s = 0.0
        si = seqs[i]
        for j in range(L):
            s += w[j] if si[j] not in gap else -w[j]
        scores.append(s)
    if normalize_to_unit:
        denom = sum(w) or 1.0
        scores = [s / denom for s in scores]  # in [-1, 1]
    return scores

# ---- Accessionâ†’Family mapping (from your TSV) ----
def load_accession_family_tsv(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t", dtype=str, keep_default_na=False)
    # try common casings
    cols = {c.lower(): c for c in df.columns}
    acc_col = cols.get("accession", "accession" if "accession" in df.columns else None)
    fam_col = cols.get("family", "family" if "family" in df.columns else None)
    if not acc_col or not fam_col:
        raise ValueError(f"TSV must contain 'family' and 'accession'; got {list(df.columns)}")
    out = df[[fam_col, acc_col]].rename(columns={fam_col: "family", acc_col: "accession"})
    out["family"] = out["family"].astype(str).str.strip()
    out["accession"] = out["accession"].astype(str).str.strip()
    out = out[(out["family"] != "") & (out["accession"] != "")]
    return out.drop_duplicates()

# ---- OPTIONAL: make accessionâ†’family mapping unique to avoid fan-out (NEW) ----
def load_accession_family_tsv_unique(tsv_path, resolve="first"):
    """
    resolve: 'first' | 'mode'
    - 'first': keep first family per accession
    - 'mode' : keep the modal family per accession (ties -> first)
    """
    base = load_accession_family_tsv(tsv_path)
    if resolve == "first":
        return base.drop_duplicates(subset=["accession"], keep="first")
    elif resolve == "mode":
        return (base.groupby("accession", as_index=False)["family"]
                    .agg(lambda s: s.mode().iat[0] if not s.mode().empty else s.iat[0]))
    else:
        raise ValueError("resolve must be 'first' or 'mode'")

# ---- How to extract accession from FASTA header (CHANGED: made robust) ----
def default_accession_from_header(header):
    """
    Rule:
      - If '|' present, take the SECOND token (common for sp|ACC|... or tr|ACC|...)
      - Else take the FIRST whitespace-delimited token, then if it contains '|', take last,
        else if it contains '/', take last, else if it contains '_', take FIRST token.
    Examples:
      >sp|P12345|...     -> P12345
      >tr|Q8ABC1|blah    -> Q8ABC1
      >XP_0123 desc      -> XP_0123
      >gi|123|ref|NP_1|  -> NP_1
    """
    h = header.strip()
    tok = None
    if "|" in h:
        parts = h.split("|")
        tok = parts[1]
    elif "_" in h:
        parts=h.split("_")
        tok = parts[1]
    return tok

# ---- Master: compute + map â†’ family (UNCHANGED behavior, but allows unique mapping if desired) ----
def coverage_with_family_df(
    fasta_path,
    accession_tsv_path,
    gap_chars="-.",
    weight_power=1.0,
    normalize_to_unit=True,
    accession_extractor=default_accession_from_header,
    score_col_name=None,
    use_unique_mapping=False,  # NEW: set True to avoid fan-out before merge
    unique_resolve="first",    # NEW: 'first' or 'mode'
):
    names, seqs = read_fasta_alignment(fasta_path)
    scores = compute_signed_coverage(seqs, gap_chars=gap_chars, weight_power=weight_power,
                                     normalize_to_unit=normalize_to_unit)
    score_col = score_col_name or ("signed_coverage_norm" if normalize_to_unit else "signed_coverage")

    # base df (one row per input sequence, aligned to names/seqs order)
    df = pd.DataFrame({"sequence": names, score_col: scores})

    # map headers -> accession -> family
    acc_map = (load_accession_family_tsv_unique(accession_tsv_path, resolve=unique_resolve)
               if use_unique_mapping else
               load_accession_family_tsv(accession_tsv_path))

    df["accession"] = [accession_extractor(h) for h in df["sequence"]]
    df = df.merge(acc_map, on="accession", how="left")

    # return df and also the raw names/seqs so downstream funcs can align/merge by key
    return df[["sequence", "accession", "family", score_col]], names, seqs

# ---- Pick best representative per family ----
def best_per_family(df, score_col):
    g = df.dropna(subset=["family"]).copy()
    if g.empty:
        return g.head(0)
    idx = g.groupby("family")[score_col].idxmax()
    return g.loc[idx.values].sort_values("family").reset_index(drop=True)

# =========================
# NEW: Consensus + scoring
# =========================
def compute_consensus_and_scores(
    seqs: List[str],
    gap_chars: str = "-.X",
    min_column_non_gap_frac: float = 0.2,
    weights: Optional[List[float]] = None,   # per-column weights (optional)
    # scoring options
    return_identity: bool = True,
    return_simple_alignment_score: bool = True,
    match_score: float = 1.0,
    mismatch_score: float = -1.0,
    gap_penalty: float = -1.0,
    ignore_consensus_gaps: bool = True,
    substitution_matrix: Optional[Dict[Tuple[str, str], float]] = None,
) -> Tuple[str, Optional[List[float]], Optional[List[float]], Optional[List[float]]]:
    """
    Build a consensus sequence and score each sequence against it.

    Consensus:
      - For each column, most frequent non-gap residue if non-gap fraction >= threshold,
        else consensus is '-' at that column.
      - Ties broken lexicographically for determinism.

    Scores (per sequence):
      - identity fraction (0..1) vs consensus (skipping consensus-gap cols if ignore_consensus_gaps).
      - simple alignment average score with match/mismatch/gap penalties.
      - optional substitution-matrix average score if a dict is provided.
    """
    if not seqs:
        raise ValueError("No sequences provided")
    L = len(seqs[0])
    if any(len(s) != L for s in seqs):
        raise ValueError("Sequences must be same length")
    gap = set(gap_chars)
    n = len(seqs)

    if weights is None:
        weights = [1.0] * L
    if len(weights) != L:
        raise ValueError("weights must be length L")

    # consensus
    cons = []
    for j in range(L):
        col = [s[j] for s in seqs if s[j] not in gap]
        if len(col) / n < min_column_non_gap_frac:
            cons.append("-"); continue
        counts = Counter(col)
        m = max(counts.values())
        candidates = sorted([aa for aa, c in counts.items() if c == m])
        cons.append(candidates[0])
    consensus = "".join(cons)

    def iter_positions():
        for j in range(L):
            cj = consensus[j]
            if ignore_consensus_gaps and cj in gap:
                continue
            yield j, cj

    identity_frac = None
    if return_identity:
        idv = []
        for si in seqs:
            match = 0
            denom = 0
            for j, cj in iter_positions():
                a = si[j]
                if cj in gap:
                    if a in gap and not ignore_consensus_gaps:
                        match += 1; denom += 1
                    elif not ignore_consensus_gaps:
                        denom += 1
                else:
                    if a in gap:
                        denom += 1
                    else:
                        match += (a == cj)
                        denom += 1
            idv.append((match / denom) if denom else 0.0)
        identity_frac = idv

    simple_avg_score = None
    if return_simple_alignment_score:
        avgs = []
        for si in seqs:
            total = 0.0
            denom = 0
            for j, cj in iter_positions():
                a = si[j]
                if cj in gap or a in gap:
                    if cj in gap and a in gap:
                        if not ignore_consensus_gaps:
                            denom += 1
                        continue
                    total += gap_penalty
                    denom += 1
                else:
                    total += match_score if a == cj else mismatch_score
                    denom += 1
            avgs.append((total / denom) if denom else 0.0)
        simple_avg_score = avgs

    matrix_avg_score = None
    if substitution_matrix is not None:
        def smat(a: str, b: str) -> float:
            if a in gap or b in gap:
                return gap_penalty
            if (a, b) in substitution_matrix:
                return substitution_matrix[(a, b)]
            if (b, a) in substitution_matrix:
                return substitution_matrix[(b, a)]
            return mismatch_score

        mavgs = []
        for si in seqs:
            total = 0.0
            denom = 0
            for j, cj in iter_positions():
                a = si[j]
                if cj in gap and a in gap:
                    if not ignore_consensus_gaps:
                        denom += 1
                    continue
                total += smat(a, cj)
                denom += 1
            mavgs.append((total / denom) if denom else 0.0)
        matrix_avg_score = mavgs

    return consensus, identity_frac, simple_avg_score, matrix_avg_score

# ---- Merge-safe: add consensus scores into df keyed by 'sequence' or 'accession' (NEW) ----
def add_consensus_scores_to_df_merged(
    df: pd.DataFrame,
    names: List[str],
    seqs: List[str],
    key: str = "sequence",      # merge key: 'sequence' or 'accession'
    gap_chars: str = "-.X",
    min_column_non_gap_frac: float = 0.2,
    add_identity_col: bool = True,
    add_simple_score_col: bool = True,
    match_score: float = 1.0,
    mismatch_score: float = -1.0,
    gap_penalty: float = -1.0,
    ignore_consensus_gaps: bool = True,
    substitution_matrix: Optional[Dict[Tuple[str, str], float]] = None,
    identity_col_name: str = "consensus_identity",
    simple_col_name: str = "consensus_simple_score",
    matrix_col_name: str = "consensus_matrix_score",
) -> Tuple[pd.DataFrame, str]:
    """
    Computes consensus on `seqs` and merges per-sequence scores into `df` by `key`.
    This avoids length mismatches when `df` has duplicates due to earlier merges.
    Returns (new_df, consensus_string).
    """
    consensus, idf, simp, mat = compute_consensus_and_scores(
        seqs,
        gap_chars=gap_chars,
        min_column_non_gap_frac=min_column_non_gap_frac,
        weights=None,
        return_identity=add_identity_col,
        return_simple_alignment_score=add_simple_score_col,
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_penalty=gap_penalty,
        ignore_consensus_gaps=ignore_consensus_gaps,
        substitution_matrix=substitution_matrix,
    )

    if key == "sequence":
        keys = names
    elif key == "accession":
        # If merging by accession, ensure df has 'accession' and compute per-name accession list
        # from the same extractor used when df was created.
        if "accession" not in df.columns:
            raise KeyError("'accession' not present in df; cannot merge by accession.")
        # Rebuild the accession mapping from the alignment names.
        # (If you want this, pass in the extractor and recompute from names as needed.)
        # Use accessions already keyed by sequence when available.
        # Better: keep using 'sequence' unless you truly need 'accession'.
        keys = None
        raise NotImplementedError("Merging by 'accession' requires a names->accession mapping.")
    else:
        raise ValueError("key must be 'sequence' or 'accession'")

    scores_cols = {key: keys}
    if add_identity_col and idf is not None:
        scores_cols[identity_col_name] = idf
    if add_simple_score_col and simp is not None:
        scores_cols[simple_col_name] = simp
    if substitution_matrix is not None and mat is not None:
        scores_cols[matrix_col_name] = mat

    scores_df = pd.DataFrame(scores_cols)
    if key not in df.columns:
        raise KeyError(f"'{key}' column not found in df to merge scores on.")
    out = df.merge(scores_df, on=key, how="left")
    return out, consensus

# =========================
# Example usage (end-to-end)
# =========================
if __name__ == "__main__":
    STAGE_DIR = os.path.dirname(os.path.abspath(__file__))
    PROJECT_DIR = os.environ.get("GPCR_ANALYSIS_DIR", os.environ.get("GPCR_ZENODO_DIR", "/path/to/GPCR_Superfamily"))
    msa = os.path.join(STAGE_DIR, "data", "alignments", "All_clust0.7_cov1_Oct25_minsize2_linsi_ep0.123_trimends_treein_ginsi_ep0.123.fasta")
    acc_path = os.environ.get("GPCR_ACCESSION_TSV", os.path.join(PROJECT_DIR, "families", "tables", "acc_taxid_map_BestEval_May22_e0.1.tsv"))
    reps_out = os.path.join(STAGE_DIR, "data", "logo_outputs", "logo_representatives_May22.tsv")

    # Build coverage + mapping; set use_unique_mapping=True to avoid fan-out at source
    df, names, seqs = coverage_with_family_df(
        fasta_path=msa,
        accession_tsv_path=acc_path,
        gap_chars="-.X",
        weight_power=1.0,
        normalize_to_unit=True,
        use_unique_mapping=False,      # set True if you prefer de-duplicated mapping
        unique_resolve="first",
    )

    # Add consensus-based scores in a merge-safe way (no length mismatch even with duplicates)
    df2, consensus = add_consensus_scores_to_df_merged(
        df=df, names=names, seqs=seqs, key="sequence",
        gap_chars="-.X", min_column_non_gap_frac=0.2,
        add_identity_col=True, add_simple_score_col=True,
        match_score=1.0, mismatch_score=-1.0, gap_penalty=-1.0,
        ignore_consensus_gaps=True, substitution_matrix=None
    )

    # Representatives by consensus identity (or use 'consensus_simple_score')
    reps_by_identity = best_per_family(df2, score_col="consensus_identity")
    reps_by_identity.to_csv(reps_out, sep="\t", index=False)
    print("Consensus (first 120):", consensus[:120], "...")
    print(reps_by_identity)
    print(f"Wrote representatives table: {reps_out}")

#%%
import pandas as pd

def _read_fasta(path):
    names, seqs = [], []
    name, chunks = None, []
    def flush():
        nonlocal name, chunks
        if name is not None:
            seq = "".join(chunks).strip()
            if seq:
                names.append(name); seqs.append(seq)
        name, chunks = None, []
    with open(path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line: continue
            if line.startswith(">"):
                flush(); name = line[1:].strip()
            else:
                chunks.append(line.strip())
        flush()
    if not seqs: raise ValueError("No sequences found")
    if len({len(s) for s in seqs}) != 1: raise ValueError("Sequences must be same length")
    return names, seqs

def _write_fasta(path, names, seqs, width=60):
    with open(path, "w") as f:
        for h, s in zip(names, seqs):
            f.write(f">{h}\n")
            for i in range(0, len(s), width):
                f.write(s[i:i+width] + "\n")

def extract_representative_alignment(in_fasta, reps, out_fasta, gap_chars="-.", keep_order=True):
    """
    reps: list of header strings, or a DataFrame containing column 'sequence' with headers.
    Produces a new FASTA with only the selected headers, after dropping columns that are all gaps across them.
    """
    all_names, all_seqs = _read_fasta(in_fasta)
    L = len(all_seqs[0]); gap = set(gap_chars)

    if isinstance(reps, pd.DataFrame):
        sel_names = reps["sequence"].astype(str).tolist()
    #     sel_names.append("sp_P28646_SSR1_RAT/60-346")
    #     sel_names.remove("tr_A0A0G4LZA0_A0A0G4LZA0_VERLO/721-1177")
    #     sel_names.append("tr|A0A0G4IQW2|A0A0G4IQW2_PLABS/463-739 G-protein coupled receptors family 3 profile domain-containing protein (Fragment) OS=Plasmodiophora brassicae OX=37360 GN=PBRA_000904 PE=4 SV=1")
    #     sel_names.remove("tr|A0A195FUD1|A0A195FUD1_9HYME/672-995 Gamma-aminobutyric acid type B receptor subunit 2 OS=Trachymyrmex septentrionalis OX=34720 GN=ALC56_01602 PE=3 SV=1")
    else:
        sel_names = list(reps)

        
    # map header -> sequence
    idx = {h: i for i, h in enumerate(all_names)}
    missing = [h for h in sel_names if h not in idx]
    if missing:
        raise ValueError(f"Selected headers not found in alignment: {missing[:5]}{' ...' if len(missing)>5 else ''}")

    # honor desired order or use original alignment order
    if keep_order:
        pick = sel_names
    else:
        # original order from the input alignment
        present = set(sel_names)
        pick = [h for h in all_names if h in present]

    sel_seqs = [all_seqs[idx[h]] for h in pick]

    # build keep-mask: drop columns that are all gaps among the selected reps
    keep = [any(s[j] not in gap for s in sel_seqs) for j in range(L)]
    if not any(keep):
        raise ValueError("All columns are gaps across selected representatives.")

    pruned = ["".join(s[j] for j in range(L) if keep[j]) for s in sel_seqs]

    _write_fasta(out_fasta, pick, pruned)
    return {"n_sequences": len(pick), "orig_len": L, "new_len": len(pruned[0]), "dropped_cols": keep.count(False)}

import pandas as pd

def _extract_acc_from_header(header: str) -> str:
    """
    Extract UniProt accession from FASTA header.
    - If '|' present, split and take 2nd field.
    - Else if '_' present, split and take 2nd field.
    - Else use full header.
    """
    if "|" in header:
        parts = header.split("|")
        if len(parts) > 1:
            return parts[1]
    if "_" in header:
        parts = header.split("_")
        if len(parts) > 1:
            return parts[1]
    return header

def extract_alignment_by_acc1(
    in_fasta,
    accessions,
    out_fasta,
    gap_chars="-.",
    keep_order=True,
):
    """
    Extract only sequences whose accessions (acc1 column) match headers in FASTA.

    Parameters
    ----------
    in_fasta : str
        Path to aligned FASTA file.
    accessions : list[str] or pandas.DataFrame
        List of accessions or DataFrame with 'acc1' column.
    out_fasta : str
        Path to write filtered alignment.
    gap_chars : str
        Characters considered as gaps.
    keep_order : bool
        If True, output sequences in the order of provided acc1 list.
    """
    all_names, all_seqs = _read_fasta(in_fasta)
    L = len(all_seqs[0])
    gap = set(gap_chars)

    # normalize accession list
    if isinstance(accessions, pd.DataFrame):
        if "acc1" not in accessions.columns:
            raise ValueError("DataFrame must have 'acc1' column")
        acc_list = accessions["acc1"].astype(str).tolist()
    else:
        acc_list = [str(a) for a in accessions]

    # build header -> accession map
    header_to_acc = {h: _extract_acc_from_header(h) for h in all_names}
    acc_to_header = {}
    for h, acc in header_to_acc.items():
        acc_to_header.setdefault(acc, []).append(h)

    # choose headers matching requested acc1
    missing = []
    pick_headers = []
    for acc in acc_list:
        hits = acc_to_header.get(acc, [])
        if not hits:
            missing.append(acc)
        else:
            pick_headers.append(hits[0])  # take first if multiple

    if not pick_headers:
        raise ValueError(f"No matches found. Missing: {missing}")

    # if keep_order = False, respect input FASTA order
    if not keep_order:
        pick_set = set(pick_headers)
        pick_headers = [h for h in all_names if h in pick_set]

    sel_seqs = [all_seqs[all_names.index(h)] for h in pick_headers]

    # prune columns that are all gaps
    keep_mask = [any(s[j] not in gap for s in sel_seqs) for j in range(L)]
    pruned = ["".join(s[j] for j in range(L) if keep_mask[j]) for s in sel_seqs]

    _write_fasta(out_fasta, pick_headers, pruned)

    return {
        "n_sequences": len(pick_headers),
        "orig_len": L,
        "new_len": len(pruned[0]),
        "dropped_cols": keep_mask.count(False),
        "missing_accessions": missing,
    }


def write_logo_representatives_from_accessions(
    accessions,
    accession_tsv_path,
    out_tsv_path,
):
    """
    Build a logo-representatives table with columns required by prepare_trim_info.py:
      - accession
      - family

    If multiple selected accessions map to the same family, keeps the first one.
    """
    if isinstance(accessions, pd.DataFrame):
        if "acc1" not in accessions.columns:
            raise ValueError("DataFrame must contain 'acc1' column")
        acc_list = accessions["acc1"].astype(str).str.strip().tolist()
    else:
        acc_list = [str(a).strip() for a in accessions]

    acc_list = [a for a in acc_list if a]
    if not acc_list:
        raise ValueError("No valid accessions provided.")

    selected = pd.DataFrame({"accession": acc_list})
    acc_family = load_accession_family_tsv(accession_tsv_path)
    reps = selected.merge(acc_family, on="accession", how="left")

    missing_family = reps["family"].isna() | (reps["family"].astype(str).str.strip() == "")
    missing_accessions = reps.loc[missing_family, "accession"].tolist()
    reps = reps.loc[~missing_family, ["accession", "family"]].copy()

    # prepare_trim_info.py expects one representative per family.
    reps = reps.drop_duplicates(subset=["family"], keep="first").reset_index(drop=True)
    reps.to_csv(out_tsv_path, sep="\t", index=False)

    return {
        "out_tsv": out_tsv_path,
        "n_rows_written": int(len(reps)),
        "missing_family_for_accessions": missing_accessions,
    }
STAGE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.environ.get("GPCR_ANALYSIS_DIR", os.environ.get("GPCR_ZENODO_DIR", "/path/to/GPCR_Superfamily"))
msa = os.path.join(STAGE_DIR, "data", "alignments", "All_clust0.7_cov1_Oct25_minsize2_linsi_ep0.123_trimends_treein_ginsi_ep0.123.fasta")
out_path = os.path.join(STAGE_DIR, "data", "alignments", "All_clust0.7_cov1_Oct25_minsize2_linsi_ep0.123_trimends_treein_ginsi_ep0.123_selected_reps.fasta")
acc_path = os.environ.get("GPCR_ACCESSION_TSV", os.path.join(PROJECT_DIR, "families", "tables", "acc_taxid_map_BestEval_May22_e0.1.tsv"))
logo_reps_out = os.path.join(STAGE_DIR, "data", "logo_outputs", "logo_representatives_May22_custom.tsv")
custom_accessions = ["A0A835NM65","A0A3M2SGC9","A0A8J4PRC8","A0A6V7VX31","A0A8C0QWZ4",
                     "A0A395RUA3","A0A2I4B258","A0A7N6AN26","A0A7J6UXT4","A0A556TZG5",
                     "A0A6J2X0B2","A0A0M9A0X0","A0A8C7JL03","A0AA97IXW4","Q646E2",
                     "A0A6P7ZWE1","A0A9P6IC23"]

stats = extract_alignment_by_acc1(
    in_fasta=msa,
    accessions=custom_accessions,
    out_fasta=out_path,
    keep_order=True
)
print(stats)

trim_ready = write_logo_representatives_from_accessions(
    accessions=custom_accessions,
    accession_tsv_path=acc_path,
    out_tsv_path=logo_reps_out,
)
print(trim_ready)
#%%
rep_msa = os.path.join(STAGE_DIR, "data", "alignments", "All_clust0.7_cov1_Oct25_minsize2_linsi_ep0.123_trimends_treein_ginsi_ep0.123_familyreps_only.fasta")

stats = extract_representative_alignment(
    in_fasta=msa,
    reps=reps_by_identity,                 # or ['header1','header2',...]
    out_fasta=rep_msa,
    gap_chars="-.X",               # adjust if needed
    keep_order=True               # keep order from reps_df; set False to keep input alignment order
)
print(stats)  
