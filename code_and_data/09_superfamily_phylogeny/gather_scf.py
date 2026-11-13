#!/usr/bin/env python3
"""
gather_scf.py - collapse a directory of IQ-TREE *.cf.tree files into per-tree
site-concordance summaries, stratified by how many whole families each branch unites.

Usage:
    python gather_scf.py  <dir_of_cf_trees>  [<dir2> ...]   -o summary.tsv

A family counts toward a branch only if the branch's smaller side contains MORE THAN
50% of that family's leaves *present in that tree*. Strata reported per tree:
    all        - every internal branch
    fam2       - branches uniting >= 2 whole families
    fam3       - >= 3
    fam4       - >= 4

Self-contained: standard library only (no biopython). The accession -> family map is
embedded below, so this runs anywhere. It is
09_superfamily_phylogeny/representatives_family_map.tsv in the archive (the 590
representatives and the 18 literature sequences, LIT_<family>) with Vomeronasal2 merged
into classC; the other classA-like lineages stay separate so discover_groupings.py can
report splits among them.

sCF has a noise floor of ~33% (three possible quartet resolutions).
"""
import os, sys, glob, argparse, statistics as st

MAJORITY = 0.50          # fraction of a family's leaves required for it to count

FAMILY_ACCESSIONS = {
    "cAMP": "A0A061EQI4 A0A078HXL8 A0A087XSI0 A0A094FQE1 A0A0D2K054 A0A0D2X067 A0A0G2EUR8 A0A0L0DTH4 A0A139AZX0 A0A152A1B2 A0A178F304 A0A1I8JHP7 A0A1J6JVK8 A0A1R2CWP8 A0A1S3KHA5 A0A1X7UNV9 A0A1Y1V7I4 A0A1Y1XT38 A0A1Y2D3Q0 A0A1Y2HQ33 A0A210PRZ4 A0A2C9JF31 A0A2G9I309 A0A2K1JDD1 A0A2P6NBR2 A0A2R6Q399 A0A2R6WMW1 A0A2T3AEB5 A0A2T4ABJ6 A0A2U1LDU6 A0A2Z6LKF1 A0A369RZ92 A0A397IIR7 A0A401SEN3 A0A433PQ83 A0A485K453 A0A4P9XTV3 A0A4V3SHI5 A0A5A8CHT9 A0A5B8MC35 A0A5D6YHP3 A0A5J4YHZ4 A0A5J5B7Q9 A0A5N5D2A2 A0A6A1W3R1 A0A6A6ZRM0 A0A6J2WDI0 A0A6P4XUW4 A0A6P6U1M5 A0A6P7TB59 A0A6P7Y615 A0A7J6CPL2 A0A7J6UXT4 A0A7M5XC19 A0A7M7NJA1 A0A830HCH9 A0A835Z4Y8 A0A8B9AM48 A0A8E2JN74 A0A8H2ZSD0 A0A8J1UF86 A0A8J1W3G8 A0A8J2SMA0 A0A8J4U4L9 A0A8J5XEQ1 A0A8S1WDB8 A0A8S3YLQ6 A0A8T2RLE7 A0A913XCC3 A0A913ZHG2 A0A9D4G800 A0A9P3F8G8 A0A9P3I6N1 A0A9Q0GNG0 A0A9Q0HS40 A0A9Q1BRQ1 A0A9W7E5P5 A0A9W8AGM6 A0A9Y3R5U7 A0AA39V9N6 A0AA88DWJ3 A0AA88KQN4 A0AA88QNZ3 A0AAD1QZX7 A0AAD2CF14 A0AAD8XXJ3 A0AAE0BDP6 A0AAI9V691 A0AB34JYC7 A9V493 D8LJZ4 D8RY89 F4PSJ3 F6TIH6 G1XU22 I7M4K7 K8EP75 L8HF38 M2X2Z8 R7YN19",
    "classA": "A0A067MZ95 A0A067PRG3 A0A094IZH7 A0A0C3EGZ0 A0A0C9YUR1 A0A0K2TJR3 A0A0L0HU57 A0A0V0RWY4 A0A137NSW2 A0A139A5U2 A0A141CIT9 A0A164SPR8 A0A164VAP4 A0A165VGM7 A0A177B8L8 A0A1B7MRK5 A0A1D2MEE0 A0A1R2C0R6 A0A1S3HSH5 A0A1X6MRU2 A0A1X7VFF6 A0A1Y1WVZ6 A0A1Y2AK03 A0A1Y2FNF1 A0A1Y2I389 A0A2H3EBV1 A0A2P4XTB2 A0A2R5GL97 A0A2Z6RS11 A0A369RWL5 A0A397ULW8 A0A3M6UW13 A0A3P7L4P6 A0A3P8VRX9 A0A401P974 A0A423SZW9 A0A433PAR7 A0A4S4LF44 A0A4S8LVH9 A0A507C1F0 A0A507FNE7 A0A550CJ50 A0A5A8CIN3 A0A5A8CMQ7 A0A662WL24 A0A662WRK0 A0A670HQH4 A0A6A4WEN0 A0A6A6GWZ6 A0A6J3FV34 A0A6J8C7K3 A0A6P5A842 A0A7J7JRN6 A0A7K5KRP8 A0A7R9ABD7 A0A811L1U0 A0A814NHY2 A0A8B7DYR5 A0A8C4T3M9 A0A8H3VW39 A0A8K0UNT1 A0A8N4QEJ1 A0A8S0TKG9 A0A8S1QXI8 A0A8S1UAW2 A0A8T2KHA8 A0A914BDG9 A0A9N8E100 A0A9N8HBH2 A0A9P5YKF0 A0A9P7G7X8 A0A9P8ZU35 A0A9Q1GZJ3 A0A9W7A5F6 A0A9W7AAI4 A0A9W7BYL4 A0A9W7DUY4 A0A9W7EU40 A0A9W7G9J4 A0A9W7L9F1 A0A9W8IHB9 A0A9X6NBY2 A0AA38CIT1 A0AAD2FVI5 A0AAD3DE31 A0AAD4BX34 A0AAD4QLR9 A0AAD5Y0R0 A0AAE1ZBT5 A0AB34JPL6 A0C0C8 A2DIV7 D2VSQ8 F8PYW2 H2YV51 H3AAD0 M5GG22 R7TCL2 T1J5B1 V4ASX2",
    "classB1": "A0A1S3K882 A0A2G8KNL6 A0A2R2MNS9 A0A3M7RHS0 A0A4Z2DGP9 A0A7D9HW74 A0A7I8V7S9 A0A7J7KSJ8 A0A7J7KT91 A0A7K9BLV0 A0A7M5V2A6 A0A7M7M809 A0A8B7XGA8 A0A8S1CXC6 A0A8S3ZGZ9 A0A915IH15 A0A9J7LTB3 A0A9X6RLR2 B3RW40 T1G164",
    "classB2": "A0A0D2UBM4 A0A0D2VQ40 A0A164RH33 A0A1V9YY12 A0A369RQL8 A0A6P3HV15 A0A6P8HDQ9 A0A7J7KG21 A0A815CAI8 A0A8S4NTU8 A0A914BCJ3 A0A914QFA0 A0A9X6RN66 A0AA35W825 A0AAE1B7R8 A9V191 A9V5K9 F2UDI1 F2UP54 J0WSA6",
    "classC": "A0A075B0P8 A0A098VSP5 A0A0D2U3H4 A0A0G4IQW2 A0A0L0HT55 A0A0L0SUZ9 A0A138ZXC1 A0A1I8N4A0 A0A1S3IFX7 A0A1S3JWS9 A0A1S3JY51 A0A1X2GMP6 A0A1Y1HVE1 A0A1Y1VIX4 A0A1Y1VQX5 A0A210PUQ9 A0A2R5GA48 A0A2R6XQI1 A0A2T7P3C3 A0A369RWK2 A0A388M4N7 A0A3M7QXM5 A0A3P9PEG3 A0A433PK10 A0A4P9XVT1 A0A4W3HJ78 A0A5B8MTU1 A0A5J4NJY0 A0A7E6ELV4 A0A7J7JZY6 A0A7M5XAT4 A0A7M7NE72 A0A8B7XHQ2 A0A8B7Y2N7 A0A8B7ZG62 A0A8C4RT37 A0A8C5PIM3 A0A8C5Q574 A0A8H7U9G9 A0A8J4PN46 A0A8J4RH77 A0A8J6AVV1 A0A8S4PJ26 A0A8W8K266 A0A914A7G3 A0A914AAA3 A0A914UAW7 A0A9D4ZEZ9 A0A9J7N6I7 A0A9N8ER37 A0A9P5S8I3 A0A9Q0YGX5 A0A9W7C6X4 A0A9W7INE9 A0A9W8B043 A0AA38GB80 A0AA88GCU7 A0AA88XKZ9 A0AAD8Y2Z1 A0AAE0LFX7 A0AAE0S8M4 A0AAE1ED28 A0AAE1J9I8 A0AAN0J1S3 A0AB34IGP0 D8SMS6 H3A402 K3X562 V4MMV8 X6NCG3",
    "classF": "A0A0R3VUA5 A0A177B8U5 A0A1D1UVD4 A0A1S3HB76 A0A1S3KAW5 A0A2G5VVN5 A0A369S4J8 A0A6P6RAQ2 A0A7J7J8I8 A0A7J7KIZ4 A0A813WNK9 A0A8J1TCG0 A0A8J9YNW6 A0A9P0G5Q8 A0A9Q0YGH6 A0AA35TF90 A0AA88H9N0 A7S829 T1FSL1 V3ZM42",
    "classT": "A0A1S2ZP50 A0A444TXK0 A0A4X2LWR5 A0A5C6N243 A0A672F6C2 A0A6P7WQS4 A0A7K7I6N6 A0A7K9UT14 A0A8B7B095 A0A8C4F3Y3 A0A8C4S1N3 A0A8C5ZPJ5 A0A8C6LDE3 A0A8D2LF50 A0A8J1JLB4 A0A9Y3VHB0 A0AAD7RQY1 H3AD47 Q58WW7 W5LU70",
    "FSLB": "A0A0L0HC06 A0A0L0SPS2 A0A0L0SZS7 A0A0U1LJ34 A0A151Z3M2 A0A1U7LQE6 A0A1Y2BSR3 A0A1Y2HNZ1 A0A397S714 A0A4P9XGJ8 A0A4P9XNW4 A0A4S2MYS9 A0A9P6RSH0 A0A9P6XGJ3 A0AAD5X2I7 B2WDA8 F4PWJ4 L8GI96 L8HDT2 Q54DP2",
    "GP143": "A0A067QTP1 A0A0D2U8X7 A0A1S3J717 A0A1X7V8B4 A0A3R7MGS6 A0A401SQU1 A0A6P7XQ45 A0A812CE34 A0A8B7XY04 A0A8B8EXZ9 A0A8C7YJ70 A0A8J1UYI6 A0A8X6FNY3 A0A913X388 A0A914BG91 A0A9Q1HF54 A0A9W9Z0P0 A0AAD9PF82 A0AAD9Q9C8 A7RKB4 C3Z4T4 E9HUQ4 R7U5Y3 R7UDL5 V4A951",
    "GP157": "A0A0D2WJG7 A0A1S3H5X8 A0A1S3IU62 A0A2U9BZ00 A0A401PMB4 A0A669QYS9 A0A6P8ILI2 A0A6S7FXV2 A0A7M5X592 A0A8B8E8F0 A0A8J1UGB3 A0A8J9YN25 A0A9D4ENU2 A0A9Q0LWX6 A0A9W7BA54 A0A9W7BJY2 A0A9W9YSL0 A0AA35TKN5 A0AA35TKV4 A0AAD8F057 A0AAE0SEM9 A0AAN0IB11 A0AAN0IB64 L8GKI9 R7ULM5",
    "GPR1": "A0A0D1XY51 A0A164NX77 A0A1E3QHL3 A0A1L7WSD3 A0A1Y2FCU1 A0A1Y2G0K7 A0A284QQB0 A0A2P6NTM4 A0A395I952 A0A397JMS8 A0A433D253 A0A4S4L5U1 A0A507CBT8 A0A507D4A0 A0A8H4AHG8 A0A8H6KYP2 A0A9N8VA50 A0A9N8WR67 A0A9N9BSC1 A0A9P7RDW1 A0A9Q3EJW1 A0AAD3HAS6 A0AAD5U6C1 G7DUL4 R9ARF1",
    "LIT_classT": "LIT0360 LIT0361",
    "LIT_GP157": "LIT0015 LIT0022 LIT0024 LIT0025 LIT0063 LIT0081 LIT0114 LIT0119 LIT0161 LIT0162 LIT0171 LIT0197",
    "LIT_GPR1like": "LIT0202 LIT0203",
    "LIT_Vomeronasal1": "LIT0217 LIT0234",
    "Mth": "A0A0R3TMP9 A0A1D1VAJ6 A0A1S3HWV8 A0A1S3K1H8 A0A369S099 A0A6P7ST12 A0A7J7JVL8 A0A7J7KG56 A0A813M4Z3 A0A835CS02 A0A8B7DN18 A0A8J1UBI2 A0A8J9VLT2 A0A8X6VY91 A0A9J7MB36 A0A9Q1BS09 A0AA35T3X5 A7RJP1 D3B489 T1EJ83",
    "Nematode": "A0A0B2W259 A0A0N4VPR2 A0A0N4WPX6 A0A1I7S9W5 A0A454XQI6 A0A813U690 A0A814MAJ6 A0A8S9ZUJ4 A0A913IBI4 A0A914EF92 A0A914HSS7 A0A914K5P8 A0A914WSC1 A0A915I390 A0AA36M610 A0AA39IF33 A0AAD4MX89 A0AAD5N0E3 Q9N2T6 W2TZC1",
    "Olfactory": "A0A091CSE8 A0A1L8HZM7 A0A287ARC1 A0A2Y9E7Y1 A0A3Q2ZR88 A0A401NGY3 A0A401PQX7 A0A401T3Q8 A0A401T3V3 A0A452GV18 A0A498M660 A0A670JGL5 A0A7K5PGR8 A0A8C4NAE1 A0A8C4QJS3 A0A8X7WW43 A0AA88M8G2 A0AA88SB49 I3KY00 W5M3L1",
    "STE2": "A0A074VLE7 A0A0E9NKH5 A0A0G2HF89 A0A165HIN9 A0A1E3PZ10 A0A1E4TIF5 A0A1U7LRS2 A0A2T0FQ04 A0A317SIP5 A0A3M2SZ46 A0A4P6XKP0 A0A5N6KSA4 A0A8H6KZQ5 A0A8T9CBW5 A0A9P4IQI3 A0AAD6ISS3 A0AAJ0FJB1 A7TJQ6 B6JZE2 M7P3B3",
    "STE3": "A0A1B7TIY0 A0A1Y2AIP4 A0A2Z6QB68 A0A370TV23 A0A397G363 A0A397TCW3 A0A397TT18 A0A3N2PWD1 A0A3N4LVA6 A0A433PVH9 A0A5C3EYI1 A0A8X7SYL2 A0A9N8YY98 A0A9N9HFT6 A0A9N9NP57 A0A9Q3CV82 A0A9W8Z3T2 A3GFV5 B6K0F1 E4V2W5 G7DZ90 I4YHN7 M5FSK7 S8AI50 S8EMI6",
    "TM116": "A0A151NSC3 A0A1S3HQT7 A0A1X7UMS0 A0A267EA86 A0A267GZC9 A0A2G8L9K8 A0A3M6UV08 A0A3Q1JAQ7 A0A401S6Z2 A0A7M5XJ91 A0A7M7PD29 A0A8B7ZXD8 A0A8J1T644 A0A8S3ZSJ3 A0A9J7LZN7 A0AA35S163 A0AAD8BDH4 A0AAE0SA04 A0AAN0J570 A7RSY4",
    "Vomeronasal1": "A0A3B5R2S4 A0A3P9AWM1 A0A401NRS8 A0A401P5Q4 A0A401SF46 A0A401SF91 A0A452DUS9 A0A4W3GVA8 A0A4W3JR94 A0A4Z2FZ25 A0A7N4PWE8 A0A8B7B4Q3 A0A8C1P5N2 A0A8C4SUP9 A0A8T0AL83 A0A8T2K0N5 A0AA35PA47 A0AA88NKN9 G1U922 H3B062",
}

ACC2FAM = {a: f for f, s in FAMILY_ACCESSIONS.items() for a in s.split()}


# ---------------------------------------------------------------- newick parser
def parse_newick(text):
    """Return (tips, internals) where internals is a list of (label, frozenset_of_tips)."""
    s = "".join(text.split())          # drop all whitespace/newlines
    i, stack, cur, internals, tips = 0, [], [], [], set()
    n = len(s)
    while i < n:
        ch = s[i]
        if ch == '(':
            stack.append(cur); cur = []; i += 1
        elif ch == ',':
            i += 1
        elif ch == ')':
            i += 1
            j = i
            while j < n and s[j] not in ',();:':
                j += 1
            label = s[i:j]
            i = j
            if i < n and s[i] == ':':
                j = i + 1
                while j < n and s[j] not in ',();':
                    j += 1
                i = j
            clade = frozenset().union(*cur) if cur else frozenset()
            internals.append((label, clade))
            if not stack:
                break
            parent = stack.pop(); parent.append(clade); cur = parent
        elif ch == ';':
            break
        else:
            j = i
            while j < n and s[j] not in ',():;':
                j += 1
            name = s[i:j]
            i = j
            if i < n and s[i] == ':':
                j = i + 1
                while j < n and s[j] not in ',();':
                    j += 1
                i = j
            if name:
                cur.append(frozenset([name])); tips.add(name)
    return tips, internals


def accession(tip):
    """tr|ACC|... or tr_ACC_... -> ACC; LIT_GP157|LIT0024|... -> LIT0024"""
    t = tip.replace('|', '_')
    parts = t.split('_')
    if parts[0] == 'LIT':
        return next((p for p in parts if p[3:].isdigit() and p.startswith('LIT')), None)
    return parts[1] if len(parts) > 1 else None


def scf_of(label):
    """Branch label is either '41.6' or 'TBE/sCF' e.g. '56/42.6'. Take the sCF."""
    if not label:
        return None
    try:
        return float(label.split('/')[-1])
    except ValueError:
        return None


# ---------------------------------------------------------------- per-tree summary
def summarise(path):
    tips, internals = parse_newick(open(path, encoding='utf-8').read())
    tip_fam = {}
    unknown = set()
    for t in tips:
        a = accession(t)
        f = ACC2FAM.get(a)
        if f is None:
            unknown.add(a)
        tip_fam[t] = f
    if unknown:
        raise SystemExit(f"ERROR {path}: {len(unknown)} tips not in the family map, "
                         f"e.g. {sorted(x for x in unknown if x)[:5]}")

    ftot = {}
    for f in tip_fam.values():
        ftot[f] = ftot.get(f, 0) + 1          # family sizes AS PRESENT IN THIS TREE
    N = len(tips)

    strata = {'all': [], 'fam2': [], 'fam3': [], 'fam4': []}
    for label, clade in internals:
        v = scf_of(label)
        if v is None:
            continue
        small = clade if len(clade) * 2 <= N else (set(tips) - clade)
        if len(small) < 2 or len(small) > N - 2:
            continue
        counts = {}
        for t in small:
            f = tip_fam[t]
            counts[f] = counts.get(f, 0) + 1
        k = sum(1 for f, c in counts.items() if c > MAJORITY * ftot[f])
        strata['all'].append(v)
        if k >= 2: strata['fam2'].append(v)
        if k >= 3: strata['fam3'].append(v)
        if k >= 4: strata['fam4'].append(v)
    return strata


def tree_index(path):
    base = os.path.basename(path)
    digits = ''.join(c if c.isdigit() else ' ' for c in base).split()
    return int(digits[0]) if digits else -1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dirs', nargs='+', help='directories containing *.cf.tree')
    ap.add_argument('-o', '--out', default='scf_per_tree_summary.tsv')
    args = ap.parse_args()

    rows = []
    for d in args.dirs:
        files = sorted(glob.glob(os.path.join(d, '*.cf.tree')), key=tree_index)
        if not files:
            print(f"  !! no *.cf.tree in {d}", file=sys.stderr); continue
        ds = os.path.basename(os.path.normpath(d))
        print(f"{ds}: {len(files)} trees")
        for p in files:
            s = summarise(p)
            for name, vals in s.items():
                if vals:
                    rows.append((ds, tree_index(p), name, len(vals),
                                 st.mean(vals), st.median(vals)))

    with open(args.out, 'w', encoding='utf-8', newline='') as fh:
        fh.write("dataset\ttree\tstratum\tn_branches\tmean_sCF\tmedian_sCF\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{r[4]:.3f}\t{r[5]:.3f}\n")
    print(f"\nwrote {args.out}  ({len(rows)} rows)")

    print("\n=== per-dataset summary across trees ===")
    print(f"{'dataset':<46} {'stratum':>7} {'trees':>6} {'branches/tree':>14} "
          f"{'mean of means':>14} {'SD':>7} {'range':>15}")
    for ds in dict.fromkeys(r[0] for r in rows):
        for stratum in ('all', 'fam2', 'fam3', 'fam4'):
            v = [r[4] for r in rows if r[0] == ds and r[2] == stratum]
            nb = [r[3] for r in rows if r[0] == ds and r[2] == stratum]
            if not v: continue
            sd = st.stdev(v) if len(v) > 1 else 0.0
            print(f"{ds[:46]:<46} {stratum:>7} {len(v):>6} {st.mean(nb):>14.1f} "
                  f"{st.mean(v):>14.2f} {sd:>7.2f} {min(v):>7.2f}-{max(v):<7.2f}")


if __name__ == '__main__':
    main()
