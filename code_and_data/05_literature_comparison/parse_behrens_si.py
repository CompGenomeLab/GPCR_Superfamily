"""Behrens et al. 2023 SI Table S1 (pages 5-7): cartilaginous fish T2R1, reference T2R and outgroup
ORA/V1R sequences, read in pypdf layout mode by the column position of the sequence field."""
import logging, re, sys
import pypdf

logging.disable(logging.WARNING)
PDF = sys.argv[1]
OUT = sys.argv[2]
SECTIONS = {4: "cartilaginous fish T2R1", 5: "reference T2R", 6: "outgroup ORA/V1R"}

r = pypdf.PdfReader(PDF)
records = []
for page, section in SECTIONS.items():
    lines = r.pages[page].extract_text(extraction_mode="layout").splitlines()
    head = next(i for i, l in enumerate(lines) if "Amino acid sequence" in l)
    s0 = lines[head].index("Amino acid sequence")
    nxt = [m.start() for m in re.finditer(r"\S+(?: \S+)*", lines[head]) if m.start() > s0 + 10]
    s1 = nxt[0] if nxt else 10**6
    block = []
    for l in lines[head + 1:] + [""]:
        if l.strip():
            block.append(l)
            continue
        if not block:
            continue
        mi = [i for i, b in enumerate(block) if b[:1].strip()]
        # rows without a blank line between them: split at the midpoint between species lines
        # labels are vertically centred on their own rows: row k starts where row k-1 ends and spans
        # twice the distance from its start to its label, plus the label line
        cuts, start = [0], 0
        for m in mi[:-1]:
            start = start + 2 * (m - start) + 1
            cuts.append(start)
        cuts.append(len(block))
        for k, m in enumerate(mi):
            part = block[cuts[k]:cuts[k + 1]]
            seq = "".join(re.sub(r"[^A-Z*!]", "", b[s0 - 2:s1 - 1]) for b in part)
            left = re.split(r"\s{2,}", block[m][:s0 - 2].strip())
            right = " ".join(re.sub(r"\s+", " ", b[s1 - 1:]).strip() for b in part).strip()
            if seq:
                records.append({"section": section, "species": left[0], "fields": left[1:], "right": right, "seq": seq})
        block = []

with open(OUT, "w", encoding="utf-8") as fh:
    for k, x in enumerate(records, 1):
        genes = [f for f in x["fields"] if re.search(r"T2R|ORA|V1R|v1r|^[A-Z][a-z]+\d", f)]
        name = (genes[0] if genes else (x["fields"][-1] if x["fields"] else "rec%d" % k)).replace(" ", "")
        fh.write(">%s|%s|%s|%s\n%s\n" % (name, x["species"].replace(" ", "_"), x["section"].replace(" ", "_"),
                                        " ".join(x["fields"][:-1]).replace(" ", "_"), x["seq"]))
for x in records:
    print(x["section"][:12], "|", x["species"], "|", x["fields"], "|", len(x["seq"]), x["seq"][:25], "|", x["right"][:60])
print(len(records), "records")
