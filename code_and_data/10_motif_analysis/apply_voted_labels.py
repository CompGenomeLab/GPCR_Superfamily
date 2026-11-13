"""Replace the single-reference GPCRdb labels of logo mapping JSONs with the majority label of the
family's human receptors (label_by_human_receptors.py) and its agreement.

Every position of a family gets `gpcrdb` = the majority label at its superfamily column,
`gpcrdbAgreement` = the majority's share and `gpcrdbVoters` = the number of receptors voting;
a position without voters has no label. Kingdom rows take their family's labels. The family
reference's residue (refResidue, refAA) is left as it is.

python apply_voted_labels.py <logo mapping dir> <labels tsv>
"""
import csv
import json
import sys
from pathlib import Path


def main(logo_dir, labels_tsv):
    votes = {}
    for r in csv.DictReader(open(labels_tsv), delimiter="\t"):
        votes[(r["family"], int(r["sup_column"]) - 1)] = (r["label"], float(r["agreement"]), int(r["n_voters"]))
    for path in sorted(Path(logo_dir).glob("*.json")):
        if path.name == "index.json":
            continue
        data = json.loads(path.read_text())
        family = data["familyKey"]
        if family.split("_")[-1] in ("Animals", "Fungi", "Plants"):
            family = family.split("_")[0]
        labelled = 0
        for c, p in enumerate(data["positions"]):
            for k in ("gpcrdb", "gpcrdbAgreement", "gpcrdbVoters"):
                p.pop(k, None)
            v = votes.get((family, c))
            if v:
                p["gpcrdb"], p["gpcrdbAgreement"], p["gpcrdbVoters"] = v
                labelled += 1
        data["labels"] = "majority of human receptors with GPCRdb numbering"
        path.write_text(json.dumps(data))
        print(f"{path.name}: {labelled} positions labelled")


if __name__ == "__main__":
    main(*sys.argv[1:3])
