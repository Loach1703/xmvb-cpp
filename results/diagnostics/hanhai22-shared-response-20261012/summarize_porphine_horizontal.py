"""Compare completed original-guess porphine jobs; parse artifacts only."""
import csv
import json
from pathlib import Path
import re

from summarize_horizontal import measure


ROOT = Path(__file__).resolve().parent


def main():
    before = measure(ROOT / "porphine-transport-1254189/singlet")
    after = measure(ROOT / "porphine-horizontal-1254797/singlet")
    reference = ROOT / "porphine-rotation-1248772/singlet/pyscf-ah"
    pyscf = json.loads((reference / "result.json").read_text())
    resource = (reference / "resource.time").read_text()
    output = (reference / "output.out").read_text()
    counts = re.search(r"CASSCF converged in (\d+) macro \((\d+) JK (\d+) micro\)", output)
    wall = re.search(r"Elapsed \(wall clock\) time.*?:\s*([\d:.]+)", resource)[1]
    wall_seconds = sum(float(part)*60**i for i, part in enumerate(wall.split(":")[::-1]))
    assert before["converged"] and after["converged"] and pyscf["converged"]
    assert before["exit_code"] == after["exit_code"] == 0
    rows = []
    for method, job, value in (("NEO previous", "1254189", before),
                               ("NEO S-horizontal", "1254797", after)):
        rows.append(dict(method=method, job=job, energy=value["energy"],
                         macros=value["iterations"], hvp=value["hvp"],
                         wall_seconds=value["wall_seconds"],
                         rss_gib=value["rss_mib"]/1024,
                         delta_pyscf_microhartree=(value["energy"]-pyscf["energy"])*1e6))
    rows.append(dict(method="PySCF default AH", job="1248772",
                     energy=pyscf["energy"], macros=int(counts[1]), hvp=None,
                     wall_seconds=wall_seconds,
                     rss_gib=int(re.search(r"Maximum resident set size \(kbytes\):\s*(\d+)", resource)[1])/1024**2,
                     delta_pyscf_microhartree=0.0))
    geometry = json.loads((ROOT / "representative-analysis-1254863/representative-analysis.json").read_text())
    summary = dict(comparison=rows,
                   pyscf_actions=dict(jk=int(counts[2]), micro=int(counts[3])),
                   terminal_geometry=geometry["between"]["neo_current_terminal_pyscf_low"],
                   qualification="Energy agreement does not establish identical active subspaces or wavefunctions.")
    target = ROOT / "porphine-horizontal-1254797"
    (target / "comparison.json").write_text(json.dumps(summary, indent=2)+"\n")
    with (target / "comparison.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
