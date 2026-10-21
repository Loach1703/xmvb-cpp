"""Parse current outputs and saved controls; never launch an executable."""
import csv
import hashlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "tests/mcscf"))
from c12_spaces import summarize


def read_csv(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def write_csv(name, rows):
    with (HERE / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def phases(path):
    with path.open() as stream:
        return {row[1]: (int(row[2]), float(row[3]))
                for row in csv.reader(stream)
                if len(row) == 4 and row[0] == "XMVB_PHASE" and row[2].isdigit()}


current = {r["space"]: r for r in summarize(HERE / "native")
           if r["space"] in ("cov", "ion02", "full")}
prior = {r["space"]: r for r in read_csv(
    HERE.parent / "c12-amplitude-20261017/comparison.csv")}
best_full = next(r for r in read_csv(
    HERE.parent / "c12-spaces-20261017/summary.csv") if r["space"] == "full")
comparison = []
for label in ("cov", "ion02", "full"):
    new, old = current[label], prior[label]
    assert new["status"] == "converged" and new["failed_hvp"] == 0
    run = HERE / "native" / label
    digest = hashlib.sha256((run / "input.xmi").read_bytes()).hexdigest()
    assert digest == old["input_sha256"]
    best = best_full if label == "full" else {
        "application_s": old["baseline_application_s"],
        "process_s": old["baseline_process_s"], "rss_GiB": old["baseline_rss_GiB"],
        "final_Eh": old["baseline_Eh"], "steps": old["baseline_steps"],
        "hvp": old["baseline_hvp"], "node": old["baseline_node"]}
    row = dict(space=label, structures=new["structures"],
               baseline_job="1265093_2" if label == "full" else old["baseline_job"],
               current_job=f"1265585_{dict(cov=6, ion02=1, full=2)[label]}",
               baseline_node=best["node"], current_node=new["node"],
               baseline_steps=int(best["steps"]), current_steps=new["steps"],
               baseline_hvp=int(best["hvp"]), current_hvp=new["hvp"],
               baseline_Eh=float(best["final_Eh"]), current_Eh=float(new["final_Eh"]),
               delta_Eh=float(new["final_Eh"])-float(best["final_Eh"]),
               baseline_application_s=float(best["application_s"]), current_application_s=new["application_s"],
               baseline_process_s=float(best["process_s"]), current_process_s=new["process_s"],
               baseline_rss_GiB=float(best["rss_GiB"]), current_rss_GiB=new["rss_GiB"],
               time_saved_percent=100*(1-new["application_s"]/float(best["application_s"])),
               parent_application_s=float(old["current_application_s"]),
               parent_time_saved_percent=100*(1-new["application_s"]/float(old["current_application_s"])),
               input_sha256=digest,
               affinity=(run / "affinity.txt").read_text().strip().split()[-1])
    comparison.append(row)
write_csv("native.csv", list(current.values()))
write_csv("comparison.csv", comparison)
before = phases(HERE.parent / "c12-amplitude-20261017/profile/run.err")
after = phases(HERE / "profile/run.err")
write_csv("phases.csv", [dict(phase=key, baseline_calls=before[key][0],
    current_calls=after[key][0], baseline_s=before[key][1], current_s=after[key][1])
    for key in ("hvp_apply", "ci_sigma", "ci_transform_left", "ci_transform_right",
                "ci_response", "ciah_macro", "ciah_keyframe", "ciah_trial",
                "ci_product_projection", "ci_product_generate")])
for row in comparison:
    print(f'{row["space"]}: {row["current_steps"]} steps / {row["current_hvp"]} HVP, '
          f'{row["current_Eh"]:.12f} Eh, {row["baseline_application_s"]:.6f} -> '
          f'{row["current_application_s"]:.6f} s, '
          f'{row["baseline_rss_GiB"]:.6f} -> {row["current_rss_GiB"]:.6f} GiB, '
          f'saved {row["time_saved_percent"]:.2f}%')
