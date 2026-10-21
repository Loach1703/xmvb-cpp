#!/usr/bin/env python3
"""Compare new native runs with saved results; never execute old binaries."""

import argparse
import csv
from pathlib import Path


def read(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("candidate",type=Path,nargs="+")
    parser.add_argument("--baseline",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    base=args.baseline
    refs={}
    for filename in ("direct-diagonal-native.csv","matched-exact-c12.csv",
                     "matched-exact-quick.csv"):
        for row in read(base/filename):
            refs[row["case"],row["mode"]]=row
    # No saved matched-node exact C14 run exists. Keep that limitation explicit.
    for row in read(base/"exact-references.csv"):
        refs.setdefault((row["case"],"exact"),row)
    rows=[]
    for path in args.candidate:
        for new in read(path):
            old=refs[new["case"],new["mode"]]
            if old["input_sha256"]!=new["input_sha256"]:
                raise RuntimeError(f"input mismatch: {new['case']}")
            same_node_mask=(old["node"]==new["node"] and
                            old.get("affinity")==new["affinity"])
            row=dict(case=new["case"],mode=new["mode"],
                     baseline_job=old["job"],candidate_job=new["job"],
                     baseline_node=old["node"],candidate_node=new["node"],
                     same_node_mask=same_node_mask)
            for key in ("steps","hvp","process_s","rss_GiB","final_Eh",
                        "ci_actions","ci_s","hvp_s"):
                row["baseline_"+key]=old.get(key,"")
                row["candidate_"+key]=new[key]
            row["speedup"]=float(old["process_s"])/float(new["process_s"])
            row["delta_Eh"]=float(new["final_Eh"])-float(old["final_Eh"])
            row["converged"]=new["converged"]
            rows.append(row)
    with args.output.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(f"{row['case']} {row['mode']}: "
              f"{row['baseline_process_s']} -> {row['candidate_process_s']} s, "
              f"{row['speedup']:.3f}x; "
              f"steps/HVP={row['candidate_steps']}/{row['candidate_hvp']}; "
              f"dE={row['delta_Eh']:.3g} Eh; matched={row['same_node_mask']}")


if __name__=="__main__":
    main()
