#!/usr/bin/env python3
"""Compare saved opt-in phase timers; nested phases are not additive."""

import argparse
import csv
from pathlib import Path


def read(path):
    result={}
    for line in path.read_text().splitlines():
        fields=line.split(",")
        if fields[0]=="XMVB_PHASE" and fields[2]!="calls":
            result.setdefault(fields[1],{}).update(
                calls=int(fields[2]),wall_s=float(fields[3]))
        elif fields[0]=="XMVB_CPU" and fields[2]!="cpu_seconds":
            result.setdefault(fields[1],{})["cpu_s"]=float(fields[2])
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline",type=Path)
    parser.add_argument("candidate",type=Path)
    parser.add_argument("output",type=Path)
    args=parser.parse_args()
    old,new=read(args.baseline),read(args.candidate)
    rows=[]
    for phase in dict.fromkeys((*old,*new)):
        row=dict(phase=phase)
        for label,source in (("baseline",old),("candidate",new)):
            values=source.get(phase,{})
            for key in ("calls","wall_s","cpu_s"):
                row[label+"_"+key]=values.get(key,"")
        row["speedup"]=(row["baseline_wall_s"]/row["candidate_wall_s"]
                        if row["baseline_wall_s"]!="" and row["candidate_wall_s"] else "")
        rows.append(row)
    with args.output.open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__=="__main__":
    main()
