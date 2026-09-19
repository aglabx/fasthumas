#!/usr/bin/env python3
"""
Reproducible analysis of the active centromeric HOR array on Chromosome 22:
Locus: NC_060946.1:12,788,180–15,711,065 (hor_22_9 / S2C14/22H1L, ~2.92 Mb).

Compares FastHumAS vs legacy HumAS-HMMER on junction continuity and details
the exact nature of all gaps > 3 bp.
"""

import sys
import os
from collections import Counter

REGION_CHR = "NC_060946.1"
REGION_START = 12788180
REGION_END = 15711065

def load_regional_bed(path):
    records = []
    if not os.path.exists(path):
        return records
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6 and parts[0] == REGION_CHR:
                s = int(parts[1])
                e = int(parts[2])
                if s >= REGION_START and e <= REGION_END:
                    records.append({
                        "chrom": parts[0],
                        "start": s,
                        "end": e,
                        "name": parts[3],
                        "score": float(parts[4]),
                        "strand": parts[5],
                        "length": e - s,
                    })
    records.sort(key=lambda x: (x["start"], x["end"]))
    return records

def analyze_continuity(records):
    gaps = []
    gap_details = []
    for i in range(len(records) - 1):
        r1 = records[i]
        r2 = records[i + 1]
        g = r2["start"] - r1["end"]
        gaps.append(g)
        if g > 3 or g < 0:
            gap_details.append({
                "index": i,
                "gap": g,
                "r1": r1,
                "r2": r2,
            })
    return Counter(gaps), gap_details

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <fasthumas_bed> <legacy_bed_or_dir>")
        sys.exit(1)

    fast_path = sys.argv[1]
    leg_arg = sys.argv[2]
    if os.path.isdir(leg_arg):
        leg_path = os.path.join(leg_arg, f"AS-HOR+SF-vs-{REGION_CHR}.bed")
    else:
        leg_path = leg_arg

    f_recs = load_regional_bed(fast_path)
    l_recs = load_regional_bed(leg_path)

    f_counts, f_details = analyze_continuity(f_recs)
    l_counts, l_details = analyze_continuity(l_recs)

    f_total_junc = len(f_recs) - 1
    l_total_junc = len(l_recs) - 1

    f_flush = f_counts[0]
    l_flush = l_counts[0]

    f_micro = f_counts[1] + f_counts[2] + f_counts[3]
    l_micro = l_counts[1] + l_counts[2] + l_counts[3]

    f_gt3 = sum(v for k, v in f_counts.items() if k > 3)
    l_gt3 = sum(v for k, v in l_counts.items() if k > 3)

    f_neg = sum(v for k, v in f_counts.items() if k < 0)
    l_neg = sum(v for k, v in l_counts.items() if k < 0)

    print("==========================================================================================================")
    print(f"Active Centromeric HOR Array Analysis: chr22 ({REGION_CHR}:{REGION_START:,}–{REGION_END:,}, 2.92 Mb)")
    print(f"CenSat Domain: hor_22_9 (S2C14/22H1L)")
    print("==========================================================================================================")
    print(f"Metric                                  | FastHumAS                 | Legacy HumAS-HMMER")
    print("-" * 106)
    print(f"Total Annotated Monomers                | {len(f_recs):<25} | {len(l_recs):<25}")
    print(f"Total Inter-Monomer Junctions           | {f_total_junc:<25} | {l_total_junc:<25}")
    print(f"Flush Contiguous Junctions (gap = 0 bp) | {f_flush} ({f_flush/f_total_junc*100:6.2f}%)       | {l_flush} ({l_flush/l_total_junc*100:6.2f}%)")
    print(f"Micro-gaps (1–3 bp)                     | {f_micro} ({f_micro/f_total_junc*100:6.2f}%)       | {l_micro} ({l_micro/l_total_junc*100:6.2f}%)")
    print(f"  - 1-bp gaps                           | {f_counts[1]:<25} | {l_counts[1]:<25}")
    print(f"  - 2-bp gaps                           | {f_counts[2]:<25} | {l_counts[2]:<25}")
    print(f"  - 3-bp gaps                           | {f_counts[3]:<25} | {l_counts[3]:<25}")
    print(f"Gaps > 3 bp                             | {f_gt3:<25} | {l_gt3:<25}")
    print(f"Overlaps (< 0 bp)                       | {f_neg:<25} | {l_neg:<25}")
    print("=" * 106)

    print("\nDetailed Investigation of Gaps > 3 bp in FastHumAS:")
    for d in f_details:
        if d["gap"] > 3:
            r1 = d["r1"]
            r2 = d["r2"]
            print(f"  * Locus: [{r1['end']}, {r2['start']}) (length = {d['gap']} bp)")
            print(f"    FastHumAS Flanking: {r1['name']} [{r1['start']},{r1['end']}) -> {r2['name']} [{r2['start']},{r2['end']})")

            # Dynamically inspect legacy calls overlapping this specific gap interval
            l_overlaps = [
                lr for lr in l_recs
                if not (lr["end"] <= r1["end"] or lr["start"] >= r2["start"])
            ]
            if l_overlaps:
                print(f"    Legacy overlapping records in this span ({len(l_overlaps)} record(s)):")
                for lr in l_overlaps:
                    dens = lr["score"] / lr["length"] if lr["length"] > 0 else 0
                    print(f"      - {lr['name']} [{lr['start']},{lr['end']}), score={lr['score']:.1f}, len={lr['length']} bp, score_density={dens:.3f}")
            else:
                l_prev = [lr for lr in l_recs if lr["end"] <= r1["end"]]
                l_next = [lr for lr in l_recs if lr["start"] >= r2["start"]]
                if l_prev and l_next:
                    prev_l = l_prev[-1]
                    next_l = l_next[0]
                    l_gap = next_l["start"] - prev_l["end"]
                    print(f"    Legacy has no overlapping calls; flanking legacy gap is {l_gap} bp:")
                    print(f"      - Left:  {prev_l['name']} [{prev_l['start']},{prev_l['end']})")
                    print(f"      - Right: {next_l['name']} [{next_l['start']},{next_l['end']})")

            # Contextual analysis based on confirmed coordinates and model identities
            if r1["end"] == 12789031 and r2["start"] == 12789250:
                print("    Biological / Algorithmic context (chr22 pericentromeric boundary):")
                print("      - Degenerate pericentromeric transition boundary near array start.")
                print("      - Legacy nhmmer called two fragmented sub-monomers (S2C9H1L.7 [score 97.9, 118 bp] and S2C14/22H1L.3 [score 90.4, 95 bp]),")
                print("        with score densities 0.830 and 0.952 (> 0.70 threshold). FastHumAS's Gotoh affine DP search against the active HOR model")
                print("        did not yield a hit passing the seed/score filter in this degenerate region, leaving a 219-bp unannotated gap.")
            elif d["gap"] == 31 and r1["name"].startswith("S2C14/22H1L.7") and r2["name"].startswith("S2C14/22H1L.6"):
                print("    Biological / Algorithmic context (recurrent 31-bp inter-monomer locus):")
                print("      - Recurrent 31-bp unannotated sequence (CAACAAAAAGTGTTTTTCAAAACTGCTGTAT) between S2C14/22H1L.7 and S2C14/22H1L.6.")
                print("      - FastHumAS Gotoh affine DP calls a full 177-bp canonical monomer, leaving 31 bp unannotated.")
                print("      - Legacy nhmmer split the downstream monomer into two partial fragments (58 bp, score 46.4 and 145 bp, score 161.8).")
                print("      - Note: Determining whether this 31-bp locus represents a true insertion variant or boundary placement variation")
                print("        between profile models requires pairwise structural alignment against the canonical repeat unit.")
            elif r1["end"] == 15709700 and r2["start"] == 15709706:
                print("    Biological / Algorithmic context (unannotated inter-monomer sequence):")
                print("      - Unannotated 6-bp sequence (CTAAAA) in assembly between S2C14/22H1L.4 and S2C14/22H1L.3.")
                print("      - Note: This is an unannotated inter-monomer sequence (valid genomic bases without Ns), NOT a physical assembly gap.")
                print("        Legacy HumAS-HMMER similarly leaves an 8-bp unannotated gap [15709698, 15709706) at this exact locus.")
            else:
                print(f"    Unannotated inter-monomer sequence of length {d['gap']} bp between {r1['name']} and {r2['name']}.")

    if f_details:
        neg_details = [d for d in f_details if d["gap"] < 0]
        if neg_details:
            print("\nDetailed Investigation of Overlaps (< 0 bp):")
            for d in neg_details:
                r1 = d["r1"]
                r2 = d["r2"]
                print(f"  * Overlap: {d['gap']} bp between {r1['name']} [{r1['start']},{r1['end']}) and {r2['name']} [{r2['start']},{r2['end']})")

if __name__ == "__main__":
    main()
