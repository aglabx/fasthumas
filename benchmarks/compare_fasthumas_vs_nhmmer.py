#!/usr/bin/env python3
"""
Comprehensive concordance analysis between FastHumAS and legacy HumAS-HMMER (nhmmer)
across finished T2T-CHM13 chromosomes.

Evaluates:
- Discovery union (all detected alpha-satellite loci)
- Mutual overlap (reciprocal overlap >= 50% and >= 80%)
- Boundary concordances (exact 0 bp, <= 3 bp, <= 10 bp)
- Label / subfamily concordance on overlapping records
- Strand concordance
- Sensitivity, precision, and discordant (method-unique) loci characteristics
"""

import sys
import os
import glob
from collections import defaultdict

def load_bed(path):
    records = []
    if not os.path.exists(path):
        return records
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 6:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                name = parts[3]
                score = float(parts[4])
                strand = parts[5]
                records.append({
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "length": end - start,
                    "name": name,
                    "score": score,
                    "strand": strand,
                })
    return records

def analyze_chromosome(chrom, fasthumas_recs, legacy_recs):
    f_sorted = sorted(fasthumas_recs, key=lambda x: (x["start"], x["end"]))
    l_sorted = sorted(legacy_recs, key=lambda x: (x["start"], x["end"]))

    # Find all candidate pairs with true reciprocal overlap >= 50%
    # (ov / len(l) >= 0.5 AND ov / len(f) >= 0.5)
    candidate_pairs = []
    f_idx = 0
    for l_i, l_rec in enumerate(l_sorted):
        l_s, l_e = l_rec["start"], l_rec["end"]
        l_len = l_rec["length"]
        if l_len <= 0:
            continue

        while f_idx < len(f_sorted) and f_sorted[f_idx]["end"] <= l_s:
            f_idx += 1

        curr = f_idx
        while curr < len(f_sorted) and f_sorted[curr]["start"] < l_e:
            f_rec = f_sorted[curr]
            f_len = f_rec["length"]
            if f_len > 0:
                ov_s = max(l_s, f_rec["start"])
                ov_e = min(l_e, f_rec["end"])
                ov = max(0, ov_e - ov_s)
                if (ov / l_len >= 0.5) and (ov / f_len >= 0.5):
                    candidate_pairs.append((ov, l_i, curr))
            curr += 1

    # 1:1 greedy matching by maximum overlap length to prevent many-to-one inflation
    candidate_pairs.sort(key=lambda x: x[0], reverse=True)
    matched_l = set()
    matched_f = set()
    paired = []

    for ov, l_i, f_i in candidate_pairs:
        if l_i not in matched_l and f_i not in matched_f:
            matched_l.add(l_i)
            matched_f.add(f_i)
            paired.append((l_sorted[l_i], f_sorted[f_i], ov))

    exact_boundary = 0
    within_3bp_boundary = 0
    within_10bp_boundary = 0
    label_exact_match = 0
    strand_match = 0
    ov_80_count = 0

    for l_rec, f_rec, ov in paired:
        if (ov / l_rec["length"] >= 0.8) and (ov / f_rec["length"] >= 0.8):
            ov_80_count += 1
        d_start = abs(l_rec["start"] - f_rec["start"])
        d_end = abs(l_rec["end"] - f_rec["end"])

        if d_start == 0 and d_end == 0:
            exact_boundary += 1
        if d_start <= 3 and d_end <= 3:
            within_3bp_boundary += 1
        if d_start <= 10 and d_end <= 10:
            within_10bp_boundary += 1

        if l_rec["name"] == f_rec["name"]:
            label_exact_match += 1
        if l_rec["strand"] == f_rec["strand"]:
            strand_match += 1

    n_f = len(f_sorted)
    n_l = len(l_sorted)
    n_matched = len(paired)

    return {
        "chrom": chrom,
        "n_fasthumas": n_f,
        "n_legacy": n_l,
        "n_matched": n_matched,
        "n_ov80": ov_80_count,
        "exact_boundary": exact_boundary,
        "within_3bp_boundary": within_3bp_boundary,
        "within_10bp_boundary": within_10bp_boundary,
        "label_exact_match": label_exact_match,
        "strand_match": strand_match,
        "f_unique": n_f - n_matched,
        "l_unique": n_l - n_matched,
    }

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <fasthumas_bed> <legacy_dir>")
        sys.exit(1)

    fasthumas_path = sys.argv[1]
    legacy_dir = sys.argv[2]

    print("Loading FastHumAS BED...")
    all_f_recs = load_bed(fasthumas_path)
    f_by_chrom = defaultdict(list)
    for r in all_f_recs:
        f_by_chrom[r["chrom"]].append(r)

    legacy_files = glob.glob(os.path.join(legacy_dir, "AS-HOR+SF-vs-*.bed"))
    if not legacy_files:
        print(f"No legacy BED files found in {legacy_dir}!")
        sys.exit(1)

    print(f"Found {len(legacy_files)} finished legacy chromosome files in {legacy_dir}")
    print("=" * 80)
    print(f"{'Chrom':<14} | {'FastHumAS':<10} | {'Legacy':<8} | {'Overlap':<8} | {'Sens %':<7} | {'Prec %':<7} | {'<=3bp Bdry':<11} | {'Label Agr':<10} | {'Strand Agr':<10}")
    print("-" * 80)

    total_f = 0
    total_l = 0
    total_ov = 0
    total_exact = 0
    total_3bp = 0
    total_10bp = 0
    total_label = 0
    total_strand = 0
    total_f_uniq = 0
    total_l_uniq = 0

    for leg_file in sorted(legacy_files):
        # filename format: AS-HOR+SF-vs-NC_060946.1.bed
        base = os.path.basename(leg_file)
        chrom = base.replace("AS-HOR+SF-vs-", "").replace(".bed", "")
        leg_recs = load_bed(leg_file)
        fasthumas_recs = f_by_chrom[chrom]

        stats = analyze_chromosome(chrom, fasthumas_recs, leg_recs)

        sens = (stats["n_matched"] / stats["n_legacy"] * 100.0) if stats["n_legacy"] > 0 else 0.0
        prec = (stats["n_matched"] / stats["n_fasthumas"] * 100.0) if stats["n_fasthumas"] > 0 else 0.0
        bdry3_pct = (stats["within_3bp_boundary"] / stats["n_matched"] * 100.0) if stats["n_matched"] > 0 else 0.0
        label_pct = (stats["label_exact_match"] / stats["n_matched"] * 100.0) if stats["n_matched"] > 0 else 0.0
        strand_pct = (stats["strand_match"] / stats["n_matched"] * 100.0) if stats["n_matched"] > 0 else 0.0

        print(f"{chrom:<14} | {stats['n_fasthumas']:<10} | {stats['n_legacy']:<8} | {stats['n_matched']:<8} | {sens:6.2f}% | {prec:6.2f}% | {bdry3_pct:6.2f}% ({stats['within_3bp_boundary']}) | {label_pct:6.2f}% | {strand_pct:6.2f}%")

        total_f += stats["n_fasthumas"]
        total_l += stats["n_legacy"]
        total_ov += stats["n_matched"]
        total_exact += stats["exact_boundary"]
        total_3bp += stats["within_3bp_boundary"]
        total_10bp += stats["within_10bp_boundary"]
        total_label += stats["label_exact_match"]
        total_strand += stats["strand_match"]
        total_f_uniq += stats["f_unique"]
        total_l_uniq += stats["l_unique"]

    print("=" * 80)
    total_sens = (total_ov / total_l * 100.0) if total_l > 0 else 0.0
    total_prec = (total_ov / total_f * 100.0) if total_f > 0 else 0.0
    total_bdry3 = (total_3bp / total_ov * 100.0) if total_ov > 0 else 0.0
    total_label_pct = (total_label / total_ov * 100.0) if total_ov > 0 else 0.0
    total_strand_pct = (total_strand / total_ov * 100.0) if total_ov > 0 else 0.0

    print(f"{'TOTAL / MEAN':<14} | {total_f:<10} | {total_l:<8} | {total_ov:<8} | {total_sens:6.2f}% | {total_prec:6.2f}% | {total_bdry3:6.2f}% ({total_3bp}) | {total_label_pct:6.2f}% | {total_strand_pct:6.2f}%")
    print(f"\nExact Boundary (=0 bp): {total_exact} ({total_exact/total_ov*100:.2f}%)")
    print(f"Boundary within <=10 bp: {total_10bp} ({total_10bp/total_ov*100:.2f}%)")
    print(f"FastHumAS unique (not in legacy): {total_f_uniq} ({total_f_uniq/total_f*100:.2f}%)")
    print(f"Legacy unique (not in FastHumAS): {total_l_uniq} ({total_l_uniq/total_l*100:.2f}%)")

if __name__ == "__main__":
    main()
