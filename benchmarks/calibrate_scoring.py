#!/usr/bin/env python3
"""
Scoring model calibration analysis: FastHumAS Gotoh PSSM log-odds vs HMMER bit scores.

Demonstrates that for 171-bp alpha-satellite monomer profiles, the empirical scaling factor
of 0.85 linearly maps pure Gotoh match log-odds sums to HMMER's Null2/transition-penalized
bit scores with high precision:
- Training split (chr1, chr3, chr8; n = 23,527): Mean 0.8495 +/- 0.0074 (SD), Median 0.8504
- Independent test split (chr10, 11, 12, 14, 22, Y; n = 34,656): Mean 0.8483 +/- 0.0085 (SD), Median 0.8500
"""

import sys
import os
import glob
import numpy as np

def load_bed(path):
    records = []
    if not os.path.exists(path):
        return records
    with open(path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 6:
                records.append({
                    "chrom": parts[0],
                    "start": int(parts[1]),
                    "end": int(parts[2]),
                    "name": parts[3],
                    "score": float(parts[4]),
                    "strand": parts[5]
                })
    return records

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <fasthumas_bed> <legacy_dir> [sample_out_tsv]")
        sys.exit(1)

    fasthumas_path = sys.argv[1]
    legacy_dir = sys.argv[2]
    sample_out_path = sys.argv[3] if len(sys.argv) > 3 else None

    f_recs = load_bed(fasthumas_path)
    f_dict = {(r["chrom"], r["start"], r["end"], r["strand"]): r for r in f_recs}

    legacy_files = sorted(glob.glob(os.path.join(legacy_dir, "AS-HOR+SF-vs-*.bed")))
    if not legacy_files:
        print("No legacy files found!")
        sys.exit(1)

    train_chroms = {"NC_060925.1", "NC_060927.1", "NC_060932.1"} # chr1, chr3, chr8
    test_chroms = {"NC_060934.1", "NC_060935.1", "NC_060936.1", "NC_060938.1", "NC_060941.1", "NC_060946.1", "NC_060948.1"}

    train_ratios = []
    test_ratios = []
    all_pairs = []

    for l_f in legacy_files:
        chrom = os.path.basename(l_f).replace("AS-HOR+SF-vs-", "").replace(".bed", "")
        l_recs = load_bed(l_f)
        for l in l_recs:
            k = (chrom, l["start"], l["end"], l["strand"])
            if k in f_dict:
                f = f_dict[k]
                # Require model name agreement for precise profile calibration
                if l["name"] != f["name"]:
                    continue
                raw_gotoh = f["score"] / 0.85
                ratio = l["score"] / raw_gotoh
                all_pairs.append({
                    "chrom": chrom,
                    "start": l["start"],
                    "end": l["end"],
                    "name": l["name"],
                    "hmmer_score": l["score"],
                    "fasthumas_score": f["score"],
                    "raw_gotoh_score": raw_gotoh,
                    "ratio": ratio,
                    "split": "train" if chrom in train_chroms else "test"
                })
                if chrom in train_chroms:
                    train_ratios.append(ratio)
                elif chrom in test_chroms:
                    test_ratios.append(ratio)

    train_r = np.array(train_ratios)
    test_r = np.array(test_ratios)

    print("================================================================================")
    print("FastHumAS Gotoh DP vs HMMER Bit Score Calibration (Name-Matched)")
    print(f"Total evaluated exact-boundary, model-matched monomers: {len(all_pairs)}")
    print("================================================================================")
    print(f"Training set (chr1, chr3, chr8): n = {len(train_r)}")
    print(f"  Mean ratio (HMMER / raw Gotoh): {train_r.mean():.4f} +/- {train_r.std():.4f} (SD)")
    print(f"  Standard Error of Mean (SEM):  {train_r.std() / np.sqrt(len(train_r)):.6f}")
    print(f"  Median ratio:                  {np.median(train_r):.4f}")
    print(f"  Interquartile Range (IQR):     [{np.percentile(train_r, 25):.4f}, {np.percentile(train_r, 75):.4f}]")
    print(f"  95% empirical range:           [{np.percentile(train_r, 2.5):.4f}, {np.percentile(train_r, 97.5):.4f}]")

    print(f"\nIndependent Validation set (chr10, 11, 12, 14, 17, 22, Y): n = {len(test_r)}")
    print(f"  Mean ratio (HMMER / raw Gotoh): {test_r.mean():.4f} +/- {test_r.std():.4f} (SD)")
    print(f"  Standard Error of Mean (SEM):  {test_r.std() / np.sqrt(len(test_r)):.6f}")
    print(f"  Median ratio:                  {np.median(test_r):.4f}")
    print(f"  Interquartile Range (IQR):     [{np.percentile(test_r, 25):.4f}, {np.percentile(test_r, 75):.4f}]")
    print(f"  95% empirical range:           [{np.percentile(test_r, 2.5):.4f}, {np.percentile(test_r, 97.5):.4f}]")

    # Stratified analysis by score range across all pairs
    all_hmmer = np.array([p["hmmer_score"] for p in all_pairs])
    all_r = np.array([p["ratio"] for p in all_pairs])
    high_mask = all_hmmer >= 150.0
    mid_mask = (all_hmmer >= 100.0) & (all_hmmer < 150.0)
    low_mask = all_hmmer < 100.0

    print("\nScore-Stratified Calibration Analysis (All Chromosomes):")
    if np.any(high_mask):
        print(f"  High score (>= 150 bits, n = {np.sum(high_mask)}): Mean ratio = {all_r[high_mask].mean():.4f} +/- {all_r[high_mask].std():.4f} (Median = {np.median(all_r[high_mask]):.4f})")
    if np.any(mid_mask):
        print(f"  Mid score  (100–150 bits, n = {np.sum(mid_mask)}): Mean ratio = {all_r[mid_mask].mean():.4f} +/- {all_r[mid_mask].std():.4f} (Median = {np.median(all_r[mid_mask]):.4f})")
    if np.any(low_mask):
        print(f"  Low score  (< 100 bits, n = {np.sum(low_mask)}):   Mean ratio = {all_r[low_mask].mean():.4f} +/- {all_r[low_mask].std():.4f} (Median = {np.median(all_r[low_mask]):.4f})")

    if sample_out_path and all_pairs:
        # Write 5000 sampled rows for reviewer inspection
        np.random.seed(42)
        indices = np.random.choice(len(all_pairs), size=min(5000, len(all_pairs)), replace=False)
        with open(sample_out_path, "w") as f:
            f.write("chrom\tstart\tend\tname\tsplit\thmmer_score\tfasthumas_scaled_score\traw_gotoh_score\tratio_hmmer_to_raw\n")
            for idx in sorted(indices):
                p = all_pairs[idx]
                f.write(f"{p['chrom']}\t{p['start']}\t{p['end']}\t{p['name']}\t{p['split']}\t{p['hmmer_score']:.2f}\t{p['fasthumas_score']:.2f}\t{p['raw_gotoh_score']:.2f}\t{p['ratio']:.4f}\n")
        print(f"\nWrote 5,000 paired sample scores to {sample_out_path}")

if __name__ == "__main__":
    main()
