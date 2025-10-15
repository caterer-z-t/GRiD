#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

def main():
    parser = argparse.ArgumentParser(description="Plot diploid CN vs Lp(a) (nmol/L)")
    parser.add_argument("--diploid_file", required=True,
                        help="Path to diploid CN estimate file (two columns: nwd_id, diploid_cn)")
    parser.add_argument("--pheno_file", required=True,
                        help="Path to phenotype file containing lpa_nmol and nwd_id (or subject_id)")
    parser.add_argument("--output_prefix", required=True,
                        help="Prefix for output files (e.g. results/diploid_lpa)")
    args = parser.parse_args()

    # === READ INPUTS ===
    diploid_df = pd.read_csv(args.diploid_file, sep=r"\s+", header=None, names=["nwd_id", "diploid_cn"])
    pheno_df = pd.read_csv(args.pheno_file, sep=None, engine="python")

    print("\n=== Diploid file preview ===")
    print(diploid_df.head())
    print("Diploid columns:", diploid_df.columns.tolist())

    print("\n=== Pheno file preview ===")
    print(pheno_df.head())
    print("Pheno columns:", pheno_df.columns.tolist())

    # === CLEAN COLUMN NAMES ===
    diploid_df.columns = diploid_df.columns.astype(str).str.lower()
    pheno_df.columns = pheno_df.columns.astype(str).str.lower()

    # === CLEAN ID FORMATTING ===
    for df in [diploid_df, pheno_df]:
        for col in ["nwd_id", "subject_id"]:
            if col in df.columns:
                df[col] = df[col].astype(str).str.strip().str.upper()

    print("\n=== After cleaning ===")
    print(f"Diploid rows: {len(diploid_df)} | Pheno rows: {len(pheno_df)}")

    # === CHECK WHICH COLUMN MATCHES BEST ===
    overlap_nwd = len(set(diploid_df["nwd_id"]) & set(pheno_df.get("nwd_id", [])))
    overlap_subject = len(set(diploid_df["nwd_id"]) & set(pheno_df.get("subject_id", [])))

    print(f"Overlap with pheno.nwd_id: {overlap_nwd}")
    print(f"Overlap with pheno.subject_id: {overlap_subject}")

    # === PICK THE RIGHT MERGE KEY ===
    if overlap_subject > overlap_nwd:
        merge_key = ("nwd_id", "subject_id")
        print("→ Using diploid.nwd_id ↔ pheno.subject_id for merge")
    else:
        merge_key = ("nwd_id", "nwd_id")
        print("→ Using diploid.nwd_id ↔ pheno.nwd_id for merge")

    # === MERGE ===
    merged = pd.merge(
        diploid_df, pheno_df,
        left_on=merge_key[0], right_on=merge_key[1],
        how="inner"
    )
    print(f"Merged {len(merged)} samples")

    # === PLOT ===
    sns.set(style="whitegrid")
    plt.figure(figsize=(8,6))
    sns.scatterplot(
        data=merged,
        x="diploid_cn",
        y="lpa_nmol",
        s=60
    )
    plt.xlabel("Diploid CN Estimate")
    plt.ylabel("Lp(a) (nmol/L)")
    plt.title("Relationship between VNTR Diploid CN and Lp(a) concentration")
    plt.tight_layout()
    plt.savefig(f"{args.output_prefix}_scatter.png", dpi=300)
    plt.close()

    # === CORRELATION ===
    corr = merged["diploid_cn"].corr(merged["lpa_nmol"], method="spearman")
    print(f"Spearman correlation between diploid CN and Lp(a): {corr:.3f}")
    with open(f"{args.output_prefix}_correlation.txt", "w") as f:
        f.write(f"Spearman correlation: {corr:.3f}\n")

if __name__ == "__main__":
    main()
