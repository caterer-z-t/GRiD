#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import argparse

def main(input_file, output_prefix):
    # Load table
    df = pd.read_csv(input_file, sep="\t")
    
    # Drop rows with NA in Percent_Diff
    df = df[df["Percent_Diff"] != "NA"].copy()
    
    # Convert to numeric
    df["Percent_Diff"] = pd.to_numeric(df["Percent_Diff"])
    
    # Plot histogram
    plt.figure(figsize=(8,6))
    plt.hist(df["Percent_Diff"], bins=30, edgecolor="black", alpha=0.7)
    plt.xlabel("Percent Difference (%)")
    plt.ylabel("Number of Samples")
    plt.title("Distribution of Percent Differences")
    plt.grid(axis="y", linestyle="--", alpha=0.6)
    
    plt.savefig(f"{output_prefix}_percent_diff_hist.png", dpi=300, bbox_inches="tight")
    plt.close()
    
    print(f"✅ Plot saved to {output_prefix}_percent_diff_hist.png")
    print(f"📊 Mean: {df['Percent_Diff'].mean():.4f}, Median: {df['Percent_Diff'].median():.4f}, Max: {df['Percent_Diff'].max():.4f}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot Percent_Diff distribution from VNTR read counts file.")
    parser.add_argument("input_file", help="Path to input .txt file")
    parser.add_argument("-o", "--output", default="percent_diff", help="Output prefix for plots")
    args = parser.parse_args()
    main(args.input_file, args.output)
