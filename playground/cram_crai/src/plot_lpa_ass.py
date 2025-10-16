#!/usr/bin/env python3
"""
Plot LPA KIV2 Copy Number vs Lp(a) Levels

This script merges KIV2 copy number estimates with Lp(a) measurements
and creates visualization plots showing the association.
"""

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from pathlib import Path
import sys

# Set plot style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 12


def load_kiv2_estimates(estimate_file):
    """
    Load KIV2 copy number estimates.
    
    Expected columns: ID, exon1A, exon1B, dip_estimate, estimate
    """
    print(f"Loading KIV2 estimates from {estimate_file}...")
    try:
        df = pd.read_csv(estimate_file, sep='\t')
        
        # Check if 'ID' column exists, if so, set it as index
        if 'ID' in df.columns:
            df.set_index('ID', inplace=True)
        
        print(f"  Loaded {len(df)} samples")
        print(f"  Columns: {list(df.columns)}")
        return df
    except Exception as e:
        print(f"ERROR: Failed to load KIV2 estimates: {e}", file=sys.stderr)
        sys.exit(1)


def load_lpa_data(lpa_file):
    """
    Load Lp(a) measurements.
    
    Expected columns: id, subject_id, nwd_id, gwas_id, bkgrd1, lpa_nmol
    """
    print(f"\nLoading Lp(a) data from {lpa_file}...")
    try:
        df = pd.read_csv(lpa_file)
        print(f"  Loaded {len(df)} samples")
        
        # Remove quotes from nwd_id if present
        if 'nwd_id' in df.columns:
            df['nwd_id'] = df['nwd_id'].astype(str).str.strip('"').str.strip()
        
        # Filter out missing NWD IDs and Lp(a) values
        original_len = len(df)
        df = df[df['nwd_id'].notna() & (df['nwd_id'] != 'NA') & (df['nwd_id'] != 'nan')]
        df = df[df['lpa_nmol'].notna()]
        
        print(f"  After filtering: {len(df)} samples with valid NWD ID and Lp(a)")
        print(f"  Removed {original_len - len(df)} samples with missing data")
        
        # Convert lpa_nmol to numeric (in case there are string issues)
        df['lpa_nmol'] = pd.to_numeric(df['lpa_nmol'], errors='coerce')
        df = df[df['lpa_nmol'].notna()]
        
        print(f"  First few nwd_ids: {list(df['nwd_id'].head())}")
        
        return df
    except Exception as e:
        print(f"ERROR: Failed to load Lp(a) data: {e}", file=sys.stderr)
        sys.exit(1)


def merge_data(kiv2_df, lpa_df):
    """
    Merge KIV2 estimates with Lp(a) data on NWD ID.
    """
    print("\nMerging KIV2 estimates with Lp(a) data...")
    
    # KIV2 estimates have NWD ID as index, reset it to a column
    kiv2_df_reset = kiv2_df.reset_index()
    
    # Rename the ID column to nwd_id for consistency
    if 'ID' in kiv2_df_reset.columns:
        kiv2_df_reset.rename(columns={'ID': 'nwd_id'}, inplace=True)
    elif 'index' in kiv2_df_reset.columns:
        kiv2_df_reset.rename(columns={'index': 'nwd_id'}, inplace=True)
    else:
        # First column is the ID column
        kiv2_df_reset.rename(columns={kiv2_df_reset.columns[0]: 'nwd_id'}, inplace=True)
    
    # Strip whitespace from nwd_id in both dataframes
    kiv2_df_reset['nwd_id'] = kiv2_df_reset['nwd_id'].astype(str).str.strip()
    lpa_df['nwd_id'] = lpa_df['nwd_id'].astype(str).str.strip()
    
    print(f"  KIV2 columns: {list(kiv2_df_reset.columns)}")
    print(f"  First few KIV2 IDs: {list(kiv2_df_reset['nwd_id'].head())}")
    print(f"  First few Lp(a) IDs: {list(lpa_df['nwd_id'].head())}")
    
    # Merge on nwd_id
    merged = pd.merge(
        kiv2_df_reset,
        lpa_df[['nwd_id', 'lpa_nmol', 'bkgrd1']],
        on='nwd_id',
        how='inner'
    )
    
    print(f"  Merged dataset: {len(merged)} samples")
    print(f"  Samples in KIV2 only: {len(kiv2_df) - len(merged)}")
    print(f"  Samples in Lp(a) only: {len(lpa_df) - len(merged)}")
    
    if len(merged) == 0:
        print("ERROR: No overlapping samples found!", file=sys.stderr)
        print("\nDebugging - checking for case sensitivity or format differences...")
        kiv2_ids = set(kiv2_df_reset['nwd_id'].str.upper())
        lpa_ids = set(lpa_df['nwd_id'].str.upper())
        overlap = kiv2_ids & lpa_ids
        print(f"  Case-insensitive overlap: {len(overlap)} samples")
        if len(overlap) > 0:
            print("  Issue: Sample IDs have different cases!")
        sys.exit(1)
    
    return merged


def compute_statistics(merged_df):
    """
    Compute correlation and regression statistics.
    """
    print("\nComputing statistics...")
    
    # Remove any remaining NaN values
    clean_df = merged_df[['estimate', 'lpa_nmol']].dropna()
    
    # Pearson correlation
    pearson_r, pearson_p = stats.pearsonr(clean_df['estimate'], clean_df['lpa_nmol'])
    
    # Spearman correlation (non-parametric)
    spearman_r, spearman_p = stats.spearmanr(clean_df['estimate'], clean_df['lpa_nmol'])
    
    # Linear regression
    slope, intercept, r_value, p_value, std_err = stats.linregress(
        clean_df['estimate'], clean_df['lpa_nmol']
    )
    
    print(f"  N = {len(clean_df)}")
    print(f"  Pearson correlation: r = {pearson_r:.3f}, p = {pearson_p:.2e}")
    print(f"  Spearman correlation: ρ = {spearman_r:.3f}, p = {spearman_p:.2e}")
    print(f"  Linear regression: slope = {slope:.3f}, R² = {r_value**2:.3f}")
    
    return {
        'n': len(clean_df),
        'pearson_r': pearson_r,
        'pearson_p': pearson_p,
        'spearman_r': spearman_r,
        'spearman_p': spearman_p,
        'slope': slope,
        'intercept': intercept,
        'r_squared': r_value**2,
        'reg_p': p_value
    }


def create_plots(merged_df, stats_dict, output_prefix):
    """
    Create visualization plots.
    """
    print("\nCreating plots...")
    
    # Remove NaN values for plotting
    plot_df = merged_df[['estimate', 'lpa_nmol', 'bkgrd1']].dropna()
    
    # Figure 1: Scatter plot with regression line
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Scatter plot
    ax.scatter(plot_df['estimate'], plot_df['lpa_nmol'], 
               alpha=0.5, s=15, edgecolors='black', linewidth=0.5)
    
    # Add regression line
    x_range = np.linspace(plot_df['estimate'].min(), plot_df['estimate'].max(), 100)
    y_pred = stats_dict['slope'] * x_range + stats_dict['intercept']
    # ax.plot(x_range, y_pred, 'r-', linewidth=2, label='Linear regression')
    
    # Labels and title
    ax.set_xlabel('KIV2 Copy Number (haploid estimate)', fontsize=14)
    ax.set_ylabel('Lp(a) [nmol/L]', fontsize=14)
    ax.set_title('KIV2 Copy Number and Lp(a) Levels', 
                 fontsize=16, fontweight='bold')
    
    # Add statistics text
    # stats_text = (f"N = {stats_dict['n']}\n"
    #               f"Pearson r = {stats_dict['pearson_r']:.3f}\n"
    #               f"p = {stats_dict['pearson_p']:.2e}\n"
    #               f"R² = {stats_dict['r_squared']:.3f}")
    # ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
    #         fontsize=12, verticalalignment='top',
    #         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Save figure 1
    plot1_path = f"{output_prefix}_scatter.png"
    plt.tight_layout()
    plt.savefig(plot1_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {plot1_path}")
    plt.close()
    
    # Figure 2: Log-scale Lp(a)
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Log transform Lp(a) (add small constant to handle zeros)
    plot_df['log_lpa'] = np.log10(plot_df['lpa_nmol'] + 0.1)
    
    ax.scatter(plot_df['estimate'], plot_df['log_lpa'], 
               alpha=0.5, s=50, edgecolors='black', linewidth=0.5)
    
    # Regression on log scale
    slope_log, intercept_log, r_log, p_log, _ = stats.linregress(
        plot_df['estimate'], plot_df['log_lpa']
    )
    x_range = np.linspace(plot_df['estimate'].min(), plot_df['estimate'].max(), 100)
    y_pred_log = slope_log * x_range + intercept_log
    ax.plot(x_range, y_pred_log, 'r-', linewidth=2, label='Linear regression')
    
    ax.set_xlabel('KIV2 Copy Number (haploid estimate)', fontsize=14)
    ax.set_ylabel('log₁₀(Lp(a) [nmol/L])', fontsize=14)
    ax.set_title('KIV2 Copy Number vs log-transformed Lp(a)', 
                 fontsize=16, fontweight='bold')
    
    stats_text_log = (f"N = {len(plot_df)}\n"
                      f"R² = {r_log**2:.3f}\n"
                      f"p = {p_log:.2e}")
    ax.text(0.05, 0.95, stats_text_log, transform=ax.transAxes,
            fontsize=12, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plot2_path = f"{output_prefix}_scatter_log.png"
    plt.tight_layout()
    plt.savefig(plot2_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {plot2_path}")
    plt.close()
    
    # Figure 3: Distribution plots
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # KIV2 distribution
    axes[0, 0].hist(plot_df['estimate'], bins=30, edgecolor='black', alpha=0.7)
    axes[0, 0].set_xlabel('KIV2 Copy Number')
    axes[0, 0].set_ylabel('Count')
    axes[0, 0].set_title('Distribution of KIV2 Copy Numbers')
    axes[0, 0].axvline(plot_df['estimate'].median(), color='r', 
                       linestyle='--', label=f'Median: {plot_df["estimate"].median():.2f}')
    axes[0, 0].legend()
    
    # Lp(a) distribution
    axes[0, 1].hist(plot_df['lpa_nmol'], bins=30, edgecolor='black', alpha=0.7)
    axes[0, 1].set_xlabel('Lp(a) [nmol/L]')
    axes[0, 1].set_ylabel('Count')
    axes[0, 1].set_title('Distribution of Lp(a) Levels')
    axes[0, 1].axvline(plot_df['lpa_nmol'].median(), color='r',
                       linestyle='--', label=f'Median: {plot_df["lpa_nmol"].median():.2f}')
    axes[0, 1].legend()
    
    # Box plot by KIV2 bins
    plot_df['kiv2_bin'] = pd.cut(plot_df['estimate'], bins=5)
    axes[1, 0].boxplot([group['lpa_nmol'].values for name, group in plot_df.groupby('kiv2_bin', observed=True)],
                       labels=[f"{interval.left:.1f}-{interval.right:.1f}" 
                              for interval in plot_df.groupby('kiv2_bin', observed=True).groups.keys()])
    axes[1, 0].set_xlabel('KIV2 Copy Number Bins')
    axes[1, 0].set_ylabel('Lp(a) [nmol/L]')
    axes[1, 0].set_title('Lp(a) Distribution by KIV2 Bins')
    axes[1, 0].tick_params(axis='x', rotation=45)
    
    # Residual plot
    y_pred_full = stats_dict['slope'] * plot_df['estimate'] + stats_dict['intercept']
    residuals = plot_df['lpa_nmol'] - y_pred_full
    axes[1, 1].scatter(plot_df['estimate'], residuals, alpha=0.5)
    axes[1, 1].axhline(0, color='r', linestyle='--')
    axes[1, 1].set_xlabel('KIV2 Copy Number')
    axes[1, 1].set_ylabel('Residuals')
    axes[1, 1].set_title('Residual Plot')
    axes[1, 1].grid(True, alpha=0.3)
    
    plot3_path = f"{output_prefix}_distributions.png"
    plt.tight_layout()
    plt.savefig(plot3_path, dpi=300, bbox_inches='tight')
    print(f"  Saved: {plot3_path}")
    plt.close()


def save_merged_data(merged_df, output_file):
    """
    Save merged dataset for downstream analysis.
    """
    print(f"\nSaving merged dataset to {output_file}...")
    merged_df.to_csv(output_file, sep='\t', index=False)
    print(f"  Saved {len(merged_df)} samples")


def main():
    parser = argparse.ArgumentParser(
        description='Plot association between KIV2 copy number and Lp(a) levels'
    )
    parser.add_argument(
        '--kiv2_file',
        required=True,
        help='KIV2 copy number estimates file (output from step9)'
    )
    parser.add_argument(
        '--lpa_file',
        required=True,
        help='Lp(a) measurements file (diploid_calls.txt)'
    )
    parser.add_argument(
        '--output_prefix',
        required=True,
        help='Output prefix for plots and merged data'
    )
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_prefix).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Load data
    kiv2_df = load_kiv2_estimates(args.kiv2_file)
    lpa_df = load_lpa_data(args.lpa_file)
    
    # Merge datasets
    merged_df = merge_data(kiv2_df, lpa_df)
    
    # Compute statistics
    stats_dict = compute_statistics(merged_df)
    
    # Create plots
    create_plots(merged_df, stats_dict, args.output_prefix)
    
    # Save merged data
    merged_file = f"{args.output_prefix}_merged_data.tsv"
    save_merged_data(merged_df, merged_file)
    
    print("\n" + "="*50)
    print("Analysis complete!")
    print("="*50)


if __name__ == '__main__':
    main()