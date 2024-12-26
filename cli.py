import pandas as pd
import numpy as np
import os
import argparse
from scripts.misc import *

def main():
    parser = argparse.ArgumentParser(description='Process gene regulation counts')
    parser.add_argument('--db_path', type=str, default="results.db",
                        help='Path to the database file')
    parser.add_argument('--threshold', type=float, default=0.35,
                        help='Threshold value')
    parser.add_argument('--additional_filter', type=bool, default=False,
                        help='Whether to apply additional filtering')
    parser.add_argument('--filter_val', type=float, default=0.25,
                        help='Filter value')
    
    args = parser.parse_args()

    real_counts = get_gene_regulation_counts(args.db_path, "real", args.threshold, 
                                           apply_additional_filter=args.additional_filter, 
                                           low_threshold=args.filter_val, 
                                           high_threshold=1-args.filter_val)
    
    synth_counts = get_gene_regulation_counts(args.db_path, "synth", args.threshold, 
                                            apply_additional_filter=args.additional_filter, 
                                            low_threshold=args.filter_val, 
                                            high_threshold=1-args.filter_val)

    counts = pd.merge(real_counts, synth_counts, how="inner", on="gene_id", suffixes=["_real", "_synth"])
    counts["upregulated_synth"] = counts["upregulated_synth"] / 10
    counts["downregulated_synth"] = counts["downregulated_synth"] / 10

    counts['log2_odds_ratio'] = counts.apply(lambda row: calculate_log2_odds_ratio(
        row['upregulated_real'], 
        row['downregulated_real'], 
        row['upregulated_synth'], 
        row['downregulated_synth']
    ), axis=1)

    counts['shrunk_log2_odds'] = shrink_log2_odds(counts)
    counts = add_z_score(counts)
    counts = perform_fisher_test_vectorized(counts, bonf_holm=False)

    counts["is_significant"] = counts['p_value'] < 0.05
    counts["is_significant_adj"] = counts['p_adj'] < 0.05

    export_path = f"results/last/{args.threshold:.2f}".replace("0.", "0")
    os.makedirs(export_path, exist_ok=True)
    count_sign = len(counts[counts.p_value < 0.05])
    count_sign_adj = len(counts[counts.p_adj < 0.05])
    filter_string = f"{args.filter_val:.2f}".replace(".", "")

    if args.additional_filter:
        counts.to_csv(f"{export_path}/counts_{filter_string}filter_sig{count_sign}_adj{count_sign_adj}.csv", index=False)
    else:
        counts.to_csv(f"{export_path}/counts_nofilter_sig{count_sign}_adj{count_sign_adj}.csv", index=False)

if __name__ == "__main__":
    main()
