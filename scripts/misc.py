import dask.dataframe as dd
import pandas as pd
import duckdb
import sqlite3
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import numpy as np

def load_and_filter_data(real_path, synth_path):

    # Import data
    real = pd.read_parquet(real_path)
    synth = dd.read_parquet(synth_path)
    
    # Apply filters
    filters = [
        ("vcf_id", "PD4120a"),
        ("gene_id", "not_found")
    ]
    
    for column, value in filters:
        real = real[real[column] != value]
        synth = synth[synth[column] != value]
    
    return real, synth

def find_common_genes(pdf, ddf, column='gene_id'):
    # Convert dask column to pandas series first
    dask_genes = set(ddf[column].compute().unique())
    pandas_genes = set(pdf[column].unique())
    
    return pandas_genes.intersection(dask_genes)

def pivot_counts(df, is_dask=False):
    """Pivot a given DataFrame with upreg counts."""
    if is_dask:
        grouped = df.groupby(['gene_id', 'is_gene_upregulated']).size().reset_index().compute()
    else:
        grouped = df.groupby(['gene_id', 'is_gene_upregulated']).size().reset_index()

    result = (grouped.pivot(index='gene_id', 
                            columns='is_gene_upregulated', 
                            values=0)
              .fillna(0)
              .reset_index())

    result.columns = ['gene_id', 'downreg_count', 'upreg_count']

    return result


def get_filtered_data(db_path, table_name, threshold=0.5, apply_additional_filter=False, 
                     low_threshold=0.1, high_threshold=0.9):
    base_query = f"""
    SELECT 
        id,
        wt_prediction,
        mut_prediction,
        pred_difference,
        vcf_id,
        mirna_accession,
        gene_id,
        is_gene_upregulated,
        mutsig,
        gene_name,
        cancer_type
    FROM {table_name}
    WHERE gene_id != 'not_found'
    AND vcf_id != 'PD4120a'
    AND ABS(wt_prediction - mut_prediction) > {threshold}
    AND is_intron = FALSE
    """

    
    if apply_additional_filter:
        additional_filter = f"""
        AND (
            (pred_difference < 0 AND mut_prediction < {low_threshold}) OR
            (pred_difference > 0 AND mut_prediction > {high_threshold})
        )
        """
        query = base_query + additional_filter
    else:
        query = base_query

    dtype_dict = {
        'id': 'int32',
        'wt_prediction': 'float32',
        'mut_prediction': 'float32',
        'pred_difference': 'float32',
        'vcf_id': 'category',
        'mirna_accession': 'category',
        'gene_id': 'category',
        'is_gene_upregulated': 'bool',
        'mutsig': 'category',
        'gene_name': 'category',
        'cancer_type': 'category'
    }

    with sqlite3.connect(db_path) as conn:
        df = pd.read_sql_query(
            query, 
            conn,
            dtype=dtype_dict
        )

    return df

def get_gene_regulation_counts(db_path, table_name, threshold=0.5, apply_additional_filter=False,
                             low_threshold=0.1, high_threshold=0.9):
    
    base_query = f"""
    SELECT 
        gene_id,
        COUNT(*) FILTER (WHERE is_gene_upregulated = TRUE) as upregulated,
        COUNT(*) FILTER (WHERE is_gene_upregulated = FALSE) as downregulated
    FROM {table_name}
    WHERE vcf_id != 'PD4120a' 
    AND gene_id != 'not_found'
    AND ABS(wt_prediction - mut_prediction) > {threshold}
    AND is_intron = FALSE
    """
    
    if apply_additional_filter:
        additional_filter = f"""
        AND (
            (pred_difference < 0 AND mut_prediction < {low_threshold}) OR
            (pred_difference > 0 AND mut_prediction > {high_threshold})
        )
        """
        query = base_query + additional_filter + " GROUP BY gene_id"
    else:
        query = base_query + " GROUP BY gene_id"
    
    with sqlite3.connect(db_path) as conn:
        counts_df = pd.read_sql_query(query, conn)
    
    print(f"Total unique genes: {len(counts_df)}")
    return counts_df

def calculate_log2_odds_ratio(a, b, c, d, k=0.5):

    odds_ratio = ((a + k) * (d + k)) / ((b + k) * (c + k))
    return np.log2(odds_ratio)

def shrink_log2_odds(values, prior_scale=1.0, min_count=10):
    """Implement empirical Bayes shrinkage for log2 odds ratios."""
    total_counts = (values['upregulated_real'] + values['downregulated_real'] +
                   values['upregulated_synth'] + values['downregulated_synth'])
    
    raw_log2_odds = values.apply(lambda row: calculate_log2_odds_ratio(
        row['upregulated_real'], 
        row['downregulated_real'], 
        row['upregulated_synth'], 
        row['downregulated_synth'],
        k=0.5
    ), axis=1)

    weights = 1 - np.exp(-total_counts / min_count)
    prior_mean = np.average(raw_log2_odds, weights=weights)
    prior_var = np.var(raw_log2_odds)
    posterior_scale = prior_scale / (1 + weights * prior_scale)
    
    return weights * raw_log2_odds + (1 - weights) * prior_mean

def perform_fisher_test_vectorized(df, pseudocount=0.01, bonf_holm=False):
    # Add pseudocount to the table
    table = np.array([
        [df['upregulated_real'] + pseudocount, df['downregulated_real'] + pseudocount],
        [df['upregulated_synth'] + pseudocount, df['downregulated_synth'] + pseudocount]
    ]).transpose((2, 0, 1))  # reshape for 2x2 tables

    p_values = np.zeros(len(df))

    for i in range(len(df)):
        _, p_values[i] = fisher_exact(table[i])

    df['p_value'] = p_values
    
    if bonf_holm:
        df['p_adj'] = multipletests(p_values, method='holm')[1]
    
    else:
        df['p_adj'] = multipletests(p_values, method='fdr_bh')[1]
    
    return df

def add_z_score(df):
    # Calculate mean and standard deviation of log2 odds ratios
    mean_log2or = df['log2_odds_ratio'].mean()
    std_log2or = df['log2_odds_ratio'].std()
    
    # Calculate Z-score
    df['z_score'] = (df['log2_odds_ratio'] - mean_log2or) / std_log2or
    
    return df

