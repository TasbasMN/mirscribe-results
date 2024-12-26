import dask.dataframe as dd
import pandas as pd
import duckdb
import sqlite3
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

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



def determine_result(row):
    """
    Determines if gene regulation supports or contradicts its cancer role
    """
    if pd.isna(row['cancer_gene_role']):
        return 'unknown'
    
    if row['cancer_gene_role'] == 'dual_role':
        return 'supports'
    elif row['cancer_gene_role'] == 'oncogene' and row['is_gene_upregulated']:
        return 'supports'
    elif row['cancer_gene_role'] == 'tumor_suppressor' and not row['is_gene_upregulated']:
        return 'supports'
    else:
        return 'contradicts'


def load_and_merge_gene_data(counts_file_path, db_path='data/mirscribe.db'):
    """
    Loads counts data from CSV and merges with gene information from SQLite database.
    
    Args:
        counts_file_path (str): Path to the counts CSV file
        db_path (str): Path to the SQLite database (default: 'data/mirscribe.db')
        
    Returns:
        pd.DataFrame: Merged dataframe with gene information
    """
    # Read counts data
    df1 = pd.read_csv(counts_file_path)
    
    # Connect to SQLite and get gene/mirna data
    with sqlite3.connect(db_path) as sqlite_conn:
        genes = pd.read_sql('SELECT * FROM genes', sqlite_conn)
        mirnas = pd.read_sql('SELECT * FROM mirnas', sqlite_conn)
    
    # Columns to merge from genes table
    cols_to_merge_genes = [
        'gene_id', 'gene_name', 
        'is_oncogene_oncokb', 'is_tsupp_oncokb',
        'is_brca_driver', 'tier_cosmic', 
        'is_hallmark_cosmic',
        'is_tsupp_cosmic', 'is_oncogene_cosmic', 
        'is_oncogene_consensus',
        'is_tsupp_consensus', 'cancer_gene_role'
    ]
    
    # Merge with genes data
    df1 = pd.merge(df1, genes[cols_to_merge_genes], how="left", on="gene_id")
    
    # Add derived columns
    df1["is_gene_upregulated"] = df1["log2_odds_ratio"] > 0
    df1["is_brca_driver"] = df1["is_brca_driver"].astype("bool")
    df1['result'] = df1.apply(determine_result, axis=1)
    
    return df1




def prepare_breast_cancer_data(
    counts_file: str,
    results_db: str,
    genes_db: str,
    threshold: float = 0.40,
    apply_additional_filter: bool = True,
    filter_val: float = 0.25
) -> pd.DataFrame:

    # Load and process gene data
    df = load_and_merge_gene_data(counts_file)
    df["is_oncogene"] = (df.cancer_gene_role == "oncogene") | (df.cancer_gene_role == "dual_role")
    df["is_tsupp"] = (df.cancer_gene_role == "tumor suppressor") | (df.cancer_gene_role == "dual_role")
    
    # Get filtered data
    real = get_filtered_data(
        results_db, 
        "real", 
        threshold=threshold,
        apply_additional_filter=apply_additional_filter,
        low_threshold=filter_val,
        high_threshold=1-filter_val
    )
    real = real[real.mutsig != "SBS0"]
    
    # Define and apply subclass mapping
    subclass_mapping = {
        'pnn': 'luminal a',
        'ppp': 'luminal b',
        'ppn': 'luminal a',
        'nnp': 'her2-enriched',
        'nnn': 'triple negative'
    }
    
    real['cancer_subclass'] = real['cancer_type'].map(subclass_mapping)
    real['cancer_subclass'].fillna("unknown", inplace=True)
    
    # Load and merge miRNA data
    with sqlite3.connect(genes_db) as sqlite_conn:
        genes = pd.read_sql('SELECT * FROM genes', sqlite_conn)
        mirnas = pd.read_sql_query("""
            SELECT mirna_accession, mirna_family, 
                   upregulated_in_brca, downregulated_in_brca 
            FROM mirnas""", sqlite_conn)
    
    # Merge with miRNA data
    real = real.merge(mirnas, how="left", on="mirna_accession")
    
    return real, df


def plot_mutsig_distribution(result_type):
    """
    Plot distribution of gene categories by mutsig for either supporting or contradicting genes
    
    Args:
        result_type (str): Either "supports" or "contradicts"
    """
    if result_type not in ["supports", "contradicts"]:
        raise ValueError("result_type must be either 'supports' or 'contradicts'")
        
    # Create flags for dual-role genes
    merged_df['is_dual'] = merged_df['is_oncogene'] & merged_df['is_tsupp']
    filtered_df = merged_df[merged_df.result == result_type]
    
    # Separate single-role genes
    upreg_oncogenes = filtered_df[(filtered_df.is_oncogene) & (~filtered_df.is_dual) & (filtered_df.log2_odds_ratio > 0)]
    downreg_tsupps = filtered_df[(filtered_df.is_tsupp) & (~filtered_df.is_dual) & (filtered_df.log2_odds_ratio < 0)]
    downreg_oncogenes = filtered_df[(filtered_df.is_oncogene) & (~filtered_df.is_dual) & (filtered_df.log2_odds_ratio < 0)]
    upreg_tsupps = filtered_df[(filtered_df.is_tsupp) & (~filtered_df.is_dual) & (filtered_df.log2_odds_ratio > 0)]

    # Handle dual-role genes separately
    upreg_dual = filtered_df[filtered_df.is_dual & (filtered_df.log2_odds_ratio > 0)]
    downreg_dual = filtered_df[filtered_df.is_dual & (filtered_df.log2_odds_ratio < 0)]

    # Get all unique mutsig values
    all_mutsigs = pd.concat([
        upreg_oncogenes['mutsig'],
        downreg_tsupps['mutsig'],
        upreg_tsupps['mutsig'],
        downreg_oncogenes['mutsig'],
        upreg_dual['mutsig'],
        downreg_dual['mutsig']
    ]).unique()
    all_mutsigs = np.sort(all_mutsigs)

    # Plotting code
    plt.figure(figsize=(12, 6))
    x = np.arange(len(all_mutsigs))
    width = 0.2

    # Stack the greens (harmful changes)
    plt.bar(x, upreg_oncogenes['mutsig'].value_counts().reindex(all_mutsigs).fillna(0), width, 
            label='Upregulated Oncogenes', color='#00ff00')
    plt.bar(x, downreg_tsupps['mutsig'].value_counts().reindex(all_mutsigs).fillna(0), width, 
            label='Downregulated Tumor Suppressors', color='#004d00',
            bottom=upreg_oncogenes['mutsig'].value_counts().reindex(all_mutsigs).fillna(0))

    # Stack the reds (beneficial changes)
    plt.bar(x + width, upreg_tsupps['mutsig'].value_counts().reindex(all_mutsigs).fillna(0), width,
            label='Upregulated Tumor Suppressors', color='#ff0000')
    plt.bar(x + width, downreg_oncogenes['mutsig'].value_counts().reindex(all_mutsigs).fillna(0), width,
            label='Downregulated Oncogenes', color='#800000',
            bottom=upreg_tsupps['mutsig'].value_counts().reindex(all_mutsigs).fillna(0))

    # Add dual-role genes with a different color (e.g., purple)
    plt.bar(x + 2*width, upreg_dual['mutsig'].value_counts().reindex(all_mutsigs).fillna(0), width,
            label='Upregulated Dual-Role', color='#800080')
    plt.bar(x + 2*width, downreg_dual['mutsig'].value_counts().reindex(all_mutsigs).fillna(0), width,
            label='Downregulated Dual-Role', color='#4B0082',
            bottom=upreg_dual['mutsig'].value_counts().reindex(all_mutsigs).fillna(0))

    plt.xlabel('Mutsig')
    plt.ylabel('Count')
    plt.title(f'Distribution of Gene Categories by Mutsig, {result_type.capitalize()} Genes')
    plt.xticks(x + width, all_mutsigs, rotation=45)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.show()
