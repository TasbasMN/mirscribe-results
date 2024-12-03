import pandas as pd
import dask.dataframe as dd
import sqlite3

def load_and_merge_data():
    # Original dtypes
    dtypes = {
        'id': 'string',
        'wt_prediction': 'float32',    # changed from float16
        'mut_prediction': 'float32',   # changed from float16
        'pred_difference': 'float32',  # changed from float16
        'vcf_id': 'category',
        'mirna_accession': 'category',
        'gene_id': 'category',
        'is_intron': 'bool',
        'mutation_context': 'category', 
        'mutsig_key': 'category',
        'is_gene_upregulated': 'bool',
        'mutsig': 'category',
        'gene_name': 'category',
        'biotype': 'category',
        'cancer_type': 'category'
    }

    # Load the combined CSV
    comb = "results/combined.csv"
    real = dd.read_csv(comb, dtype=dtypes)

    real = real.drop(columns=["experiment", "mutsig_key"])

    db_path = "data/mirscribe_backup.db"
    cols_to_merge_genes = ['gene_id', 'is_oncogene_consensus', 'is_tsupp_consensus']
    cols_to_merge_mirna = ['mirna_accession', 'mirna_family', 'is_mirna_oncogene', 'is_mirna_tumor_suppressor']

    # Load gene and miRNA data from SQLite
    with sqlite3.connect(db_path) as conn:
        genes = pd.read_sql_query("SELECT gene_id, is_oncogene_consensus, is_tsupp_consensus FROM genes", conn)
        mirnas = pd.read_sql_query("SELECT mirna_accession, mirna_family, is_mirna_oncogene, is_mirna_tumor_suppressor FROM mirnas", conn)

    # Convert to Dask DataFrames
    genes_dd = dd.from_pandas(genes[cols_to_merge_genes], npartitions=1)
    genes_dd['gene_id'] = genes_dd['gene_id'].astype('category')

    mirnas_dd = dd.from_pandas(mirnas[cols_to_merge_mirna], npartitions=1)
    mirnas_dd['mirna_accession'] = mirnas_dd['mirna_accession'].astype('category')

    # Merge data
    real = real.merge(genes_dd, how="left", on="gene_id")
    real = real.merge(mirnas_dd, how="left", on="mirna_accession")

    # Fix dtypes
    real.is_oncogene_consensus = real.is_oncogene_consensus.astype('bool')
    real.is_tsupp_consensus = real.is_tsupp_consensus.astype('bool')
    real.is_mirna_oncogene = real.is_mirna_oncogene.astype('bool')
    real.is_mirna_tumor_suppressor = real.is_mirna_tumor_suppressor.astype('bool')
    real.mirna_family = real.mirna_family.astype('category')

    # Load final data CSV
    df = pd.read_csv("results/the_final_data.csv")

    # Add filter column
    filter_logic = (abs(df['log2_odds_ratio']) > 0.32) & (df['p_adj'] < 0.05)
    df["is_significant"] = filter_logic
    df['is_gene_upregulated'] = df['log2_odds_ratio'] >= 0

    # Add genes table details
    sqlite_conn = sqlite3.connect('data/mirscribe_backup.db')
    genes = pd.read_sql('SELECT * FROM genes', sqlite_conn)
    sqlite_conn.close()

    # Define columns to merge
    cols_to_merge = ['gene_id', 'is_oncogene_oncokb', 'is_tsupp_oncokb',
                     'is_driver_intogen', 'tier_cosmic', 'is_hallmark_cosmic',
                     'is_tsupp_cosmic', 'is_oncogene_cosmic', 'is_oncogene_consensus',
                     'is_tsupp_consensus', 'is_gene_of_interest', 'cancer_gene_role']

    # Merge with final DataFrame
    df = pd.merge(df, genes[cols_to_merge], how="left", on="gene_id")

    return real, df

