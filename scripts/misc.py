import dask.dataframe as dd
import pandas as pd
import duckdb


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