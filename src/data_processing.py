# data_processing.py
import pandas as pd
from tqdm.auto import tqdm
from scripts.pyensembl_operations import initialize_genome, cached_pyensembl_call
from scripts.sequence_utils import get_nucleotide_at_position
from scripts.globals import MUTSIG_PROBABILITIES

def load_and_validate_data(file_path: str) -> pd.DataFrame:
    """Load and validate data from CSV."""
    df = pd.read_csv(file_path)
    required_columns = ['id', 'wt_prediction', 'mut_prediction']
    if not all(col in df.columns for col in required_columns):
        raise ValueError("Missing required columns")
    return df

def process_file(file_path: str, assembly) -> pd.DataFrame:
    """Process a single VCF file."""
    df = load_and_validate_data(file_path)
    df = split_id_column(df)
    df = optimize_dtypes(df)
    df = generate_mutation_context_column(df)
    df = add_mutsig_probabilities(df, MUTSIG_PROBABILITIES)
    df = generate_gene_id_column(df, assembly)
    df = generate_is_intron_column(df, assembly)
    return df

def split_id_column(df):
    df[['vcf_id', 'chr', 'pos', 'ref', 'alt', 'mirna_accession']] = df['id'].str.split('_', expand=True)
    df['pos'] = pd.to_numeric(df['pos'], downcast='integer')
    return df

def optimize_dtypes(df):
    """Optimize DataFrame memory usage."""
    float_cols = ['wt_prediction', 'mut_prediction']
    df[float_cols] = df[float_cols].astype(np.float32)
    cat_cols = ['vcf_id', 'chr', 'pos', 'ref', 'alt', 'mirna_accession']
    df[cat_cols] = df[cat_cols].astype('category')
    return df

def generate_mutation_context_column(df):
    """Generate mutation context information."""
    df['before'] = df.apply(lambda x: get_nucleotide_at_position(x['chr'], x['pos']-1), axis=1)
    df['after'] = df.apply(lambda x: get_nucleotide_at_position(x['chr'], x['pos']+1), axis=1)
    
    def create_mutation_context_string(row):
        ref, alt, before, after = row['ref'], row['alt'], row['before'], row['after']
        return f"{before}[{ref}>{alt}]{after}"  # Simplified for example
    df['mutation_context'] = df.apply(create_mutation_context_string, axis=1)
    df['mutsig_key'] = df['vcf_id'].astype(str) + '_' + df['mutation_context']
    df.drop(columns=['before', 'after', 'ref', 'alt'], inplace=True)
    return df

def add_mutsig_probabilities(df, mutsig_file):
    """Add mutational significance probabilities to the DataFrame."""
    mutsig_df = pd.read_csv(mutsig_file)
    mutsig_df['mutsig_key'] = mutsig_df['Sample Names'] + '_' + mutsig_df['MutationTypes']
    mutsig_dict = dict(zip(mutsig_df['mutsig_key'], mutsig_df['mutsig']))
    df['mutsig'] = df['mutsig_key'].map(mutsig_dict)
    return df

def generate_is_intron_column(df, assembly):
    """Generate a column to indicate if a gene is in an intron."""
    unique_loci = df['locus'].unique()
    mask = {locus: cached_pyensembl_call(locus, assembly, True, 'exons_at_locus') for locus in unique_loci}
    df['is_intron'] = df['locus'].map(lambda x: np.nan if mask[x] else True)
    return df

def generate_gene_id_column(df, assembly):
    """Generate gene ID column from assembly data."""
    unique_loci = df['locus'].unique()
    gene_id_dict = {locus: cached_pyensembl_call(locus, assembly, True, 'gene_ids_at_locus') for locus in unique_loci}
    df['gene_id'] = df['locus'].map(gene_id_dict)
    return df
