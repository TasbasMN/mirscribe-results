from socket import gethostname

def get_pyensembl_cache_location():
    hostname = gethostname()
    if hostname == "nazo":
        return "/home/nazif/thesis/data"
    elif hostname == "Minerva":
        return "/home/yamak/Code/nazif/data"
    else:
        return "/truba/home/mtasbas/data"


PYENSEMBL_CACHE_DIR = get_pyensembl_cache_location()

GRCH37_DIR = "data/fasta/grch37"
MIRNA_COORDS_DIR = "data/mirna_coordinates"
TA_SPS_CSV = "data/ta_sps/ta_sps.csv"
MIRNA_CSV = "data/mirna/mirna.csv"
XGB_MODEL = "misc/models/model_with_no_close_proximity.json"

NUCLEOTIDE_OFFSET = 30


AWK_SCRIPT_PATH = "scripts/rnaduplex_to_csv.awk"
MUTSIG_PROBABILITIES = "data/mutsig_probabilities/probabilities.csv"


from scripts.globals import PYENSEMBL_CACHE_DIR
from functools import lru_cache
import numpy as np

def import_pyensembl(grch):
    if grch not in [37, 38]:
        raise ValueError("grch must be either 37 or 38")

    ens_release = 75 if grch == 37 else 111
    
    import os

    from pyensembl import EnsemblRelease
    os.environ['PYENSEMBL_CACHE_DIR'] = PYENSEMBL_CACHE_DIR
    assembly = EnsemblRelease(ens_release)
    assembly.download()
    assembly.index()
    
    return assembly


def generate_transcript_id_and_gene_name_columns(df, assembly):
    
    df['transcript_id'] = df.apply(lambda x: assembly.transcript_ids_at_locus(x['chr'], x['pos']), axis=1)
    # df["transcript_name"] = df.apply(lambda x: assembly.transcript_names_at_locus(x['chr'], x['pos']), axis=1)
    df["gene_name"] = df.apply(lambda x: assembly.gene_names_at_locus(x['chr'], x['pos']), axis=1)
    
    

@lru_cache(maxsize=None)
def cached_pyensembl_call(locus, assembly, canonical_only, function_name):
    """
    Call a function from the assembly object with the given locus and parameters.

    Args:
        locus (str): A string in the format "chr:pos" representing the chromosome and position.
        assembly (object): An object containing functions to be called.
        canonical_only (bool): If True, return only the first result. If False, return all results.
        function_name (str): The name of the function to call from the assembly object.

    Returns:
        str: The result of the function call, formatted as a string.
            If the result is an empty list, returns np.nan.
            If canonical_only is True, returns the first element of the result list.
            If canonical_only is False, returns a comma-separated string of all elements in the result list.

    Example:
        >>> cached_pyensembl_call('chr1:12345', assembly_obj, True, 'gene_ids_at_locus')
        'ENSGxxxxxxxxx'
        >>> cached_pyensembl_call('chr2:67890', assembly_obj, False, 'transcript_ids_at_locus')
        'ENSTxxxxxxxxx,ENSTxxxxxxxxx,ENSTxxxxxxxxx'
    """
    chrom, pos = locus.split(':')
    
    # Get the function from the assembly object using its name
    func = getattr(assembly, function_name)
    
    result = func(chrom, int(pos))
    
    if len(result) == 0:
        return np.nan
    elif canonical_only:
        return result[0]
    else:
        return ','.join(result)




def get_gene_names(df, assembly, canonical_only=True):
    
    unique_loci = df['locus'].unique()

    gene_name_dict = {
        locus: cached_pyensembl_call(locus, assembly, canonical_only, 'gene_names_at_locus')
        for locus in unique_loci
    }
    
    df['gene_name'] = df['locus'].map(gene_name_dict)

    return df

def get_transcript_ids(df, assembly, canonical_only=False):

    # Get unique locus values
    unique_loci = df['locus'].unique()

    transcript_id_dict = {
        locus: cached_pyensembl_call(locus, assembly, canonical_only, "transcript_ids_at_locus")
        for locus in unique_loci
    }
    
    if canonical_only:
        df['transcript_id'] = df['locus'].map(transcript_id_dict)
    
    else:
        # Create a new DataFrame with expanded transcript IDs
        expanded_df = df['locus'].map(transcript_id_dict).str.split(',', expand=True).stack().reset_index(level=1, drop=True).to_frame('transcript_id')
        
        # Merge the expanded DataFrame with the original DataFrame
        df = df.join(expanded_df)
   

    return df


def get_gene_biotypes(df, assembly):

    # Get unique locus values
    unique_genes = df['gene_id'].unique()

    biotype_dict = {
        gene: 'NA' if gene == 'NA' else assembly.gene_by_id(gene).biotype
        for gene in unique_genes
    }
    # Create a new column 'biotype' by mapping the 'locus' column to the biotype dictionary
    df['biotype'] = df['gene_id'].map(biotype_dict)
    
    return df



from functools import lru_cache
import pandas as pd
from scripts.globals import *


def calculate_au_content(sequence):
    au_count = sequence.count('A') + sequence.count('T') + sequence.count('U')
    total_length = len(sequence)
    return au_count / total_length if total_length > 0 else None


@lru_cache(maxsize=None)
def get_nucleotides_in_interval(chrom, start, end):
    """
    Given a chromosome name, start and end positions, this function reads the DNA sequence from the corresponding FASTA file and returns the nucleotides in the specified interval.

    Parameters:
    - chrom (str): The name of the chromosome.
    - start (int): The starting position of the interval.
    - end (int): The ending position of the interval.

    Returns:
    - nucleotides (str): The nucleotides in the specified interval.
    """

    # change chrom into str
    chrom = str(chrom)

    file_path = f"{GRCH37_DIR}/Homo_sapiens.GRCh37.dna.chromosome.{chrom}.fa"
    with open(file_path, 'r') as file:
        file.readline()
        byte_position = file.tell()
        line_length = len(file.readline().strip())
        start_offset = start - 1
        end_offset = end - 1
        num_start_new_lines = start_offset // line_length
        num_end_new_lines = end_offset // line_length
        start_byte_position = byte_position + start_offset + num_start_new_lines
        end_byte_position = byte_position + end_offset + num_end_new_lines
        file.seek(start_byte_position)

        # Read the nucleotides in the interval
        nucleotides = file.read(end_byte_position - start_byte_position + 1)

    # Remove newlines from the nucleotides
    nucleotides = nucleotides.replace('\n', '')

    return nucleotides


@lru_cache(maxsize=None)
def get_nucleotide_at_position(chrom, position):
    """
    Given a chromosome name and a position, this function reads the DNA sequence from the corresponding FASTA file and returns the nucleotide at the specified position.

    Parameters:
    - chrom (str): The name of the chromosome.
    - position (int): The position of the nucleotide.

    Returns:
    - nucleotide (str): The nucleotide at the specified position.
    """
    file_path = f"{GRCH37_DIR}/Homo_sapiens.GRCh37.dna.chromosome.{chrom}.fa"
    with open(file_path, 'r') as file:
        file.readline()
        byte_position = file.tell()
        line_length = len(file.readline().strip())
        offset = position - 1
        num_new_lines = offset // line_length
        byte_position = byte_position + offset + num_new_lines
        file.seek(byte_position)

        # Read the nucleotide at the position
        nucleotide = file.read(1)
    return nucleotide


@lru_cache(maxsize=None)
def get_upstream_sequence(chrom, pos, n=30):
    """
    Get the upstream sequence of length n from the given position.

    Args:
        row (pandas.Series): A row from the DataFrame containing the 'chr', 'pos', and 'ref_len' columns.
        n (int, optional): The length of the upstream sequence. Defaults to 30.

    Returns:
        str: The upstream sequence.
    """
    int_pos = int(pos)
    upstream_start = max(1, int_pos - n)
    upstream_end = int_pos - 1
    return get_nucleotides_in_interval(chrom, upstream_start, upstream_end)


@lru_cache(maxsize=None)
def get_downstream_sequence(chrom, pos, ref, n=30):
    """
    Get the downstream sequence of length n from the given position.

    Args:
        chrom (str): The chromosome name.
        pos (int): The position.
        ref (str): The reference allele.
        n (int, optional): The length of the downstream sequence. Defaults to 30.

    Returns:
        str: The downstream sequence.
    """
    int_pos = int(pos)
    ref_len = len(ref)
    downstream_start = int_pos + ref_len
    downstream_end = downstream_start + n - 1
    return get_nucleotides_in_interval(chrom, downstream_start, downstream_end)


@lru_cache(maxsize=None)
def get_mre_sequence(mrna_sequence, mrna_end, mirna_start, mirna_length):
    mre_end = mrna_end + mirna_start
    # Ensure MRE start is not negative
    mre_start = max(mre_end - mirna_length, 0)
    return mrna_sequence[mre_start:mre_end]


import pandas as pd
from scripts.sequence_utils import get_nucleotide_at_position
from scripts.pyensembl_operations import *
from scripts.globals import *
import os
import sqlite3
from scripts.globals import MUTSIG_PROBABILITIES
from tqdm.auto import tqdm  # Change this import


def filter_rows_with_same_prediction(df, threshold=0.5):
    mask = (df['wt_prediction'] >= threshold) != (df['mut_prediction'] >= threshold)
    return df[mask].reset_index(drop=True)



def split_id_column(df):
    df[['vcf_id', 'chr', 'pos', 'ref', 'alt', 'mirna_accession']] = df['id'].str.split('_', expand=True)
    df['pos'] = pd.to_numeric(df['pos'], downcast='integer')
    return df



COMPLEMENT = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

def create_mutation_context_string(row):
    ref, alt, before, after = row['ref'], row['alt'], row['before'], row['after']
    if ref in 'CT':  # Pyrimidine
        return f"{before}[{ref}>{alt}]{after}"
    elif ref in 'AG':  # Purine
        ref_complement = COMPLEMENT[ref]
        alt_complement = COMPLEMENT[alt]
        before_complement = ''.join(COMPLEMENT[base] for base in before[::-1])
        after_complement = ''.join(COMPLEMENT[base] for base in after[::-1])
        return f"{after_complement}[{ref_complement}>{alt_complement}]{before_complement}"
    else:
        raise ValueError(f"Invalid reference nucleotide: {ref}")



def generate_mutation_context_column(df):
    df['before'] = df.apply(lambda x: get_nucleotide_at_position(x['chr'], x['pos']-1), axis=1)
    df['after'] = df.apply(lambda x: get_nucleotide_at_position(x['chr'], x['pos']+1), axis=1)
    df['mutation_context'] = df.apply(create_mutation_context_string, axis=1)
    df['mutsig_key'] = df['vcf_id'].astype(str) + '_' + df['mutation_context']
    df.drop(columns=['before', 'after', 'ref', 'alt'], inplace=True)
    return df



def add_mutsig_probabilities(df, mutsig_file):
    mutsig_df = pd.read_csv(mutsig_file)
    mutsig_df['mutsig_key'] = mutsig_df['Sample Names'] + '_' + mutsig_df['MutationTypes']
    mutsig_dict = dict(zip(mutsig_df['mutsig_key'], mutsig_df['mutsig']))
    df['mutsig'] = df['mutsig_key'].map(mutsig_dict)
    return df




def generate_is_intron_column(df, assembly):

    @lru_cache(maxsize=None)
    def cached_pyensembl_call(locus):

        chrom, pos = locus.split(':')
        result = assembly.exons_at_locus(chrom, int(pos))

        return np.nan if len(result) == 0 else result[0].exon_id

    unique_loci = df['locus'].unique()
    map_dict = {locus: cached_pyensembl_call(locus) for locus in unique_loci}
    df['is_exon'] = df['locus'].map(map_dict)
    
    mask = (df.is_exon.isna()) & (~df.gene_id.isna())

    df["is_intron"] = mask

    df.drop(columns=["is_exon"], inplace=True)

    return df



###################################################


def generate_gene_id_column(df, assembly):
    unique_loci = df['locus'].unique()
    
    @lru_cache(maxsize=None)
    def cached_pyensembl_call(locus):

        chrom, pos = locus.split(':')
        result = assembly.gene_ids_at_locus(chrom, int(pos))

        return np.nan if len(result) == 0 else result[0]

    gene_id_dict = {locus: cached_pyensembl_call(locus) for locus in unique_loci}
    df['gene_id'] = df['locus'].map(gene_id_dict)

    return df


def apply_step_5(file_path, assembly, mutsig_probabilities):
    """Process a single VCF file with detailed progress tracking."""
    
    # Define processing steps for progress tracking
    steps = [
        ("Reading CSV", lambda _: pd.read_csv(file_path)),
        ("Optimizing dtypes", optimize_dtypes),
        ("Filtering rows", filter_rows_with_same_prediction),
        ("Splitting ID column", split_id_column),
        ("Converting categories", convert_categories),
        ("Processing gene IDs", lambda df: process_gene_ids(df, assembly)),
        ("Generating contexts", generate_contexts),
        ("Finalizing", finalize_processing)
    ]
    
    df = None
    # Create progress bar for sub-steps
    progress_bar = tqdm(total=len(steps), 
                       desc=f"Processing {os.path.basename(file_path)}", 
                       leave=False, 
                       position=1)
    
    try:
        for step_name, step_func in steps:
            progress_bar.set_description(f"{step_name: <20}")
            try:
                df = step_func(df) if df is not None else step_func(None)
                progress_bar.update(1)
            except Exception as e:
                progress_bar.write(f"Error in {step_name}: {str(e)}")
                raise
    finally:
        progress_bar.close()
    
    return df

def optimize_dtypes(df):
    """Optimize datatypes for predictions."""
    pred_columns = ['wt_prediction', 'mut_prediction']
    for col in pred_columns:
        df[col] = df[col].astype(np.float32)
    return df

def convert_categories(df):
    """Convert columns to category type."""
    cat_columns = ['vcf_id','chr', 'pos', 'ref', 'alt', 'mirna_accession']
    for col in cat_columns:
        df[col] = df[col].astype('category')
    df["locus"] = df["chr"].astype(str) + ":" + df["pos"].astype(str)
    return df

def process_gene_ids(df, assembly):
    """Process gene IDs and introns."""
    # Import the specific functions from step5 module
    from scripts.step5 import generate_gene_id_column, generate_is_intron_column
    
    df = generate_gene_id_column(df, assembly)
    df['gene_id'] = df['gene_id'].astype('category')
    df = generate_is_intron_column(df, assembly)
    df.drop(columns=["locus"], inplace=True)
    return df

def generate_contexts(df):
    """Generate mutation contexts."""
    # Import the specific function from step5 module
    from scripts.step5 import generate_mutation_context_column
    
    df = generate_mutation_context_column(df)
    df.drop(columns=["chr", "pos"], inplace=True)
    df["mutation_context"] = df["mutation_context"].astype('category')
    return df

def finalize_processing(df):
    """Perform final processing steps."""
    df["is_gene_upregulated"] = (df.wt_prediction > df.mut_prediction).astype('bool')
    
    for col in ["gene_id"]:
        existing_categories = df[col].cat.categories
        new_categories = existing_categories.append(pd.Index(['not_found']))
        df[col] = df[col].cat.set_categories(new_categories)
        df[col] = df[col].fillna('not_found')
    
    return df



def crawl_and_import_results(folder_path, ending_string, db_path, table_name, assembly):
    csv_files = []

    # Find CSV files
    for root, _, files in os.walk(folder_path):
        csv_files.extend(
            os.path.join(root, file)
            for file in files
            if file.endswith(f"{ending_string}.csv")
        )

    # Connect to SQLite database
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    # Create table if it doesn't exist
    cursor.execute(f"""
        CREATE TABLE IF NOT EXISTS {table_name} (
            id TEXT PRIMARY KEY,
            wt_prediction REAL,
            mut_prediction REAL,
            pred_difference REAL,
            vcf_id TEXT,
            mirna_accession TEXT,
            gene_id TEXT,
            mutation_context TEXT,
            mutsig TEXT,
            is_intron BOOLEAN,
            is_gain BOOLEAN,
            is_gene_upregulated BOOLEAN
        )
    """)

    # Import CSV files into the table
    for file_path in csv_files:
        
        df = apply_step_5(file_path, assembly, MUTSIG_PROBABILITIES)

        columns = ', '.join(df.columns)
        placeholders = ', '.join(['?'] * len(df.columns))

        # Insert data into the table
        for row in df.values:
            cursor.execute(f"INSERT INTO {table_name} ({columns}) VALUES ({placeholders})", row)

    # Commit changes and close the connection
    conn.commit()
    conn.close()



def crawl_and_import_results_into_df(folder_path, ending_string, assembly):
    csv_files = []

    # Find CSV files
    for root, _, files in os.walk(folder_path):
        csv_files.extend(
            os.path.join(root, file)
            for file in files
            if file.endswith(f"{ending_string}.csv")
        )

    # Create a list to store DataFrames
    dataframes = []

    # Process each CSV file
    for file_path in csv_files:
        df = apply_step_5(file_path, assembly, MUTSIG_PROBABILITIES)
        dataframes.append(df)

    return pd.concat(dataframes, ignore_index=True)




import os
import sys
import argparse
import time
from typing import List, Tuple, Optional

import pandas as pd
from tqdm.auto import tqdm
from scripts import step5 as s5
from scripts.globals import *


from scripts.step5 import apply_step_5


def find_csv_files(folder_path: str) -> Tuple[List[str], List[str]]:
    """Finds case 1 and case 2 CSV files in the given folder and subfolders.

    Args:
        folder_path: Path to the root folder to search.

    Returns:
        A tuple containing two lists:
        - List of paths to case 1 CSV files
        - List of paths to case 2 CSV files
    """
    case_1s = []
    case_2s = []
    
    for root, _, files in os.walk(folder_path):
        for file in files:
            file_path = os.path.join(root, file)
            if file.endswith(".csv"):
                case_1s.append(file_path)
            elif file.endswith("case_2.csv"):
                case_2s.append(file_path)
    
    return case_1s, case_2s


def save_single_result(df: pd.DataFrame, file_path: str, target_folder: str) -> None:
    """Saves a single processed DataFrame to a parquet file.

    Args:
        df: DataFrame to save.
        file_path: Original CSV file path.
        target_folder: Target folder for output.
    """
    # Create output directory if it doesn't exist
    output_dir = os.path.join("results", target_folder)
    os.makedirs(output_dir, exist_ok=True)
    
    # Generate output filename from input filename
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    output_path = os.path.join(output_dir, f"{base_name}.parquet")
    
    # Save the DataFrame
    df.to_parquet(output_path)



def process_single_file(file_path: str) -> Optional[pd.DataFrame]:
    """Processes a single CSV file and returns the resulting DataFrame.

    Args:
        file_path: Path to the CSV file to process.

    Returns:
        Processed DataFrame or None if processing fails.
    """
    try:
        # Call apply_step_5 with the correct parameters
        processed_df = s5.apply_step_5(
            file_path=file_path,
            assembly="37",
            mutsig_probabilities=MUTSIG_PROBABILITIES
        )
        return processed_df
    except Exception as e:
        tqdm.write(f"Error processing {os.path.basename(file_path)}: {str(e)}")
        return None



def main(target_folder: str) -> None:
    """Main function to process VCF files individually."""
    folder_name = os.path.basename(target_folder)
    case_1s, _ = find_csv_files(target_folder)
    
    # Import and initialize pyensembl
    from scripts.step5 import import_pyensembl
    genome = import_pyensembl(37)
    if not genome:
        print("Failed to initialize pyensembl")
        return
    
    # Process and save files individually with progress bar
    successful_files = 0
    progress_bar = tqdm(total=len(case_1s), 
                       desc="Overall progress", 
                       position=0, 
                       leave=True)
    
    try:
        for file_path in case_1s:
            try:
                # Process single file with genome object
                df = s5.apply_step_5(
                    file_path=file_path,
                    assembly=genome,
                    mutsig_probabilities=MUTSIG_PROBABILITIES
                )
                
                if df is not None:
                    save_single_result(df, file_path, folder_name)
                    successful_files += 1
                
                progress_bar.update(1)
                progress_bar.set_postfix(
                    processed=f"{successful_files}/{len(case_1s)}",
                    success_rate=f"{successful_files/len(case_1s)*100:.1f}%"
                )
                
            except Exception as e:
                progress_bar.write(f"Failed to process {os.path.basename(file_path)}: {str(e)}")
                continue
    finally:
        progress_bar.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Process VCF files and generate a combined DataFrame.")
    parser.add_argument(
        "target_folder",
        help="Path to the target folder containing VCF files")
    
    args = parser.parse_args()
    main(args.target_folder)
    os.system("paplay /usr/share/sounds/gnome/default/alerts/drip.ogg")
