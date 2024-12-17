import os
import pandas as pd
from tqdm.auto import tqdm
from scripts.pyensembl_operations import import_pyensembl
from scripts.globals import *
from scripts.step5 import *


def find_csv_files(folder_path):
    """Find all CSV files in the given folder."""
    csv_files = []
    for root, _, files in os.walk(folder_path):
        csv_files.extend(
            os.path.join(root, f) for f in files 
            if f.endswith('.csv')
        )
    return csv_files

def load_reference_data(genome, mutsig_file, clinpath_file):
    """Load all reference data including gene mappings."""
    gene_names = {g.gene_id: g.gene_name for g in genome.genes()}
    
    mutsig_df = pd.read_csv(mutsig_file)
    mutsig_df['mutsig_key'] = mutsig_df['Sample Names'] + '_' + mutsig_df['MutationType']
    mutsig_dict = dict(zip(mutsig_df['mutsig_key'], mutsig_df['mutsig']))
    
    clinpath_df = pd.read_csv(clinpath_file, usecols=["vcf_id", "initials"])
    clinpath_dict = dict(zip(clinpath_df["vcf_id"], clinpath_df["initials"]))
    
    return {
        'gene_names': gene_names,
        'mutsig': mutsig_dict,
        'clinpath': clinpath_dict
    }

def process_single_file(file_path, genome, ref_data, output_dir):
    """Process a single file with all steps in one flow."""
    try:
        # Read and basic preprocessing
        df = pd.read_csv(file_path)
        df = optimize_dtypes(df)
        df = filter_rows_with_same_prediction(df)
        
        # ID and category processing
        df = split_id_column(df)
        df = convert_categories(df)
        df = process_gene_ids(df, genome)
        df = generate_contexts(df)
        df = finalize_processing(df)
        
        # Add reference data columns
        df['mutsig'] = df['mutsig_key'].map(ref_data['mutsig'])
        df['gene_name'] = df['gene_id'].map(ref_data['gene_names']).fillna('not_found')
        df['cancer_type'] = df['vcf_id'].map(ref_data['clinpath'])
        
        # Save processed file
        os.makedirs(output_dir, exist_ok=True)
        base_name = os.path.splitext(os.path.basename(file_path))[0]
        output_path = os.path.join(output_dir, f"{base_name}_processed.csv")
        df.to_csv(output_path, index=False)
        return True
        
    except Exception as e:
        print(f"Error processing {os.path.basename(file_path)}: {str(e)}")
        return False

def main(target_folder, output_dir, genome_version=37):
    # Initialize
    genome = import_pyensembl(genome_version)
    ref_data = load_reference_data(genome, PROBABILITIES_LATEST, CLINPATH_FILE)
    
    # Get files
    csv_files = find_csv_files(target_folder)
    
    # Process files with progress bar
    successful = 0
    with tqdm(total=len(csv_files), desc="Processing files") as pbar:
        for file_path in csv_files:
            if process_single_file(file_path, genome, ref_data, output_dir):
                successful += 1
            pbar.update(1)
            pbar.set_postfix(success_rate=f"{successful/len(csv_files)*100:.1f}%")

if __name__ == "__main__":
    TARGET_FOLDER = "downloads/real"
    OUTPUT_DIR = os.path.join('results', os.path.basename(TARGET_FOLDER))
    main(TARGET_FOLDER, OUTPUT_DIR)
