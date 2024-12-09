import os
import pandas as pd
from tqdm.auto import tqdm
from scripts import step5 as s5
from scripts.globals import *
from scripts.pyensembl_operations import import_pyensembl


def load_mutsig_dict(mutsig_file):
    mutsig_df = pd.read_csv(mutsig_file)
    mutsig_df['mutsig_key'] = mutsig_df['Sample Names'] + '_' + mutsig_df['MutationTypes']
    return dict(zip(mutsig_df['mutsig_key'], mutsig_df['mutsig']))

def load_clinpath_dict(clinpath_file):
    clinpath_df = pd.read_csv(clinpath_file, usecols=["vcf_id", "initials"])
    return dict(zip(clinpath_df["vcf_id"], clinpath_df["initials"]))

def find_csv_files(folder_path):
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

def save_single_result(df, file_path, target_folder):
    output_dir = os.path.join("results", target_folder)
    os.makedirs(output_dir, exist_ok=True)
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    output_path = os.path.join(output_dir, f"{base_name}.parquet")
    df.to_parquet(output_path)

def process_single_file(file_path):
    try:
        return s5.apply_step_5(
            file_path=file_path,
            assembly="37",
            mutsig_probabilities=MUTSIG_PROBABILITIES,
        )
    except Exception:
        tqdm.write(f"Error processing {os.path.basename(file_path)}")
        return None

def process_parquet_files(input_folder, output_folder, mutsig_file, clinpath_file, experiment_value=1, genome_version=37):
    os.makedirs(output_folder, exist_ok=True)
    g37 = import_pyensembl(genome_version)
    mutsig_dict = load_mutsig_dict(mutsig_file)
    clinpath_dict = load_clinpath_dict(clinpath_file)

    gene_names = {gene.gene_id: gene.gene_name for gene in g37.genes()}
    gene_biotypes = {gene.gene_id: gene.biotype for gene in g37.genes()}

    parquet_files = [f for f in os.listdir(input_folder) if f.endswith('.parquet')]
    
    for filename in tqdm(parquet_files, desc="Processing files", unit="file"):
        filepath = os.path.join(input_folder, filename)
        try:
            df = pd.read_parquet(filepath)
        except Exception:
            continue  # Skip to the next file
        
        df['experiment'] = experiment_value
        df['mutsig'] = df['mutsig_key'].map(mutsig_dict)
        df['gene_name'] = df['gene_id'].map(gene_names).fillna('not_found')
        df['biotype'] = df['gene_id'].map(gene_biotypes).fillna('not_found')
        df['cancer_type'] = df['vcf_id'].map(clinpath_dict)

        output_filepath = os.path.join(output_folder, f'processed_{filename}')
        try:
            df.to_parquet(output_filepath, engine='pyarrow', compression='snappy', index=False)
        except Exception:
            continue  # Skip saving if it fails