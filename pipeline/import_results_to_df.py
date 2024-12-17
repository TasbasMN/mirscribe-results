import os
import sys
import argparse
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
        processed_df = s5.apply_step_5(
            file_path=file_path,
            assembly="37",
            mutsig_probabilities=MUTSIG_PROBABILITIES
        )
        return processed_df
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

def main(target_folder, mutsig_file, clinpath_file, experiment_value):
    folder_name = os.path.basename(target_folder)
    case_1s, _ = find_csv_files(target_folder)
    
    genome = import_pyensembl(37)
    if not genome:
        print("Failed to initialize pyensembl")
        return
    
    successful_files = 0
    progress_bar = tqdm(total=len(case_1s), desc="Overall progress", position=0, leave=True)
    
    try:
        for file_path in case_1s:
            df = process_single_file(file_path)
            if df is not None:
                save_single_result(df, file_path, folder_name)
                successful_files += 1
            progress_bar.update(1)
            progress_bar.set_postfix(processed=f"{successful_files}/{len(case_1s)}", success_rate=f"{successful_files/len(case_1s)*100:.1f}%")
    finally:
        progress_bar.close()

    # Now process Parquet results
    input_folder = os.path.join("results", folder_name)
    output_folder = f'results/processed/{experiment_value}'
    process_parquet_files(input_folder, output_folder, mutsig_file, clinpath_file, experiment_value)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process VCF files and generate a combined DataFrame.")
    parser.add_argument("target_folder", help="Path to the target folder containing VCF files")
    parser.add_argument("mutsig_file", help="Path to the Mutsig CSV file")
    parser.add_argument("clinpath_file", help="Path to the ClinPath CSV file")
    parser.add_argument("--experiment_value", type=int, default=1, help="Value to add to the 'experiment' column")
    
    args = parser.parse_args()
    main(args.target_folder, args.mutsig_file, args.clinpath_file, args.experiment_value)
    os.system("paplay /usr/share/sounds/gnome/default/alerts/drip.ogg")
