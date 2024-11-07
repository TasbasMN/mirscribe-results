# main.py
import os
import argparse
from tqdm.auto import tqdm
from scripts.pyensembl_operations import initialize_genome
from scripts.data_processing import process_file

def find_input_files(folder_path: str) -> list:
    """Find all CSV files in the specified folder."""
    return [os.path.join(root, file) for root, _, files in os.walk(folder_path)
            for file in files if file.endswith('.csv') and not file.endswith('case_2.csv')]

def main(target_folder: str) -> None:
    """Main function to process files."""
    genome = initialize_genome(37)
    input_files = find_input_files(target_folder)
    successful_files = 0
    
    progress_bar = tqdm(total=len(input_files), desc="Overall progress")
    for file_path in input_files:
        try:
            df = process_file(file_path, genome)
            # Save the DataFrame to a parquet file
            output_dir = f"results/{os.path.basename(target_folder)}"
            os.makedirs(output_dir, exist_ok=True)
            output_path = os.path.join(output_dir, f"{os.path.splitext(os.path.basename(file_path))[0]}.parquet")
            df.to_parquet(output_path)
            successful_files += 1
            progress_bar.update(1)
        except Exception as e:
            progress_bar.write(f"Failed to process {os.path.basename(file_path)}: {str(e)}")
    progress_bar.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process VCF files.")
    parser.add_argument("target_folder", help="Path to the target folder containing VCF files.")
    args = parser.parse_args()
    main(args.target_folder)
