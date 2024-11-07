"""Module for processing large batches of VCF files with progress tracking."""

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
