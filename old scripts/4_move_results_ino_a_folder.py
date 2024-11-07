import os
import shutil

def move_csv_files(source_folder, target_folder):
    # Create the target folder if it doesn't exist
    if not os.path.exists(target_folder):
        os.makedirs(target_folder)

    # Walk through the source folder and its subfolders
    for root, dirs, files in os.walk(source_folder, topdown=False):
        for file in files:
            # Check if the file has a .csv extension
            if file.lower().endswith('.csv'):
                # Construct the full file path
                source_path = os.path.join(root, file)
                
                # Construct the target file path
                target_path = os.path.join(target_folder, file)
                
                # Move the file to the target folder
                shutil.move(source_path, target_path)
                print(f"Moved: {source_path} -> {target_path}")

    # Remove empty directories
    for root, dirs, files in os.walk(source_folder, topdown=False):
        for dir in dirs:
            dir_path = os.path.join(root, dir)
            if not os.listdir(dir_path):
                os.rmdir(dir_path)
                print(f"Removed empty directory: {dir_path}")

# Example usage
source_folder = "downloads/560_real/"
target_folder = "crawled/"

move_csv_files(source_folder, target_folder)
