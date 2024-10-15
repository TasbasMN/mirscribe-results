import argparse
import paramiko
from src.scripts import *

def main():
    parser = argparse.ArgumentParser(description="Traverse remote directory for .err files")
    parser.add_argument("target_folder", nargs='?', default="slurm_logs", help="Target folder relative to home directory")
    parser.add_argument("-l", "--local_path", default="logs", help="Local path to download the files (default: downloads)")
    args = parser.parse_args()
    
    

    TARGET_FOLDER = os.path.join(HOME_FOLDER, args.target_folder)
    TARGET_FOLDER_LAST = os.path.basename(TARGET_FOLDER)
    
    
    try:
        ssh = setup_ssh_connection(HOSTNAME, USERNAME, PASSWORD)
        if ssh is None:
            print("Failed to establish SSH connection.")
            return

        
        successful_jobs = []
        error_jobs = []        
        
        
        with ssh.open_sftp() as sftp:
        
            print(f"Traversing the directory: {TARGET_FOLDER}")
            traverse_directory_for_reports(sftp, TARGET_FOLDER, successful_jobs, error_jobs)
            
            
        
        
        
        
        
        
        
        
        
        
        

        # # Your SFTP operations here
        # with ssh.open_sftp() as sftp:
        #     # Example: List files in the target folder
        #     print(f"Files in {TARGET_FOLDER}:")
        #     for file in sftp.listdir(TARGET_FOLDER):
        #         print(file)

        # # Example: Run a command to find .err files
        # cmd(ssh, f"find {TARGET_FOLDER} -name '*.err'")

        ssh.close()
        print("SSH connection closed.")
    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    main()