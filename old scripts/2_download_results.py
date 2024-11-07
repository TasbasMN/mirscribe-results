import os
import argparse
import paramiko
from tqdm import tqdm
import logging
import stat
import re

logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')

HOSTNAME = "mtasbas@172.16.7.1"
HOME_FOLDER = "/truba/home/mtasbas/"

def get_password():
    password = os.environ.get('TRUBA_PWD')
    if not password:
        raise ValueError("TRUBA_PWD environment variable not set")
    return password

def create_ssh_client():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    username, ip = HOSTNAME.split('@')
    ssh.connect(ip, username=username, password=get_password())
    return ssh

def download_file(sftp, remote_file, local_file):
    try:
        remote_file_size = sftp.stat(remote_file).st_size
        
        if os.path.exists(local_file) and os.path.getsize(local_file) == remote_file_size:
            logging.info(f"File {local_file} already exists and has the correct size. Skipping download.")
            return True

        os.makedirs(os.path.dirname(local_file), exist_ok=True)
        with tqdm(total=remote_file_size, unit='B', unit_scale=True, desc=os.path.basename(remote_file)) as pbar:
            def callback(bytes_transferred, total_bytes):
                pbar.update(bytes_transferred - pbar.n)
            sftp.get(remote_file, local_file, callback=callback)
        
        if os.path.getsize(local_file) == remote_file_size:
            logging.info(f"Successfully downloaded {remote_file}")
            return True
        else:
            logging.warning(f"Downloaded file size mismatch for {remote_file}. Retrying...")
            os.remove(local_file)
            return download_file(sftp, remote_file, local_file)
    except Exception as e:
        logging.error(f"Failed to download {remote_file}: {str(e)}")
        return False

def download_folder(sftp, remote_path, local_path, successful_downloads):
    try:
        items = sftp.listdir_attr(remote_path)
        for item in items:
            remote_item_path = os.path.join(remote_path, item.filename)
            local_item_path = os.path.join(local_path, item.filename)
            
            if stat.S_ISDIR(item.st_mode):
                logging.info(f"Entering directory: {remote_item_path}")
                download_folder(sftp, remote_item_path, local_item_path, successful_downloads)
            elif item.filename.endswith('.csv'):
                if download_file(sftp, remote_item_path, local_item_path):
                    successful_downloads.append(local_item_path)
    except IOError as e:
        logging.error(f"Error accessing {remote_path}: {str(e)}")

def parse_arguments():
    parser = argparse.ArgumentParser(description="Download CSV files from remote HPC")
    parser.add_argument("remote_path", help="Remote path to folder containing CSV files (relative to home folder)")
    parser.add_argument("-l", "--local_path", default="downloads",
                        help="Local path to download the files (default: downloads)")
    return parser.parse_args()

def extract_pd_name(filepath):
    filename = os.path.basename(filepath)
    match = re.search(r'(PD\d+a)', filename)
    if match:
        return match.group(1)
    return filename  # Return the full filename if no match is found

def write_files(file_list, output_path):
    with open(output_path, 'w') as output_file:
        pd_names = map(extract_pd_name, file_list)
        output_file.write('\n'.join(pd_names))
    print(f"Successfully wrote {len(file_list)} entries to {output_path}")

def main():
    args = parse_arguments()
    remote_path = os.path.join(HOME_FOLDER, args.remote_path.lstrip('/'))
    
    # Extract the PD folder name from the remote_path
    pd_folder = os.path.basename(args.remote_path)
    local_path = os.path.join(args.local_path, pd_folder)

    os.makedirs(local_path, exist_ok=True)

    successful_downloads = []

    try:
        ssh = create_ssh_client()

        with ssh.open_sftp() as sftp:
            download_folder(sftp, remote_path, local_path, successful_downloads)

        ssh.close()

        # Use the modified write_files function
        output_file = os.path.join(local_path, "download_successful.txt")
        write_files(successful_downloads, output_file)

        logging.info(f"All CSV files downloaded successfully to {local_path}")
    except Exception as e:
        logging.error(f"An error occurred: {str(e)}", exc_info=True)

if __name__ == "__main__":
    main()
