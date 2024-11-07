import os
import argparse
import paramiko
from paramiko import SSHException
from stat import S_ISDIR

def get_password():
    password = os.environ.get('TRUBA_PWD')
    if not password:
        raise ValueError("TRUBA_PWD environment variable not set")
    return password

def traverse_directory(sftp, remote_path, successful_files, error_files):
    try:
        for entry in sftp.listdir_attr(remote_path):
            mode = entry.st_mode
            filename = entry.filename
            filepath = os.path.join(remote_path, filename)
            
            if S_ISDIR(mode):
                traverse_directory(sftp, filepath, successful_files, error_files)
            elif filename.endswith('.err'):
                if entry.st_size == 0:
                    successful_files.append(filepath)
                else:
                    error_info = f"{filepath} (Size: {entry.st_size} bytes)"
                    error_files.append(error_info)
    except IOError as e:
        print(f"Could not access {remote_path}: {e}")


def extract_filename(file_path):
    return os.path.splitext(os.path.basename(file_path))[0].split('_')[-1]

def write_files(file_list, output_path):
    with open(output_path, 'w') as output_file:
        filenames = map(extract_filename, file_list)
        output_file.write('\n'.join(filenames))
    print(f"Successfully wrote {len(file_list)} entries to {output_path}")

def main():
    parser = argparse.ArgumentParser(description="Traverse remote directory for .err files")
    parser.add_argument("target_folder", nargs='?', default="slurm_logs", help="Target folder relative to home directory")
    parser.add_argument("-l", "--local_path", default="logs", help="Local path to download the files (default: downloads)")
    args = parser.parse_args()

    HOSTNAME = "172.16.7.1"
    USERNAME = "mtasbas"
    HOME_FOLDER = "/truba/home/mtasbas/"
    TARGET_FOLDER = os.path.join(HOME_FOLDER, args.target_folder)
    TARGET_FOLDER_LAST = os.path.basename(TARGET_FOLDER)
    
    password = get_password()
    
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    
    try:
        print("Connecting to the remote server...")
        ssh.connect(hostname=HOSTNAME, username=USERNAME, password=password)
        sftp = ssh.open_sftp()
        print("Connection established.")
        
        successful_files = []
        error_files = []
        
        print(f"Traversing the directory: {TARGET_FOLDER}")
        traverse_directory(sftp, TARGET_FOLDER, successful_files, error_files)
        
        # Ensure the local directory exists
        os.makedirs(args.local_path, exist_ok=True)
        
        success_path = os.path.join(args.local_path, f"{TARGET_FOLDER_LAST}_run_successful.txt")
        error_path = os.path.join(args.local_path, f"{TARGET_FOLDER_LAST}_error_files.txt")
        
        # with open(success_path, 'w') as success_file:
        #     for file in successful_files:
        #         filename_to_write = file.split("/")[-1]
        #         filename_to_write = filename_to_write.split(".err")[0]
        #         filename_to_write = filename_to_write.split("_")[-1] 
        #         success_file.write(f"{filename_to_write}\n")
        # print(f"Successfully wrote {len(successful_files)} entries to {success_path}")
        
        write_files(successful_files, success_path)
        write_files(error_files, error_path)
        
        # with open(error_path, 'w') as error_file:
        #     for error in error_files:
        #         error_file.write(f"{error}\n")
        # print(f"Successfully wrote {len(error_files)} entries to {error_path}")
        
    except SSHException as e:
        print(f"SSH connection failed: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        if 'sftp' in locals():
            sftp.close()
        ssh.close()
        print("Connection closed.")

if __name__ == "__main__":
    main()
