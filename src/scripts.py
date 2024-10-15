import os
import paramiko
from paramiko import SSHException
from src.globals import *
from stat import S_ISDIR


def create_ssh_client():
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    username, ip = HOSTNAME.split('@')
    ssh.connect(ip, username=username, password=PASSWORD)
    return ssh


def setup_ssh_connection(hostname, username, password):
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        print("Connecting to the remote server...")
        ssh.connect(hostname=hostname, username=username, password=password)
        print("Connection established.")
        return ssh
    except SSHException as e:
        print(f"Failed to establish SSH connection: {e}")
        return None

def cmd(ssh, command, print_output=True):
    stdin, stdout, stderr = ssh.exec_command(command)
    output = stdout.read().decode().strip()
    error = stderr.read().decode().strip()
    
    print(f"$ {command}")
    if error:
        print(f"Error: {error}")
    if print_output:
        print(output)
    return output



def traverse_directory_for_reports(sftp, remote_path, successful_files, error_files):
    try:
        for entry in sftp.listdir_attr(remote_path):
            mode = entry.st_mode
            filename = entry.filename
            filepath = os.path.join(remote_path, filename)
            
            if S_ISDIR(mode):
                traverse_directory_for_reports(sftp, filepath, successful_files, error_files)
            elif filename.endswith('.err'):
                if entry.st_size == 0:
                    successful_files.append(filepath)
                else:
                    error_info = f"{filepath} (Size: {entry.st_size} bytes)"
                    error_files.append(error_info)
    except IOError as e:
        print(f"Could not access {remote_path}: {e}")

    print(f"Success count: {len(successful_files)}")
    print(f"Error count: {len(error_files)}")
