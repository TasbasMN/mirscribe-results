#!/usr/bin/env python3

import sys
import subprocess

def generate_and_run_commands(PDx):
    commands = [
        f"python 1_crawl_reports.py slurm_logs/{PDx}",
        f"python 2_download_results.py mirscribe-vcf/results/{PDx}",
        f"python 3_compare_logs.py downloads/{PDx}/download_successful.txt logs/{PDx}_run_successful.txt"
    ]

    for cmd in commands:
        print(f"Executing: {cmd}")
        try:
            # Run the command without capturing output
            process = subprocess.Popen(cmd, shell=True)
            process.wait()

            if process.returncode == 0:
                print("Command executed successfully")
            else:
                print("Command failed")
        except Exception as e:
            print(f"An error occurred: {str(e)}")

        print("---")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py PDx")
        sys.exit(1)

    PDx = sys.argv[1]
    generate_and_run_commands(PDx)
