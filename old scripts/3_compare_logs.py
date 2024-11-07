import argparse
import re

def extract_pd_names(file_path):
    pd_names = set()
    with open(file_path, 'r') as file:
        for line in file:
            match = re.search(r'(PD\d+a)', line.strip())
            if match:
                pd_names.add(match.group(1))
    return pd_names

def compare_pd_files(file1_path, file2_path):
    pd_names1 = extract_pd_names(file1_path)
    pd_names2 = extract_pd_names(file2_path)

    common_pds = pd_names1.intersection(pd_names2)
    only_in_file1 = pd_names1 - pd_names2
    only_in_file2 = pd_names2 - pd_names1

    print(f"Total PDs in {file1_path}: {len(pd_names1)}")
    print(f"Total PDs in {file2_path}: {len(pd_names2)}")
    print(f"PDs present in both files: {len(common_pds)}")
    print(f"PDs only in {file1_path}: {len(only_in_file1)}")
    print(f"PDs only in {file2_path}: {len(only_in_file2)}")

    if only_in_file1:
        print(f"\nPDs only in {file1_path}:")
        for pd in sorted(only_in_file1):
            print(pd)

    if only_in_file2:
        print(f"\nPDs only in {file2_path}:")
        for pd in sorted(only_in_file2):
            print(pd)

def main():
    parser = argparse.ArgumentParser(description="Compare PD names in two text files")
    parser.add_argument("file1", help="Path to the first text file")
    parser.add_argument("file2", help="Path to the second text file")
    args = parser.parse_args()

    compare_pd_files(args.file1, args.file2)

if __name__ == "__main__":
    main()
