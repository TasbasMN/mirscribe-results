# sequence_utils.py
import numpy as np

def get_nucleotides_in_interval(chrom, start, end):
    """Read nucleotides from a FASTA file."""
    chrom = str(chrom)
    file_path = f"data/fasta/grch37/Homo_sapiens.GRCh37.dna.chromosome.{chrom}.fa"
    with open(file_path, 'r') as file:
        file.readline()  # Skip header
        byte_position = file.tell()
        line_length = len(file.readline().strip())
        start_offset = start - 1
        end_offset = end - 1
        num_start_new_lines = start_offset // line_length
        num_end_new_lines = end_offset // line_length
        start_byte_position = byte_position + start_offset + num_start_new_lines
        end_byte_position = byte_position + end_offset + num_end_new_lines
        file.seek(start_byte_position)
        nucleotides = file.read(end_byte_position - start_byte_position + 1).replace('\n', '')
        return nucleotides

def get_nucleotide_at_position(chrom, position):
    """Get a nucleotide at a specific position."""
    file_path = f"data/fasta/grch37/Homo_sapiens.GRCh37.dna.chromosome.{chrom}.fa"
    with open(file_path, 'r') as file:
        file.readline()  # Skip header
        byte_position = file.tell()
        line_length = len(file.readline().strip())
        offset = position - 1
        num_new_lines = offset // line_length
        byte_position = byte_position + offset + num_new_lines
        file.seek(byte_position)
        return file.read(1)
