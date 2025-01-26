# FILE : file_parsers.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs FinalProject 2025
# DESCRIPTION: This module contains all file parsing in the project - FASTA, FASTQ, KDB, ALN
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
from typing import Generator

from genome import RawGenome


def parse_fasta_file(fasta_file_path: str) -> Generator[RawGenome, None, None]:
    """
    This function handles the parsing of the fasta file - header line starts with >
    and all subsequent lines contain DNA samples
    :param fasta_file_path: the path to the fasta file
    :return: Generator of each reference in the fasta file (reference - single genome)
    """
    try:
        with open(fasta_file_path, 'r') as fasta_file:
            current_line = fasta_file.readline()
            if not current_line or not current_line.startswith(">"):
                raise ValueError("Fasta file does not start with '>' line")
            # TODO check if file is empty without any genome
            while True:
                current_genome = [current_line[1:], ""]
                current_genome[0] = current_line.rstrip("\n")
                current_line = fasta_file.readline()
                while current_line and not current_line.startswith(">"):
                    current_line = "".join(current_line.split())
                    current_genome[1] += current_line
                    current_line = fasta_file.readline()
                yield RawGenome(current_genome[0], current_genome[1])
                if not current_line:
                    break

    except ValueError:
        print(f'file not in fasta format: {fasta_file_path}')
        raise
    except FileNotFoundError:
        print(f'Fasta file not found: {fasta_file_path}')
        raise
    except PermissionError:
        print(f'Permission denied: {fasta_file_path}')
        raise
    except IOError:
        print(f'I/O error: {fasta_file_path}')
        raise
    except:
        print(f'Unexpected error: {fasta_file_path}')
        raise
