# FILE : file_handlers.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs FinalProject 2025
# DESCRIPTION: This module contains all file parsing in the project - FASTA, FASTQ, KDB, ALN
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
import gzip
import pickle
from pickle import UnpicklingError
from typing import Generator, Optional, List, Union

import program_constants
from genome import ReferencedGenome
from kmer_reference import KmerReference
from pseudo_aligner import AlnFileDataObject
from read import Read
from validators import validate_file_type


def parse_fastq_file(fastq_file_path: str) -> Generator[Read, None, None]:
    """
    This function handles the parsing of the fastq file.
    First line - header starts with @
    Second line - Read
    Third line - +
    Fourth line - Quality sequence
    :param fastq_file_path: the path to the fastq file
    :return: Generator of each Read in the fastq file
    """
    try:
        validate_file_type(fastq_file_path, program_constants.FASTQ_FILE_TYPES)
    except TypeError as e:
        print(e)
        raise
    try:
        with open(fastq_file_path, 'r') as fastq_file:
            header = fastq_file.readline()
            while True:
                if not header or not header.startswith("@"):
                    raise ValueError("Fastq file does not start with '@' line")
                header = str.strip(header, '\n')[1:]
                read_value = fastq_file.readline()
                if not read_value:
                    raise ValueError(
                        "No sequence in Fastq file for the read {}".format(
                            header))
                read_value = str.rstrip(read_value, "\n")
                if not fastq_file.readline().startswith("+"):
                    raise ValueError(
                        "Fastq file does not contain a third + line")
                quality = fastq_file.readline()
                if not quality:
                    raise ValueError(
                        "No quality sequence in Fastq file for the read {}".format(
                            header))
                quality = str.rstrip(quality, "\n")
                yield Read(header, read_value, quality)
                header = fastq_file.readline()
                if not header:
                    break

    except ValueError:
        raise
    except FileNotFoundError:
        raise FileNotFoundError(f'Fastq file not found: {fastq_file_path}')
    except PermissionError:
        raise PermissionError(f'Permission denied: {fastq_file_path}')
    except IOError:
        raise IOError(f'I/O error: {fastq_file_path}')
    except Exception as e:
        raise Exception(f'Unexpected error: {fastq_file_path}')


def parse_fasta_file(fasta_file_path: str) -> Generator[
    ReferencedGenome, None, None]:
    """
    This function handles the parsing of the fasta file - header line starts with >
    and all subsequent lines contain DNA samples
    :param fasta_file_path: the path to the fasta file
    :return: Generator of each reference in the fasta file (reference - single genome)
    """
    try:
        validate_file_type(fasta_file_path, program_constants.FASTA_FILE_TYPES)
    except TypeError as e:
        print(e)
        raise
    try:
        with open(fasta_file_path, 'r') as fasta_file:
            current_line = fasta_file.readline()
            if not current_line or not current_line.startswith(">"):
                raise ValueError("Fasta file does not start with '>' line")
            current_genome_index = 0
            while True:
                current_line = str.rstrip(current_line, "\n")
                current_genome = [current_line[1:], ""]
                current_line = fasta_file.readline()
                while current_line and not current_line.startswith(">"):
                    current_line = "".join(current_line.split())
                    current_genome[1] += current_line
                    current_line = fasta_file.readline()
                if current_genome[1] == "":
                    raise ValueError(
                        "No genome data found for the header line {}".format(
                            current_genome[0]))
                yield ReferencedGenome(current_genome[0], current_genome[1],
                                       current_genome_index)
                current_genome_index += 1
                if not current_line:
                    break

    except ValueError:
        raise
    except FileNotFoundError:
        raise FileNotFoundError(f'Fasta file not found: {fasta_file_path}')
    except PermissionError:
        print(f'Permission denied: {fasta_file_path}')
        raise
    except IOError:
        raise IOError(f'I/O error: {fasta_file_path}')
    except Exception as e:
        raise Exception(f'Unexpected error: {fasta_file_path}')


def write_to_pickle_file(pickle_file_path: str,
                      object_to_write: object, supported_file_types:List[str]) -> bool:
    """
    This function handles the writing of a python object to a pickle file
    The python object will be pickled and its content stored in a gzip file
    :param supported_file_types: the file types that are supported for this object
    :param pickle_file_path: the path to the pickle file
    :param object_to_write: the  object to write
    :return: True if the file was written, False otherwise
    """
    try:
        validate_file_type(pickle_file_path, supported_file_types)
    except TypeError as e:
        print(e)
        return False
    pickled_object = pickle.dumps(object_to_write)
    try:
        with open(pickle_file_path, 'wb') as kdb_file:
            kdb_file.write(gzip.compress(pickled_object))
    except (FileNotFoundError, PermissionError, IOError) as e:
        print("An error occurred while writing to the given file {}.".format(pickle_file_path))
        return False
    except Exception as e:
        print(
            "An unexpected error occurred while writing to the given file {}.".format(pickle_file_path))
        return False
    return True


def decompress_pickle_file(pickle_file_path: str, supported_file_types:List[str]) -> Union[KmerReference, AlnFileDataObject, None]:
    """
    This function receives a path to a pickled file, checks its integrity and decompresses the file
    as well as unpickles it
    :param supported_file_types: the file types that are supported for this object
    :param pickle_file_path: the path to the kdb file
    :return: the object within the pickle_file
    """
    try:
        validate_file_type(pickle_file_path, supported_file_types)
    except TypeError as e:
        print(e)
        return None
    try:
        with open(pickle_file_path, 'rb') as kdb_file:
            pickled_object = gzip.decompress(kdb_file.read())
    except (FileNotFoundError, PermissionError, IOError) as e:
        print("An error occurred while reading the given file {}.".format(pickle_file_path))
        return None
    except (ValueError, OSError) as e:
        print("There was an error while decompressing the file {}".format(pickle_file_path))
        return None
    except Exception as e:
        print(
            "An unexpected error occurred while reading the given file.".format(pickle_file_path))
        return None

    try:
        returned_object = pickle.loads(pickled_object)
    except (UnpicklingError, TypeError) as e:
        print(
            "There was a problem unpickling the file's decompressed content. File {}".format(pickle_file_path))
        print(e)
        return None
    except Exception as e:
        print(
            "There was an unexpected problem unpickling the file's decompressed content. File {}".format(pickle_file_path))
        return None
    return returned_object
