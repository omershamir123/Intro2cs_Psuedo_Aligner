# FILE : test_file_handlers.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: This file tests the file handlers module
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
import os
import tempfile
from typing import List, Callable, Generator

import pytest

import file_handlers
from genome import ReferencedGenome
from kmer_reference import KmerReference
from program_constants import FASTQ_FILE_TYPES, FASTA_FILE_TYPES, \
    ALN_FILE_TYPES, KDB_FILE_TYPES
from pseudo_aligner import AlnFileDataObject, PseudoAlignerOutput
from read import Read

EMPTY_FILE = b"\n\n\n\n\n\n\n\n\r\t"
FASTQ_FILE_NOT_FOUND_PATH = "PATH_NOT_FOUND_UXWX.fq"
FASTA_FILE_NOT_FOUND_PATH = "PATH_NOT_FOUND_UXWX.fa"

VALID_FASTQ_FILE = b"@NAME\nATCGC\n+\nABCDE"
NO_HEADER_FASTQ_FILE = b"NAME"
NO_READ_FASTQ_FILE = b"@NAME"
NO_PLUS_LINE_FASTQ_FILE = b"@NAME\nATCGC\n"
NO_QUALITY_SEQUENCE_FASTQ_FILE = b"@NAME\nATCGC\n+\n"

VALID_FASTA_FILE = b">GENOME1\nATCGNCG\n>GENOME2\nATCGNNNN"
NO_GENOME_FASTA_FILE = b"@>NAME\nATCGC\n+\n"
NO_SEQUENCE_FASTA_FILE = b">GENOME1\n"
NO_SEQUENCE_SECOND_GENOME_FASTA_FILE = b">GENOME1\nATCGNCG\n>GENOME2\n\r\n"


def create_temporary_file(content: bytes, wanted_suffix: str) -> str:
    """
    This function creates a temporary file with the content given.
    It returns the filename of the temporary file.
    Note - the temporary file is not deleted and the user must remove it by calling os.remove().
    :param wanted_suffix:  the suffix of the temporary file.
    :param content: the content of the temporary file.
    :return: the filename of the temporary file (full path)
    """
    with tempfile.NamedTemporaryFile(delete=False,
                                     suffix=wanted_suffix) as temporary_file:
        temporary_file.write(content)
        file_name = temporary_file.name
    return file_name


def check_outputs_in_file(file_path: str,
                          expected_outputs: List[object], generator_from_file:Callable[[str], Generator]) -> None:
    outputs_from_file = []
    try:
        for genome in generator_from_file(file_path):
            outputs_from_file.append(genome)
    except ValueError:
        pass
    os.remove(file_path)
    assert len(outputs_from_file) == len(expected_outputs) and all(
        expected_outputs[i].__dict__ == outputs_from_file[i].__dict__ for i in
        range(len(outputs_from_file)))

def test_parse_fastq_file_valid():
    expected_read_output = [Read("NAME", "ATCGC", "ABCDE")]
    fastq_file_path = create_temporary_file(VALID_FASTQ_FILE,
                                            FASTQ_FILE_TYPES[0])
    check_outputs_in_file(fastq_file_path, expected_read_output, file_handlers.parse_fastq_file)


def test_parse_fastq_file_invalid():
    fastq_file_path = create_temporary_file(NO_HEADER_FASTQ_FILE,
                                            FASTQ_FILE_TYPES[0])
    check_outputs_in_file(fastq_file_path, [], file_handlers.parse_fastq_file)

    fastq_file_path = create_temporary_file(NO_READ_FASTQ_FILE,
                                            FASTQ_FILE_TYPES[0])
    check_outputs_in_file(fastq_file_path, [], file_handlers.parse_fastq_file)

    fastq_file_path = create_temporary_file(NO_PLUS_LINE_FASTQ_FILE,
                                            FASTQ_FILE_TYPES[0])
    check_outputs_in_file(fastq_file_path, [], file_handlers.parse_fastq_file)

    fastq_file_path = create_temporary_file(NO_QUALITY_SEQUENCE_FASTQ_FILE,
                                            FASTQ_FILE_TYPES[0])
    check_outputs_in_file(fastq_file_path, [], file_handlers.parse_fastq_file)

    fastq_file_path = create_temporary_file(EMPTY_FILE, FASTQ_FILE_TYPES[0])
    check_outputs_in_file(fastq_file_path, [], file_handlers.parse_fastq_file)

    with pytest.raises(FileNotFoundError):
        for read in file_handlers.parse_fastq_file(FASTQ_FILE_NOT_FOUND_PATH):
            print(read)


def test_parse_fasta_file_valid():
    expected_genome_output = [ReferencedGenome("GENOME1", "ATCGNCG", 0),
                              ReferencedGenome("GENOME2", "ATCGNNNN", 1)]
    fasta_file_path = create_temporary_file(VALID_FASTA_FILE,
                                            FASTA_FILE_TYPES[0])
    check_outputs_in_file(fasta_file_path, expected_genome_output,
                          file_handlers.parse_fasta_file)


def test_parse_fasta_file_invalid():
    fasta_file_path = create_temporary_file(NO_GENOME_FASTA_FILE,
                                            FASTA_FILE_TYPES[0])
    check_outputs_in_file(fasta_file_path, [], file_handlers.parse_fasta_file)

    fasta_file_path = create_temporary_file(NO_GENOME_FASTA_FILE,
                                            FASTA_FILE_TYPES[0])
    check_outputs_in_file(fasta_file_path, [], file_handlers.parse_fasta_file)

    expected_genomes = [ReferencedGenome("GENOME1", "ATCGNCG", 0)]
    fasta_file_path = create_temporary_file(NO_SEQUENCE_SECOND_GENOME_FASTA_FILE,
                                            FASTA_FILE_TYPES[0])
    check_outputs_in_file(fasta_file_path, expected_genomes, file_handlers.parse_fasta_file)

    fasta_file_path = create_temporary_file(EMPTY_FILE,
                                            FASTA_FILE_TYPES[0])
    check_outputs_in_file(fasta_file_path, [], file_handlers.parse_fasta_file)

    with pytest.raises(FileNotFoundError):
        for genome in file_handlers.parse_fasta_file(FASTA_FILE_NOT_FOUND_PATH):
            print(genome)


def test_write_and_read_pickle_file():
    reference = KmerReference(31)
    align_output = AlnFileDataObject(PseudoAlignerOutput(reference))
    aln_file_path = create_temporary_file(b"",ALN_FILE_TYPES[0])
    reference_file_path = create_temporary_file(b"",KDB_FILE_TYPES[0])
    assert file_handlers.write_to_pickle_file(aln_file_path, align_output, ALN_FILE_TYPES) == True
    assert file_handlers.write_to_pickle_file(reference_file_path, reference, KDB_FILE_TYPES) == True

    assert file_handlers.decompress_pickle_file(aln_file_path, ALN_FILE_TYPES).__dict__ == align_output.__dict__
    assert file_handlers.decompress_pickle_file(reference_file_path, KDB_FILE_TYPES).__dict__ == reference.__dict__

    os.remove(aln_file_path)
    os.remove(reference_file_path)

    unpickled_file = create_temporary_file(VALID_FASTQ_FILE, KDB_FILE_TYPES[0])
    unpickled_object = file_handlers.decompress_pickle_file(unpickled_file, KDB_FILE_TYPES[0])
    assert unpickled_object is None
