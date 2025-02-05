# FILE : FILE_NAME.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final project 2025
# DESCRIPTION: 
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
import json
import os

from kmer_reference import extract_kmers_from_string, KmerReference
from program_constants import FASTA_FILE_TYPES
from test_file_handlers import create_temporary_file

FULL_FASTA_FILE = (b">Mouse\n"
                   b"AAAAAAAA\n\n\n"
                   b">Duck\n"
                   b"ATCG\n\n"
                   b">Otter\n"
                   b"TTTTTTTTNNNNNNNNNNNNN\n\n"
                   b">Turtle\n"
                   b"TTTTTTTTAAAA\n\n"
                   b">Moose\n"
                   b"NACACACACACAC\n"
                   b">Virus\n"
                   b"AAT")

FULL_FASTA_EXPECTED_KMER_DB = {
    "AAAA": {"Mouse": [0, 1, 2, 3, 4], "Turtle": [8]},
    "ATCG": {"Duck": [0]},
    "TTTT": {"Otter": [0, 1, 2, 3, 4],
             "Turtle": [0, 1, 2, 3, 4]},
    "TTTA": {"Turtle": [5]},
    "TTAA": {"Turtle": [6]},
    "TAAA": {"Turtle": [7]},
    "ACAC": {"Moose": [1, 3, 5, 7, 9]},
    "CACA": {"Moose": [2, 4, 6, 8]}
    }
FULL_FASTA_EXPECTED_GENOME_DICT = {
    "Mouse": {"total_bases": 8, "unique_kmers": 0, "multi_mapping_kmers": 5},
    "Duck": {"total_bases": 4, "unique_kmers": 1, "multi_mapping_kmers": 0},
    "Otter": {"total_bases": 21, "unique_kmers": 0, "multi_mapping_kmers": 5},
    "Turtle": {"total_bases": 12, "unique_kmers": 3, "multi_mapping_kmers": 6},
    "Moose": {"total_bases": 13, "unique_kmers": 9, "multi_mapping_kmers": 0},
    "Virus": {"total_bases": 3, "unique_kmers": 0, "multi_mapping_kmers": 0}}

FULL_FASTA_EXPECTED_OUTPUT = {"Kmers": FULL_FASTA_EXPECTED_KMER_DB,
                              "Summary": FULL_FASTA_EXPECTED_GENOME_DICT}

SIMILARITY_OUTPUT = {"Similarity":
                         {"Mouse":
                              {"kept": "no", "unique_kmers": 0,
                               "total_kmers": 1,
                               "genome_length": 8, "similar_to": "Turtle",
                               "similarity_score": 1.0},
                          "Duck": {"kept": "yes", "unique_kmers": 1,
                                   "total_kmers": 1,
                                   "genome_length": 4, "similar_to": "NA",
                                   "similarity_score": "NA"},
                          "Otter": {"kept": "no", "unique_kmers": 0,
                                    "total_kmers": 1,
                                    "genome_length": 21, "similar_to": "Turtle",
                                    "similarity_score": 1.0},
                          "Turtle": {"kept": "yes", "unique_kmers": 3,
                                     "total_kmers": 5,
                                     "genome_length": 12, "similar_to": "NA",
                                     "similarity_score": "NA"},
                          "Moose": {"kept": "yes", "unique_kmers": 2,
                                    "total_kmers": 2,
                                    "genome_length": 13, "similar_to": "NA",
                                    "similarity_score": "NA"},
                          "Virus": {"kept": "yes", "unique_kmers": 0,
                                    "total_kmers": 0,
                                    "genome_length": 3, "similar_to": "NA",
                                    "similarity_score": "NA"}}}

DUPLICATE_GENOMES_FASTA_FILE = (b">Mouse\n"
                                b"AAAAAAAA\n\n\n"
                                b">Mouse\n"
                                b"AAAAAAAATAAA\n\n\n")

KMER_SIZE = 4

def build_testing_reference()->KmerReference:
    """
    This function is used to build the testing KmerReference used in most tetss
    Its structure is seen above
    :return:
    """
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    reference = KmerReference(KMER_SIZE)
    build_successful = reference.build_kmer_reference(genome_fasta_file)
    reference.calculate_kmers_type()
    os.remove(genome_fasta_file)

    assert build_successful

    return reference


def test_extract_kmers():
    """
    This function checks the kmer extraction functionality.
    It checks sequences with one kmer, no kmers, with and without wildcards
    :return:
    """
    sequence_without_wildcard = "ATCGATCGATCG"

    kmer_size = 8
    expected_kmers = ["ATCGATCG", "TCGATCGA", "CGATCGAT", "GATCGATC",
                      "ATCGATCG"]
    assert [kmer_tuple[0] for kmer_tuple in
            extract_kmers_from_string(sequence_without_wildcard, kmer_size,
                                      True)] == expected_kmers
    assert [kmer_tuple[0] for kmer_tuple in
            extract_kmers_from_string(sequence_without_wildcard, kmer_size,
                                      False)] == expected_kmers

    sequence_with_wildcard = "ATCGATCGATCGNNN"
    expected_kmers_with_wildcard = ["ATCGATCG", "TCGATCGA", "CGATCGAT",
                                    "GATCGATC", "ATCGATCG", "TCGATCGN",
                                    "CGATCGNN", "GATCGNNN"]

    assert [kmer_tuple[0] for kmer_tuple in
            extract_kmers_from_string(sequence_with_wildcard, kmer_size,
                                      True)] == expected_kmers
    assert [kmer_tuple[0] for kmer_tuple in
            extract_kmers_from_string(sequence_with_wildcard, kmer_size,
                                      False)] == expected_kmers_with_wildcard

    sequence_below_kmer_size = "ACTGCCC"
    assert [kmer_tuple for kmer_tuple in
            extract_kmers_from_string(sequence_below_kmer_size, kmer_size,
                                      True)] == []

    sequence_with_one_kmer = "ACTGCCCG"
    expected_kmers = ["ACTGCCCG"]
    assert [kmer_tuple[0] for kmer_tuple in
            extract_kmers_from_string(sequence_with_one_kmer, kmer_size,
                                      True)] == expected_kmers


def test_build_reference_valid():
    """
    This function makes sure the reference building step is as expected
    :return:
    """
    reference = build_testing_reference()

    assert reference.kmer_db == FULL_FASTA_EXPECTED_KMER_DB
    assert reference.genomes_db["Moose"].unique_kmers == 9
    assert reference.genomes_db["Moose"].total_bases == 13
    assert reference.genomes_db["Turtle"].unique_kmers == 3
    assert reference.genomes_db["Turtle"].multi_mapping_kmers == 6
    assert reference.genomes_db["Mouse"].unique_kmers == 0 and \
           reference.genomes_db["Mouse"].kmers_set == {"AAAA"}
    assert reference.genomes_db["Virus"].unique_kmers == 0
    assert reference.genomes_db["Virus"].multi_mapping_kmers == 0
    assert reference.genomes_db["Virus"].total_bases == 3

    # Checking the final json output of the reference_file
    reference_output = reference.to_json()
    assert reference_output == json.dumps(FULL_FASTA_EXPECTED_OUTPUT)


def test_build_reference_invalid():
    """
    This function checks that the reference building step is stopped since
    there are two genomes of the same name
    :return:
    """
    # Duplicate_genomes
    kmer_size = 4
    genome_fasta_file = create_temporary_file(DUPLICATE_GENOMES_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    reference = KmerReference(kmer_size)
    build_successful = reference.build_kmer_reference(genome_fasta_file)
    reference.calculate_kmers_type()
    os.remove(genome_fasta_file)

    assert build_successful == False


def test_filter_similar():
    """
    This function checks the functionality of the filter similar extension
    With the given fasta file it shows which genomes were filtered out due to similarity
    :return:
    """
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    reference = KmerReference(KMER_SIZE)
    build_successful = reference.build_kmer_reference(genome_fasta_file)
    reference.calculate_kmers_type(count_duplicates=False)
    os.remove(genome_fasta_file)

    assert build_successful

    reference.filter_genomes_logic(0.95)
    assert reference.print_similarity_results() == json.dumps(SIMILARITY_OUTPUT)
    assert len(reference.genomes_db) == 4
