# FILE : test_facade.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: This module checks the facade and all the possible inputs
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
import json
import os
import sys

import facade
import file_handlers
import main
from facade import reference_command, dumpref_command, align_command, \
    dumpalign_command, start_program
from file_handlers import decompress_pickle_file
from program_constants import FASTA_FILE_TYPES, KDB_FILE_TYPES, ALN_FILE_TYPES, \
    FASTQ_FILE_TYPES
from pseudo_aligner import align_algorithm
from test_file_handlers import create_temporary_file, EMPTY_FILE
from test_kmer_reference import build_testing_reference, FULL_FASTA_FILE
from test_pseudo_aligner import FASTQ_FILE_SINGLE_READ

DUMPALIGN_OUTPUT_OF_FULL_FASTA_FILE = {
    "Statistics": {"unique_mapped_reads": 0, "ambiguous_mapped_reads": 1,
                   "unmapped_reads": 0},
    "Summary": {"Mouse": {"unique_reads": 0, "ambiguous_reads": 0},
                "Duck": {"unique_reads": 0, "ambiguous_reads": 1},
                "Otter": {"unique_reads": 0, "ambiguous_reads": 0},
                "Turtle": {"unique_reads": 0, "ambiguous_reads": 1},
                "Moose": {"unique_reads": 0, "ambiguous_reads": 1},
                "Virus": {"unique_reads": 0, "ambiguous_reads": 0}}}

COVERAGE_EXAMPLE_READS = (b"@READ1\nATCGAAAATTTTACAC\n+\n1111111111111111\n"
                          b"@READ2\nGTGTAAAATTTTCGAT\n+\n1111111111111111")

COVERAGE_OUTPUT_FOR_READS = {"Coverage": {
    "Mouse": {"covered_bases_unique": 0, "covered_bases_ambiguous": 0,
              "mean_coverage_unique": 0.0, "mean_coverage_ambiguous": 0.0},
    "Duck": {"covered_bases_unique": 0, "covered_bases_ambiguous": 0,
             "mean_coverage_unique": 1.0, "mean_coverage_ambiguous": 0.0},
    "Otter": {"covered_bases_unique": 0, "covered_bases_ambiguous": 0,
              "mean_coverage_unique": 0.0, "mean_coverage_ambiguous": 0.0},
    "Turtle": {"covered_bases_unique": 12, "covered_bases_ambiguous": 0,
               "mean_coverage_unique": 2.0, "mean_coverage_ambiguous": 0.0},
    "Moose": {"covered_bases_unique": 0, "covered_bases_ambiguous": 0,
              "mean_coverage_unique": 0.9, "mean_coverage_ambiguous": 0.0},
    "Virus": {"covered_bases_unique": 0, "covered_bases_ambiguous": 0,
              "mean_coverage_unique": 0.0, "mean_coverage_ambiguous": 0.0}},
    "Details": {"Mouse": {
        "unique_cov": [0, 0, 0, 0, 0, 0, 0, 0],
        "ambiguous_cov": [0, 0, 0, 0, 0, 0, 0, 0]},
        "Duck": {"unique_cov": [1, 1, 1, 1],
                 "ambiguous_cov": [0, 0, 0,
                                   0]},
        "Otter": {
            "unique_cov": [0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0, 0,
                           0, 0, 0, 0, 0, 0,
                           0],
            "ambiguous_cov": [0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0, 0,
                              0, 0, 0]},
        "Turtle": {
            "unique_cov": [2, 2, 2, 2, 2, 2, 2,
                           2, 2, 2, 2, 2],
            "ambiguous_cov": [0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0,
                              0]}, "Moose": {
            "unique_cov": [0, 1, 1, 1, 1, 1, 1, 1, 1,
                           1, 1, 1, 1],
            "ambiguous_cov": [0, 0, 0, 0, 0, 0, 0, 0,
                              0, 0, 0, 0, 0]},
        "Virus": {"unique_cov": [0, 0, 0],
                  "ambiguous_cov": [0, 0, 0]}}}


def test_build_reference():
    """
    This basic test checks that the building reference function works correctly
    :return: None
    """
    reference = build_testing_reference()
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    facade_reference = facade.build_reference(reference.kmer_size,
                                              genome_fasta_file)
    os.remove(genome_fasta_file)
    assert reference.kmer_db == facade_reference.kmer_db


def test_reference_command_3_1(monkeypatch):
    """
    This test checks the 3.1 command - the reference command
    :param monkeypatch: mock used to convert argv
    :return: None
    """
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    kdb_file = create_temporary_file(EMPTY_FILE, KDB_FILE_TYPES[0])
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "reference", "--genomefile",
                         genome_fasta_file, "-k", "4", "--referencefile",
                         kdb_file])
    args = main.readargs(sys.argv[1:])
    reference = build_testing_reference()

    reference_command(args)
    reference_from_kdb_file = file_handlers.decompress_pickle_file(kdb_file,
                                                                   KDB_FILE_TYPES)
    os.remove(genome_fasta_file)
    os.remove(kdb_file)

    assert reference.kmer_db == reference_from_kdb_file.kmer_db


def test_reference_insufficient_arguments(monkeypatch):
    """
    This function checks the call for a reference command without a genome file
    :param monkeypatch: mock used to convert argv
    :return: None
    """
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    kdb_file = create_temporary_file(EMPTY_FILE, KDB_FILE_TYPES[0])
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "reference",
                         "-k", "4", "--referencefile",
                         kdb_file])
    args = main.readargs(sys.argv[1:])

    reference_command(args)

    with (open(kdb_file, "rb") as kdb_file_handler):
        kdb_file_content = kdb_file_handler.read()

    assert kdb_file_content == EMPTY_FILE
    os.remove(genome_fasta_file)
    os.remove(kdb_file)


def test_dumpref_command_3_2(monkeypatch, capsys):
    """
    This function tests the 3.2 command - the dumpref with a given kdb file
    :param monkeypatch: mock used to convert argv
    :param capsys: used to capture stdout
    :return: None
    """
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    kdb_file = create_temporary_file(EMPTY_FILE, KDB_FILE_TYPES[0])
    reference = build_testing_reference()
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "reference", "--genomefile",
                         genome_fasta_file, "-k", "4", "--referencefile",
                         kdb_file])
    args = main.readargs(sys.argv[1:])
    reference_command(args)

    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "dumpref", "--referencefile",
                         kdb_file])
    args = main.readargs(sys.argv[1:])
    dumpref_command(args)

    assert json.loads(capsys.readouterr().out) == json.loads(
        reference.to_json())

    os.remove(kdb_file)
    os.remove(genome_fasta_file)


def test_dumpref_command_3_3(monkeypatch, capsys):
    """
    This test checks the output of the dumpref command with  a genome file a kmersize
    :param monkeypatch: mock used to convert argv
    :param capsys: used to capture stdout
    :return: None
    """
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    reference = build_testing_reference()
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "dumpref", "--genomefile",
                         genome_fasta_file, "-k", "4"])
    args = main.readargs(sys.argv[1:])
    dumpref_command(args)
    assert json.loads(capsys.readouterr().out) == json.loads(
        reference.to_json())

    os.remove(genome_fasta_file)


def test_dumpref_command_no_ref(monkeypatch, capsys):
    """
    This test checks the output of the dumpref command without any parameter
    :param monkeypatch: mock used to convert argv
    :param capsys: used to capture stdout
    :return: None
    """
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "dumpref"])
    args = main.readargs(sys.argv[1:])
    dumpref_command(args)
    assert capsys.readouterr().out.strip() == "No reference file given - invalid input for the dumpref command"


def test_dumpref_command_wrong_params(monkeypatch, capsys):
    """
    This test checks the output of the dumpref command with the wrong parameters
    genome file kmer size and kdb file given - invalid input for the dumpref command
    :param monkeypatch: mock used to convert argv
    :param capsys: used to capture stdout
    :return: None
    """
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "dumpref", "--genomefile",
                         "Gsdfgsd", "-k", "4", "--referencefile",
                         "fsdgfdsfgsd"])
    args = main.readargs(sys.argv[1:])
    dumpref_command(args)
    assert capsys.readouterr().out.strip() == "Both a genome file and a reference file were provided - invalid input for the dumpref command"


def test_align_insuffiecient_params(monkeypatch, capsys):
    """
    This test checks the output of the align command with insufficient parameters
    No parameters given
    :param monkeypatch: mock used to convert argv
    :param capsys: used to capture stdout
    :return: None
    """
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "align"])
    args = main.readargs(sys.argv[1:])
    align_command(args)
    assert capsys.readouterr().out.strip() == "For the alignment command - an align file path and a reads file must be provided"


def test_align_3_4(monkeypatch):
    """
    This function tests the 3.4 command - the align command which creates an align file
    :param monkeypatch: mock used to convert argv
    :return: None
    """
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    align_file = create_temporary_file(EMPTY_FILE, ALN_FILE_TYPES[0])
    fastq_file = create_temporary_file(FASTQ_FILE_SINGLE_READ,
                                       FASTQ_FILE_TYPES[0])

    kdb_file = create_temporary_file(EMPTY_FILE, KDB_FILE_TYPES[0])
    reference = build_testing_reference()
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "reference", "--genomefile",
                         genome_fasta_file, "-k", "4", "--referencefile",
                         kdb_file])
    args = main.readargs(sys.argv[1:])
    reference_command(args)

    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "align", "--referencefile",
                         kdb_file, "-a", align_file,
                         "--reads",
                         fastq_file, "-m", "1", "-p", "1"])
    args = main.readargs(sys.argv[1:])
    align_command(args)
    expected_align_output = align_algorithm(fastq_file, reference, 1,
                                            1).convert_to_aln_object(False)
    align_command_output = decompress_pickle_file(align_file, ALN_FILE_TYPES)
    key_to_exclude = "genomes_db"  # the genomes all contain different objects hence a simple == operator won't work
    assert {k: v for k, v in expected_align_output.__dict__.items() if
            k != key_to_exclude} == {k: v for k, v in
                                     align_command_output.__dict__.items() if
                                     k != key_to_exclude}

    os.remove(align_file)
    os.remove(genome_fasta_file)
    os.remove(fastq_file)


def test_align_3_5(monkeypatch, capsys):
    """
    This function tests the 3.5 command - align command with a given genome file and kmer zise
    :param monkeypatch: mock used to convert argv
    :param capsys: used to capture stdout
    :return: None
    """
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    align_file = create_temporary_file(EMPTY_FILE, ALN_FILE_TYPES[0])
    fastq_file = create_temporary_file(FASTQ_FILE_SINGLE_READ,
                                       FASTQ_FILE_TYPES[0])
    reference = build_testing_reference()
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "align", "--genomefile",
                         genome_fasta_file, "-k", "4", "-a", align_file,
                         "--reads",
                         fastq_file, "-m", "1", "-p", "1"])
    args = main.readargs(sys.argv[1:])
    align_command(args)
    expected_align_output = align_algorithm(fastq_file, reference, 1,
                                            1).convert_to_aln_object(False)
    align_command_output = decompress_pickle_file(align_file, ALN_FILE_TYPES)
    key_to_exclude = "genomes_db"  # the genomes all contain different objects hence a simple == operator won't work
    assert {k: v for k, v in expected_align_output.__dict__.items() if
            k != key_to_exclude} == {k: v for k, v in
                                     align_command_output.__dict__.items() if
                                     k != key_to_exclude}

    os.remove(align_file)
    os.remove(genome_fasta_file)
    os.remove(fastq_file)


def test_dumpalign_invalid_input(monkeypatch, capsys):
    """
    This test checks the output of the dumpalign command with invalid parameters
    There are various types of invalid input - each is explained in the response the user gets
    :param monkeypatch: mock used to convert argv
    :param capsys: used to capture stdout
    :return: None
    """
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "dumpalign"])
    args = main.readargs(sys.argv[1:])
    dumpalign_command(args)
    assert capsys.readouterr().out.strip() == "No reads and no aln file given - invalid input for the dumpalign command"

    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "dumpalign", "-a", "njl", "-r",
                         "afsfrsd"])
    args = main.readargs(sys.argv[1:])
    dumpalign_command(args)
    assert capsys.readouterr().out.strip() == "For the dumpalign command - if -a is a given flag, no other flags can be inputted"

    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "dumpalign", "--genomefile",
                         genome_fasta_file, "-k",
                         "4", "--reads", "reads.fq", "-m", "-1"])
    args = main.readargs(sys.argv[1:])
    dumpalign_command(args)
    assert capsys.readouterr().out.strip() == "Invalid unique threshold value provided"


def test_dumpalign_3_6(monkeypatch, capsys):
    """
    This function tests the 3.6 command - dumpalign with a given aln file
    :param monkeypatch: mock used to convert argv
    :param capsys: used to capture stdout
    :return: None
    """
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    align_file = create_temporary_file(EMPTY_FILE, ALN_FILE_TYPES[0])
    fastq_file = create_temporary_file(FASTQ_FILE_SINGLE_READ,
                                       FASTQ_FILE_TYPES[0])

    kdb_file = create_temporary_file(EMPTY_FILE, KDB_FILE_TYPES[0])
    reference = build_testing_reference()
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "reference", "--genomefile",
                         genome_fasta_file, "-k", "4", "--referencefile",
                         kdb_file])
    args = main.readargs(sys.argv[1:])
    reference_command(args)

    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "align", "--referencefile",
                         kdb_file, "-a", align_file,
                         "--reads",
                         fastq_file, "-m", "1", "-p", "1"])
    args = main.readargs(sys.argv[1:])
    align_command(args)

    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "dumpalign", "-a", align_file, "-m",
                         "1", "-p", "1"])
    args = main.readargs(sys.argv[1:])

    dumpalign_command(args)
    print(capsys.readouterr().out)
    assert json.loads(capsys.readouterr().out) == json.loads(
        json.dumps(DUMPALIGN_OUTPUT_OF_FULL_FASTA_FILE))

    os.remove(align_file)
    os.remove(fastq_file)
    os.remove(kdb_file)
    os.remove(genome_fasta_file)


def test_dumpalign_3_7(monkeypatch, capsys):
    """
    This function tests the 3.7 command which receives a reference and a reads file
    :param monkeypatch: mock used to convert argv
    :param capsys: used to capture stdout
    :return: None
    """
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    fastq_file = create_temporary_file(FASTQ_FILE_SINGLE_READ,
                                       FASTQ_FILE_TYPES[0])

    kdb_file = create_temporary_file(EMPTY_FILE, KDB_FILE_TYPES[0])
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "reference", "--genomefile",
                         genome_fasta_file, "-k", "4", "--referencefile",
                         kdb_file])
    args = main.readargs(sys.argv[1:])
    reference_command(args)

    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "dumpalign", "--referencefile",
                         kdb_file, "--reads", fastq_file, "-m",
                         "1", "-p", "1"])
    args = main.readargs(sys.argv[1:])

    dumpalign_command(args)

    assert json.loads(capsys.readouterr().out) == json.loads(
        json.dumps(DUMPALIGN_OUTPUT_OF_FULL_FASTA_FILE))

    os.remove(fastq_file)
    os.remove(kdb_file)
    os.remove(genome_fasta_file)


def test_dumpalign_3_8(monkeypatch, capsys):
    """
    This function tests the 3.8 command (dumpalign) which receives a genome file and a reads file as well as a kmer size
    :param monkeypatch: mock used to convert argv
    :param capsys: used to capture stdout
    :return: None
    """
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    fastq_file = create_temporary_file(FASTQ_FILE_SINGLE_READ,
                                       FASTQ_FILE_TYPES[0])

    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "dumpalign", "--genomefile",
                         genome_fasta_file, "-k", "4", "--reads", fastq_file,
                         "-m",
                         "1", "-p", "1"])
    args = main.readargs(sys.argv[1:])

    dumpalign_command(args)

    assert json.loads(capsys.readouterr().out) == json.loads(
        json.dumps(DUMPALIGN_OUTPUT_OF_FULL_FASTA_FILE))

    os.remove(fastq_file)
    os.remove(genome_fasta_file)


def test_4_3_coverage_dumplalign(monkeypatch, capsys):
    """
    This function tests the 4.3 extension of the coverage.
    :param monkeypatch: mock used to convert argv
    :param capsys: used to capture stdout
    :return: None
    """
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    fastq_file = create_temporary_file(COVERAGE_EXAMPLE_READS,
                                       FASTQ_FILE_TYPES[0])

    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "dumpalign", "--genomefile",
                         genome_fasta_file, "-k", "4", "--reads", fastq_file,
                         "-m",
                         "1", "-p", "1", "--coverage", "--full-coverage",
                         "--min-coverage", "2"])
    args = main.readargs(sys.argv[1:])

    dumpalign_command(args)
    assert json.loads(str.split(capsys.readouterr().out, '\n')[
                          1].strip()) == json.loads(
        json.dumps(COVERAGE_OUTPUT_FOR_READS))


def test_4_3_coverage_dumpalign_specific_genome(monkeypatch, capsys):
    """
    This function tests the 4.3 extension of the coverage with a parameter of specific genomes to check their coverage.
    :param monkeypatch: mock used to convert argv
    :param capsys: used to capture stdout
    :return: None
    """
    genome_file = b">Cat\nAACCTTTTTNTNTGT\n>Mouse\nAAAGTTGGTTGCAAACCC\n"

    fastq_file = (
        b"@READ1\nAAAAAAAAAAA\n+\nABCDEFGHIJK\n@READ2\nGGCCCAAAG\n+\n123456789\n@Read3\n"
        b"AACCTTTTT\n+\n986533332")
    reads_path = create_temporary_file(fastq_file, FASTQ_FILE_TYPES[0])
    genome_path = create_temporary_file(genome_file, FASTA_FILE_TYPES[0])

    monkeypatch.setattr(sys, 'argv',
                        ["main.py", "-t", "dumpalign", '-g', genome_path, '-k',
                         "3",
                         '--coverage', '--genomes',
                         "Mouse", "--reads",
                         reads_path, "--full-coverage"])
    args = main.readargs()
    start_program(args)
    os.remove(genome_path)
    os.remove(reads_path)

    coverage_output = str.split(capsys.readouterr().out.strip(), '\n', 1)[1]

    expected_coverage_result = {
        "Coverage": {
            "Mouse": {
                "covered_bases_unique": 11,
                "covered_bases_ambiguous": 0,
                "mean_coverage_unique": 0.9,
                "mean_coverage_ambiguous": 0.0
            }
        },
        "Details": {
            "Mouse": {
                "unique_cov": [2, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 1, 2, 2, 2, 1, 1, 1],
                "ambiguous_cov": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            }
        }
    }

    assert json.loads(coverage_output) == json.loads(
        json.dumps(expected_coverage_result))


def test_function_calls(monkeypatch):
    """
    This function tests the function calls - calls the start_program method instead of the
    specific funtion of the given command
    :param monkeypatch: mock used to convert argv
    :return: None
    """
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    kdb_file = create_temporary_file(EMPTY_FILE, KDB_FILE_TYPES[0])
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "reference", "--genomefile",
                         genome_fasta_file, "-k", "4", "--referencefile",
                         kdb_file])
    args = main.readargs(sys.argv[1:])
    reference = build_testing_reference()

    start_program(args)
    reference_from_kdb_file = file_handlers.decompress_pickle_file(kdb_file,
                                                                   KDB_FILE_TYPES)
    os.remove(genome_fasta_file)
    os.remove(kdb_file)

    assert reference.kmer_db == reference_from_kdb_file.kmer_db


def test_function_calls_invalid_task(monkeypatch, capsys):
    """
    This funciton checks the case in which a given task isnn't a valid task
    :param monkeypatch: mock used to convert argv
    :param capsys: used to capture stdout
    :return: None
    """
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "re1"])
    args = main.readargs(sys.argv[1:])
    start_program(args)
    assert capsys.readouterr().out == "Invalid task\n"
