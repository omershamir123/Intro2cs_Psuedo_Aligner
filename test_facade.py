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
    reference = build_testing_reference()
    genome_fasta_file = create_temporary_file(FULL_FASTA_FILE,
                                              FASTA_FILE_TYPES[0])
    facade_reference = facade.build_reference(reference.kmer_size,
                                              genome_fasta_file)
    os.remove(genome_fasta_file)
    assert reference.kmer_db == facade_reference.kmer_db


def test_reference_command_3_1(monkeypatch):
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
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "dumpref"])
    args = main.readargs(sys.argv[1:])
    dumpref_command(args)
    assert capsys.readouterr().out.strip() == "No reference file given - invalid input for the dumpref command"


def test_dumpref_command_wrong_params(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "dumpref", "--genomefile",
                         "Gsdfgsd", "-k", "4", "--referencefile",
                         "fsdgfdsfgsd"])
    args = main.readargs(sys.argv[1:])
    dumpref_command(args)
    assert capsys.readouterr().out.strip() == "Both a genome file and a reference file were provided - invalid input for the dumpref command"


def test_align_insuffiecient_params(monkeypatch, capsys):
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "align"])
    args = main.readargs(sys.argv[1:])
    align_command(args)
    assert capsys.readouterr().out.strip() == "For the alignment command - an align file path and a reads file must be provided"


def test_align_3_4(monkeypatch):
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
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "dumpalign"])
    args = main.readargs(sys.argv[1:])
    dumpalign_command(args)
    assert capsys.readouterr().out.strip() == "No reads file given - invalid input for the dumpalign command"

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
    assert json.loads(str.split(capsys.readouterr().out, '{\n    "Statistics')[
                          0].strip()) == json.loads(
        json.dumps(COVERAGE_OUTPUT_FOR_READS))


def test_function_calls(monkeypatch):
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
    monkeypatch.setattr(sys, "argv",
                        ["main.py", "-t", "re1"])
    args = main.readargs(sys.argv[1:])
    start_program(args)
    assert capsys.readouterr().out == "Invalid task\n"