# FILE : test_pseudo_aligner.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: 
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
import json
import os

from kmer_reference import KmerReference
from program_constants import UNMAPPED_READ, UNIQUE_READ, AMBIGUOUS_READ, \
    FASTQ_FILE_TYPES
from pseudo_aligner import PseudoAlignerOutput, should_filter_read, map_read, \
    extract_and_map_kmers_from_read, _map_read_using_specific_kmers, \
    validate_uniqueness_using_unspecific, align_algorithm, extract_genomes_list, \
    initialize_filtering_stats
from read import Read, ReadKmerMapping
from test_file_handlers import create_temporary_file
from test_kmer_reference import build_testing_reference

FASTQ_FILE_SINGLE_READ = (b"@READ1\n"
                          b"ATCGAAAATTTTACAC\n"
                          b"+\n"
                          b"1111111111111111")

FASTQ_FILE_WITH_REVERSE_READS = (b"@READ1\n"
                                  b"GTGTAAAATTTTCGAT\n"
                                  b"+\n"
                                  b"1111111111111111\n"
                                  b"@READ2\n"
                                  b"GTGTGAAAATTTTCGAT\n"
                                  b"+\n"
                                  b"11111111111111111"
                                 )


def test_add_read_update_stats():
    """
    This function tests the update stats + add read functionality of the alignerOutput
    :return: None
    """
    read1 = Read("READ1", "ATCG", "1111")
    read2 = Read("READ2", "ATCGG", "11111")
    read3 = Read("READ3", "AATCG", "11111")
    read4 = Read("READ4", "AATCGG", "111111")
    kmer_size = 4
    reference = KmerReference(4)
    aligner_output = PseudoAlignerOutput(reference)

    aligner_output.add_read(read1)
    aligner_output.add_read(read2)
    aligner_output.add_read(read3)
    aligner_output.add_read(read4)

    read1.read_status = UNMAPPED_READ
    read2.read_status = UNIQUE_READ
    read3.read_status = UNIQUE_READ
    read4.read_status = AMBIGUOUS_READ

    aligner_output.update_reads_stats(read1.identifier)
    aligner_output.update_reads_stats(read2.identifier)
    aligner_output.update_reads_stats(read3.identifier)
    aligner_output.update_reads_stats(read4.identifier)

    assert len(aligner_output.reads) == 4
    assert aligner_output.unmapped_reads == 1
    assert aligner_output.unique_mapped_reads == 2
    assert aligner_output.ambiguous_mapped_reads == 1


def test_should_filter_read():
    """
    This function checks the functionality of should_filter_read
    :return:
    """
    read1 = Read("READ1", "ATCG", "1111")
    read2 = Read("READ2", "ATCGG", "00000")
    min_read_quality = 16

    assert should_filter_read(read1, min_read_quality) == False
    assert should_filter_read(read2, min_read_quality) == True


def test_map_read():
    """
    This function checks the functionality of map_read - maps the reads and updates the unique_reads + ambiguous reads
    :return:
    """
    read1 = Read("READ1", "ATCG", "1111")
    read2 = Read("READ2", "CCCC", "1111")
    read3 = Read("READ3", "TTTT", "1111")

    reference = build_testing_reference()
    map_read(read1, UNIQUE_READ, reference, "Duck")
    map_read(read2, UNMAPPED_READ, reference)

    map_read(read3, AMBIGUOUS_READ, reference, "Otter")
    map_read(read3, AMBIGUOUS_READ, reference, "Turtle")

    assert read1.read_status == UNIQUE_READ
    assert read2.read_status == UNMAPPED_READ
    assert read3.read_status == AMBIGUOUS_READ
    assert reference.genomes_db["Duck"].unique_reads == 1
    assert reference.genomes_db["Otter"].ambiguous_reads == 1
    assert reference.genomes_db["Turtle"].ambiguous_reads == 1


def test_extract_and_map_kmers_from_read():
    """
    This function checks the functionality of extract_and_map_kmers_from_read
    as well as them being either specific or unspecific
    :return:
    """
    read1 = Read("READ1", "ATCGAAAATTTTACAC", "ATCGAAAATTTTACAC")
    reference = build_testing_reference()
    aligner_output = PseudoAlignerOutput(reference)
    expected_kmer_mapping = ReadKmerMapping()

    expected_kmer_mapping.specific_kmers.append("ATCG")
    expected_kmer_mapping.specific_kmers.append("TTTA")
    expected_kmer_mapping.specific_kmers.append("ACAC")

    expected_kmer_mapping.specific_kmers_in_genomes["Duck"] = ["ATCG"]
    expected_kmer_mapping.specific_kmers_in_genomes["Moose"] = ["ACAC"]
    expected_kmer_mapping.specific_kmers_in_genomes["Turtle"] = ["TTTA"]

    expected_kmer_mapping.unspecific_kmers.extend(["AAAA", "TTTT"])
    expected_kmer_mapping.unspecific_kmers_in_genomes["Turtle"] = ["AAAA",
                                                                   "TTTT"]
    expected_kmer_mapping.unspecific_kmers_in_genomes["Mouse"] = ["AAAA"]
    expected_kmer_mapping.unspecific_kmers_in_genomes["Otter"] = ["TTTT"]

    read_kmer_mapping = extract_and_map_kmers_from_read(read1, reference,
                                                        aligner_output)

    assert read_kmer_mapping.__dict__ == expected_kmer_mapping.__dict__


def test_extract_and_map_kmers_from_read_with_filter_kmer_quality():
    """
    This function checks the functionality of extract_and_map_kmers_from_read
    but with a filter quality threshold of an entire kmer
    :return:
    """
    read1 = Read("READ1", "ATCGAAAATTTTACAC", "0111111111111111")
    reference = build_testing_reference()
    aligner_output = PseudoAlignerOutput(reference)
    expected_kmer_mapping = ReadKmerMapping()

    expected_kmer_mapping.specific_kmers.append("TTTA")
    expected_kmer_mapping.specific_kmers.append("ACAC")

    expected_kmer_mapping.specific_kmers_in_genomes["Moose"] = ["ACAC"]
    expected_kmer_mapping.specific_kmers_in_genomes["Turtle"] = ["TTTA"]

    expected_kmer_mapping.unspecific_kmers.extend(["AAAA", "TTTT"])
    expected_kmer_mapping.unspecific_kmers_in_genomes["Turtle"] = ["AAAA",
                                                                   "TTTT"]
    expected_kmer_mapping.unspecific_kmers_in_genomes["Mouse"] = ["AAAA"]
    expected_kmer_mapping.unspecific_kmers_in_genomes["Otter"] = ["TTTT"]

    initialize_filtering_stats(aligner_output, min_kmer_quality=16)
    read_kmer_mapping = extract_and_map_kmers_from_read(read1, reference,
                                                        aligner_output,
                                                        min_kmer_quality=16)
    assert read_kmer_mapping.__dict__ == expected_kmer_mapping.__dict__
    # One kmer was filtered - ATCG due to poor quality
    assert aligner_output.filtered_quality_kmers == 1


def test_extract_and_map_kmers_from_read_with_full_filter():
    """
    This function tests that all the read's kmers are filtered out
    :return:
    """
    read1 = Read("READ1", "ATCGAAAATTTTACACNNNNNN", "0111111111111111AAAAAA")
    reference = build_testing_reference()
    aligner_output = PseudoAlignerOutput(reference)
    expected_kmer_mapping = ReadKmerMapping()

    initialize_filtering_stats(aligner_output, min_kmer_quality=17)

    # no kmer in the read will be with sufficient quality, though the read itself might be
    # due to the presence of N in the sequence
    read_kmer_mapping = extract_and_map_kmers_from_read(read1, reference,
                                                        aligner_output,
                                                        min_kmer_quality=17)

    assert read_kmer_mapping.__dict__ == expected_kmer_mapping.__dict__
    # All kmers extracted from the read were filtered
    assert aligner_output.filtered_quality_kmers == 13

    reference = build_testing_reference()
    aligner_output = PseudoAlignerOutput(reference)

    expected_kmer_mapping = ReadKmerMapping()

    expected_kmer_mapping.specific_kmers.append("TTTA")
    expected_kmer_mapping.specific_kmers.append("ACAC")

    expected_kmer_mapping.specific_kmers_in_genomes["Moose"] = ["ACAC"]
    expected_kmer_mapping.specific_kmers_in_genomes["Turtle"] = ["TTTA"]

    initialize_filtering_stats(aligner_output, min_kmer_quality=16, max_genomes=1)
    # The first kmer will be filtered out due to quality, two kmers will be filtered out due to max genomes
    read_kmer_mapping = extract_and_map_kmers_from_read(read1, reference,
                                                        aligner_output,
                                                        min_kmer_quality=16,
                                                        max_genomes=1)
    assert read_kmer_mapping.__dict__ == expected_kmer_mapping.__dict__
    assert aligner_output.filtered_quality_kmers == 1
    assert aligner_output.filtered_hr_kmers == 2


def test_map_read_using_specific_kmers_one_specific_unique():
    """
    This function checks the mapping of a specific kmer to a unique reads
    :return:
    """
    read1 = Read("READ1", "ACACACAC", "11111111")
    reference = build_testing_reference()
    aligner_output = PseudoAlignerOutput(reference)
    read_kmer_mapping = extract_and_map_kmers_from_read(read1, reference,
                                                        aligner_output)

    # This read should be mapped uniquely mapped to Moose
    _map_read_using_specific_kmers(read_kmer_mapping, read1, reference, 1)

    assert read1.read_status == UNIQUE_READ
    assert reference.genomes_db["Moose"].unique_reads == 1


def test_map_read_using_specific_kmers_unique():
    """
    This function tests the case of a read being mapped as unique using the m parameter
    :return:
    """
    read1 = Read("READ1", "ACACAATCG", "111111111")
    reference = build_testing_reference()
    aligner_output = PseudoAlignerOutput(reference)
    read_kmer_mapping = extract_and_map_kmers_from_read(read1, reference,
                                                        aligner_output)

    assert len(read_kmer_mapping.specific_kmers) == 3
    # This read should be mapped uniquely mapped to Moose
    _map_read_using_specific_kmers(read_kmer_mapping, read1, reference, 1)

    assert read1.read_status == UNIQUE_READ
    assert reference.genomes_db["Moose"].unique_reads == 1


def test_map_read_using_specific_kmers_ambiguous():
    """
    This function tests a case where a read is ambiguously mapped before checking
    unspecific kmers
    :return:
    """
    read1 = Read("READ1", "ATCGAAAATTTTACAC", "ATCGAAAATTTTACAC")
    reference = build_testing_reference()
    aligner_output = PseudoAlignerOutput(reference)
    read_kmer_mapping = extract_and_map_kmers_from_read(read1, reference,
                                                        aligner_output)

    # This read should be mapped ambiguously to Duck, Moose, Turtle
    _map_read_using_specific_kmers(read_kmer_mapping, read1, reference, 1)

    assert read1.read_status == AMBIGUOUS_READ
    assert reference.genomes_db["Duck"].ambiguous_reads == 1
    assert reference.genomes_db["Moose"].ambiguous_reads == 1
    assert reference.genomes_db["Turtle"].ambiguous_reads == 1


def test_validate_uniqueness_using_unspecific_kmers_changes_to_ambiguous():
    """
    This function tests a case where a read is uniquely mapped and doesn't stay that way after the validation
    using unspecific kmers
    :return:
    """
    read1 = Read("READ1", "ATCGAAAATTTTACACA", "ATCGAAAATTTTACACA")
    reference = build_testing_reference()
    aligner_output = PseudoAlignerOutput(reference)

    read_kmer_mapping = extract_and_map_kmers_from_read(read1, reference,
                                                        aligner_output)
    # In this case - the read will be mapped uniquely to Moose due to specific kmers
    # However, once running the validation of unspecific kmers - the read will be mapped ambiguosly
    # to both moose and turtle

    _map_read_using_specific_kmers(read_kmer_mapping, read1, reference, 1)
    assert read1.read_status == UNIQUE_READ
    validate_uniqueness_using_unspecific(read_kmer_mapping, read1, reference, 0)

    assert read1.read_status == AMBIGUOUS_READ
    assert reference.genomes_db["Moose"].ambiguous_reads == 1
    assert reference.genomes_db["Turtle"].ambiguous_reads == 1


def test_validate_uniqueness_using_unspecific_kmers_stays_unique():
    """
    This function tests a case where a read is uniquely mapped and stays that way after the validation
    using unspecific kmers
    :return:
    """
    read1 = Read("READ1", "ATCGAAAATTTTACACACA", "ATCGAAAATTTTACACACA")
    reference = build_testing_reference()
    aligner_output = PseudoAlignerOutput(reference)

    read_kmer_mapping = extract_and_map_kmers_from_read(read1, reference,
                                                        aligner_output)

    # In this case - the read will be mapped uniquely to Moose due to specific kmers
    # And it stays that way after the validation

    _map_read_using_specific_kmers(read_kmer_mapping, read1, reference, 1)
    assert read1.read_status == UNIQUE_READ
    validate_uniqueness_using_unspecific(read_kmer_mapping, read1, reference, 0)

    assert read1.read_status == UNIQUE_READ
    assert reference.genomes_db["Moose"].unique_reads == 1


def test_align_algorithm():
    """
    This funciton checks the functionality of a simple align algorithm (no extensions)
    :return:
    """
    fastq_file_path = create_temporary_file(FASTQ_FILE_SINGLE_READ,
                                            FASTQ_FILE_TYPES[0])
    reference = build_testing_reference()
    aligner_output = align_algorithm(fastq_file_path, reference, 1, 1)
    os.remove(fastq_file_path)

    expected_reads_stats = {"unique_mapped_reads": 0,
                            "ambiguous_mapped_reads": 1,
                            "unmapped_reads": 0}
    expected_mapped_genome_stats = {
        "Mouse": {"unique_reads": 0, "ambiguous_reads": 0},
        "Duck": {"unique_reads": 0, "ambiguous_reads": 1},
        "Otter": {"unique_reads": 0, "ambiguous_reads": 0},
        "Turtle": {"unique_reads": 0, "ambiguous_reads": 1},
        "Moose": {"unique_reads": 0, "ambiguous_reads": 1},
        "Virus": {"unique_reads": 0, "ambiguous_reads": 0}}

    expected_align_output = {"Statistics": expected_reads_stats,
                             "Summary": expected_mapped_genome_stats}

    aln_object = aligner_output.convert_to_aln_object(False)
    assert json.dumps(expected_align_output) == aln_object.to_json()


def test_align_algorithm_reverse():
    """
    This function checks the functionality of the align algorithm with the reverse complement extension
    :return:
    """
    fastq_file_path = create_temporary_file(FASTQ_FILE_WITH_REVERSE_READS,
                                            FASTQ_FILE_TYPES[0])
    reference = build_testing_reference()
    aligner_output = align_algorithm(fastq_file_path, reference, 1, 1,
                                     reverse_complement=True)
    os.remove(fastq_file_path)

    expected_reads_stats = {"unique_mapped_reads": 2,
                            "ambiguous_mapped_reads": 0,
                            "unmapped_reads": 0}
    expected_mapped_genome_stats = {
        "Mouse_F": {"unique_reads": 0, "ambiguous_reads": 0},
        "Duck_F": {"unique_reads": 0, "ambiguous_reads": 0},
        "Otter_F": {"unique_reads": 0, "ambiguous_reads": 0},
        "Turtle_F": {"unique_reads": 1, "ambiguous_reads": 0},
        "Moose_F": {"unique_reads": 0, "ambiguous_reads": 0},
        "Virus_F": {"unique_reads": 0, "ambiguous_reads": 0},
        "Mouse_R": {"unique_reads": 0, "ambiguous_reads": 0},
        "Duck_R": {"unique_reads": 0, "ambiguous_reads": 0},
        "Otter_R": {"unique_reads": 0, "ambiguous_reads": 0},
        "Turtle_R": {"unique_reads": 0, "ambiguous_reads": 0},
        "Moose_R": {"unique_reads": 1, "ambiguous_reads": 0},
        "Virus_R": {"unique_reads": 0, "ambiguous_reads": 0}
    }

    expected_align_output = {"Statistics": expected_reads_stats,
                             "Summary": expected_mapped_genome_stats}

    aln_object = aligner_output.convert_to_aln_object(True)
    assert json.dumps(expected_align_output) == aln_object.to_json()


def test_extract_genomes_list():
    """
    This function checks the functionality of a extraction of the genomes list for
    the coverage extension
    :return:
    """
    reference = build_testing_reference()
    expected_genome_list = ["Mouse", "Duck", "Otter", "Turtle", "Moose",
                            "Virus"]
    genomes_list = extract_genomes_list(None, reference)
    assert genomes_list == expected_genome_list

    expected_genome_list = ["Mouse", "Otter"]
    genomes_list = extract_genomes_list("Mouse,Otter", reference)
    assert genomes_list == expected_genome_list


def test_align_algorithm_read_filtered_due_to_quality():
    """
    This function checks the align algorithm with the case that a read is completely filtered out
    :return:
    """
    fastq_file_path = create_temporary_file(FASTQ_FILE_SINGLE_READ,
                                            FASTQ_FILE_TYPES[0])
    reference = build_testing_reference()
    aligner_output = align_algorithm(fastq_file_path, reference, 1, 1,
                                     min_read_quality=17)
    os.remove(fastq_file_path)

    expected_reads_stats = {"unique_mapped_reads": 0,
                            "ambiguous_mapped_reads": 0,
                            "unmapped_reads": 0,
                            "filtered_quality_reads": 1}
    expected_mapped_genome_stats = {
        "Mouse_F": {"unique_reads": 0, "ambiguous_reads": 0},
        "Duck_F": {"unique_reads": 0, "ambiguous_reads": 0},
        "Otter_F": {"unique_reads": 0, "ambiguous_reads": 0},
        "Turtle_F": {"unique_reads": 0, "ambiguous_reads": 0},
        "Moose_F": {"unique_reads": 0, "ambiguous_reads": 0},
        "Virus_F": {"unique_reads": 0, "ambiguous_reads": 0},
        "Mouse_R": {"unique_reads": 0, "ambiguous_reads": 0},
        "Duck_R": {"unique_reads": 0, "ambiguous_reads": 0},
        "Otter_R": {"unique_reads": 0, "ambiguous_reads": 0},
        "Turtle_R": {"unique_reads": 0, "ambiguous_reads": 0},
        "Moose_R": {"unique_reads": 0, "ambiguous_reads": 0},
        "Virus_R": {"unique_reads": 0, "ambiguous_reads": 0}
    }

    expected_align_output = {"Statistics": expected_reads_stats,
                             "Summary": expected_mapped_genome_stats}

    aln_object = aligner_output.convert_to_aln_object(True)
    assert json.dumps(expected_align_output) == aln_object.to_json()

