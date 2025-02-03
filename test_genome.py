# FILE : test_genome.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs FInal Project 2025
# DESCRIPTION: This module handles all the tests for the ReferencedGenome class
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:

import pytest

from genome import ReferencedGenome


def test_init():
    """
    This function tests the initialization of the reference genome
    :return:
    """
    genome1 = ReferencedGenome("genome1", "ATCNG", 2)
    assert genome1 is not None and genome1.sequence == "ATCNG"
    with pytest.raises(ValueError):
        ReferencedGenome("genome1", "ATCfNG", 3)

def test_kmer_add_to_genome():
    """
    This function tests that the kmer addition to a genome works fine,
    and makes sure it's a set so no duplicates added
    :return:
    """
    genome1 = ReferencedGenome("genome1", "ATCNG", 2)

    genome1.add_kmer_to_genome_mapping("KMER1")
    assert genome1.kmers_set == {"KMER1"}

    genome1.add_kmer_to_genome_mapping("KMER2")
    genome1.add_kmer_to_genome_mapping("KMER1")
    assert "KMER1" in genome1.kmers_set and "KMER2" in genome1.kmers_set and len(genome1.kmers_set) == 2


def test_genome_to_dict_format():
    """
    This function tests that the genome_to_dict function works fine, with all the dictionaries
    coming out in the right format
    :return:
    """
    genome1 = ReferencedGenome("genome1", "ATCNG", 2)
    genome1.unique_reads += 1
    genome1.ambiguous_reads += 10
    genome1.unique_kmers += 1
    genome1.multi_mapping_kmers += 1
    assert genome1.genome_ref_to_dict() == {"total_bases": 5,
                "unique_kmers": 1,
                "multi_mapping_kmers": 1}

    assert genome1.genome_mapped_to_dict() == {"unique_reads": 1,
                "ambiguous_reads": 10}




