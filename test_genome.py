# FILE : test_genome.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs FInal Project 2025
# DESCRIPTION: This module handles all the tests for the ReferencedGenome class
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:

import pytest

from genome import ReferencedGenome
from program_constants import UNIQUE_KMER, MULTI_MAP_KMER


def test_init():
    genome1 = ReferencedGenome("genome1", "ATCNG", 2)
    assert genome1 is not None and genome1.sequence == "ATCNG"
    with pytest.raises(ValueError):
        ReferencedGenome("genome1", "ATCfNG", 3)

def test_kmer_add_to_genome():
    genome1 = ReferencedGenome("genome1", "ATCNG", 2)

    genome1.add_kmer_to_genome_mapping("KMER1", UNIQUE_KMER)
    assert genome1.unique_kmers == 1 and genome1.unique_kmers_set == {"KMER1"}

    genome1.add_kmer_to_genome_mapping("KMER2", UNIQUE_KMER)
    genome1.add_kmer_to_genome_mapping("KMER1", MULTI_MAP_KMER)
    assert genome1.unique_kmers == 1 and genome1.unique_kmers_set == {"KMER2"}
    assert genome1.multi_mapping_kmers == 1 and genome1.multi_mapping_kmers_set == {"KMER1"}

    genome1.add_kmer_to_genome_mapping("KMER2", MULTI_MAP_KMER)
    genome1.add_kmer_to_genome_mapping("KMER2", UNIQUE_KMER)
    assert genome1.unique_kmers == 1 and genome1.unique_kmers_set == {"KMER2"}
    assert genome1.multi_mapping_kmers == 1 and genome1.multi_mapping_kmers_set == {"KMER1"}


def test_genome_to_dict_format():
    genome1 = ReferencedGenome("genome1", "ATCNG", 2)
    genome1.add_kmer_to_genome_mapping("KMER1", UNIQUE_KMER)
    genome1.add_kmer_to_genome_mapping("KMER2", MULTI_MAP_KMER)
    genome1.unique_reads += 1
    genome1.ambiguous_reads += 10

    assert genome1.genome_ref_to_dict() == {"total_bases": 5,
                "unique_kmers": 1,
                "multi_mapping_kmers": 1}

    assert genome1.genome_mapped_to_dict() == {"unique_reads": 1,
                "ambiguous_reads": 10}




