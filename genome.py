# FILE : genome.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: this file handles the class of the individual genomes
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
from typing import Dict, Set

from program_constants import ALLOWED_DNA_VALUES, KMER_TYPE, UNIQUE_KMER, \
    MULTI_MAP_KMER
from validators import validate_values_in_given_list


class ReferencedGenome:
    """
    This class represents a raw genome - right after it's read from the fasta file
    """

    def __init__(self, identifier: str, sequence: str,
                 index_in_fasta: int) -> None:
        validate_values_in_given_list(sequence, ALLOWED_DNA_VALUES)
        self._identifier = identifier
        self._sequence = sequence
        self._index_in_fasta = index_in_fasta
        self._total_bases = len(sequence)
        self._unique_kmers: int = 0
        self._multi_mapping_kmers: int = 0
        self._unique_reads: int = 0
        self._ambiguous_reads: int = 0
        self._unique_kmers_set: Set[str] = set()
        self._multi_mapping_kmers_set: Set[str] = set()

    @property
    def identifier(self) -> str:
        return self._identifier

    @property
    def sequence(self) -> str:
        return self._sequence

    @property
    def index_in_fasta(self) -> int:
        return self._index_in_fasta

    @property
    def total_bases(self) -> int:
        return self._total_bases

    @property
    def unique_kmers(self) -> int:
        return self._unique_kmers

    @property
    def multi_mapping_kmers(self) -> int:
        return self._multi_mapping_kmers

    @unique_kmers.setter
    def unique_kmers(self, value: int) -> None:
        self._unique_kmers = value

    @multi_mapping_kmers.setter
    def multi_mapping_kmers(self, value: int) -> None:
        self._multi_mapping_kmers = value

    @property
    def unique_reads(self) -> int:
        return self._unique_reads

    @unique_reads.setter
    def unique_reads(self, value: int):
        self._unique_reads = value

    @property
    def ambiguous_reads(self) -> int:
        return self._ambiguous_reads

    @ambiguous_reads.setter
    def ambiguous_reads(self, value: int) -> None:
        self._ambiguous_reads = value

    @property
    def unique_kmers_set(self) -> Set[str]:
        return self._unique_kmers_set

    @property
    def multi_mapping_kmers_set(self) -> Set[str]:
        return self._multi_mapping_kmers_set

    def add_kmer_to_genome_mapping(self, kmer_value: str, kmer_type: KMER_TYPE) -> None:
        """
        This function updates the genome's kmer mapping - adds kmer_value to the KMER_TYPE mapping
        Can be either UNIQUE_KMER or MULTI_MAP_KMER
        :param kmer_value: the kmer_value
        :param kmer_type: the kmer_type
        :return: None
        """
        if kmer_type == UNIQUE_KMER:
            if kmer_value in self.multi_mapping_kmers_set:
                self.multi_mapping_kmers_set.remove(kmer_value)
                self.multi_mapping_kmers -= 1
            if kmer_value not in self.unique_kmers_set:
                self.unique_kmers_set.add(kmer_value)
                self.unique_kmers += 1

        if kmer_type == MULTI_MAP_KMER:
            if kmer_value in self.unique_kmers_set:
                self.unique_kmers_set.remove(kmer_value)
                self.unique_kmers -= 1
            if kmer_value not in self.multi_mapping_kmers_set:
                self.multi_mapping_kmers_set.add(kmer_value)
                self.multi_mapping_kmers += 1

    def genome_ref_to_dict(self) -> Dict[str, int]:
        return {"total_bases": self._total_bases,
                "unique_kmers": self._unique_kmers,
                "multi_mapping_kmers": self._multi_mapping_kmers}

    def genome_mapped_to_dict(self) -> Dict[str, int]:
        return {"unique_reads": self._unique_reads,
                "ambiguous_reads": self._ambiguous_reads}


