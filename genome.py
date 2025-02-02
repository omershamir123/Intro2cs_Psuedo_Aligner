# FILE : genome.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: this file handles the class of the individual genomes
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
from typing import Dict, Set

import numpy as np
from numpy import ndarray

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
        self._unique_reads_reverse: int = 0
        self._ambiguous_reads_reverse: int = 0
        self._unique_kmers_set: Set[str] = set()
        self._multi_mapping_kmers_set: Set[str] = set()
        self._kmers_set: Set[str] = set()
        self._unique_coverage_positions: ndarray = np.zeros(1)
        self._ambiguous_coverage_positions: ndarray = np.zeros(1)

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
    @unique_kmers.setter
    def unique_kmers(self, value: int) -> None:
        self._unique_kmers = value

    @property
    def multi_mapping_kmers(self) -> int:
        return self._multi_mapping_kmers

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
    def unique_reads_reverse(self) -> int:
        return self._unique_reads_reverse

    @unique_reads_reverse.setter
    def unique_reads_reverse(self, value: int):
        self._unique_reads_reverse = value

    @property
    def ambiguous_reads_reverse(self) -> int:
        return self._ambiguous_reads_reverse

    @ambiguous_reads_reverse.setter
    def ambiguous_reads_reverse(self, value: int) -> None:
        self._ambiguous_reads_reverse = value

    @property
    def kmers_set(self) -> Set[str]:
        return self._kmers_set

    @property
    def unique_kmers_set(self) -> Set[str]:
        return self._unique_kmers_set

    @property
    def multi_mapping_kmers_set(self) -> Set[str]:
        return self._multi_mapping_kmers_set

    @property
    def unique_coverage_positions(self) -> ndarray:
        return self._unique_coverage_positions

    @unique_coverage_positions.setter
    def unique_coverage_positions(self, value: ndarray) -> None:
        self._unique_coverage_positions = value

    @property
    def ambiguous_coverage_positions(self) -> ndarray:
        return self._ambiguous_coverage_positions

    @ambiguous_coverage_positions.setter
    def ambiguous_coverage_positions(self, value: ndarray) -> None:
        self._ambiguous_coverage_positions = value

    def add_kmer_to_genome_mapping(self, kmer_value: str) -> None:
        """
        This function updates the genome's kmer mapping - adds kmer_value to the KMER_TYPE mapping
        Can be either UNIQUE_KMER or MULTI_MAP_KMER
        :param kmer_value: the kmer_value
        :return: None
        """
        self.kmers_set.add(kmer_value)

    def genome_ref_to_dict(self) -> Dict[str, int]:
        return {"total_bases": self._total_bases,
                "unique_kmers": self._unique_kmers,
                "multi_mapping_kmers": self._multi_mapping_kmers}

    def genome_mapped_to_dict(self, is_reversed: bool = False) -> Dict[
        str, int]:
        if is_reversed:
            return {"unique_reads": self._unique_reads_reverse,
                    "ambiguous_reads": self._ambiguous_reads_reverse}
        else:
            return {"unique_reads": self._unique_reads,
                    "ambiguous_reads": self._ambiguous_reads}

    def initialize_coverage_arrays(self) -> None:
        self._unique_coverage_positions = np.zeros(self._total_bases, dtype=int)
        self._ambiguous_coverage_positions = np.zeros(self._total_bases,
                                                      dtype=int)

    def genome_coverage_stats(self, min_coverage: int) -> dict:
        covered_bases_unique = int(
            (self._unique_coverage_positions >= min_coverage).sum())
        covered_bases_ambiguous = int(
            (self._ambiguous_coverage_positions >= min_coverage).sum())
        mean_coverage_unique = round(self._unique_coverage_positions.mean(),
                                     1)
        mean_coverage_ambiguous = round(
            self._ambiguous_coverage_positions.mean(), 1)
        return {"covered_bases_unique": covered_bases_unique,
                "covered_bases_ambiguous": covered_bases_ambiguous,
                "mean_coverage_unique": mean_coverage_unique,
                "mean_coverage_ambiguous": mean_coverage_ambiguous}


    def genome_full_coverage_stats(self) -> dict:
        return {"unique_cov": self._unique_coverage_positions.tolist(),
                "ambiguous_cov": self._ambiguous_coverage_positions.tolist()}
