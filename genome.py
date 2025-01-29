# FILE : genome.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: this file handles the class of the individual genomes
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
from typing import Dict

from program_constants import ALLOWED_DNA_VALUES
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

    def genome_ref_to_dict(self) -> Dict[str, int]:
        return {"total_bases": self._total_bases,
                "unique_kmers": self._unique_kmers,
                "multi_mapping_kmers": self._multi_mapping_kmers}

    def genome_mapped_to_dict(self) -> Dict[str, int]:
        return {"unique_reads": self._unique_reads,"ambiguous_reads": self._ambiguous_reads}
