# FILE : genome.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: this file handles the class of the individual genomes
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
import json

class ReferencedGenome:
    """
    This class represents a raw genome - right after it's read from the fasta file
    """
    def __init__(self, identifier:str, sequence:str, index_in_fasta:int) -> None:
        self._identifier = identifier
        self._sequence = sequence
        self._index_in_fasta = index_in_fasta
        self._total_bases = len(sequence)
        self._unique_kmers = 0
        self._multi_mapping_kmers = 0
        self._unique_reads = 0
        self._ambiguous_reads = 0

    @property
    def identifier(self):
        return self._identifier

    @property
    def sequence(self):
        return self._sequence

    @property
    def index_in_fasta(self):
        return self._index_in_fasta

    @property
    def total_bases(self):
        return self._total_bases

    @property
    def unique_kmers(self):
        return self._unique_kmers

    @property
    def multi_mapping_kmers(self):
        return self._multi_mapping_kmers

    @unique_kmers.setter
    def unique_kmers(self, value):
        self._unique_kmers = value

    @multi_mapping_kmers.setter
    def multi_mapping_kmers(self, value):
        self._multi_mapping_kmers = value

    def genome_ref_to_dict(self):
        return {"total_bases": self._total_bases,
                           "unique_kmers": self._unique_kmers,
                           "multi_mapping_kmers": self._multi_mapping_kmers}



