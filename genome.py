# FILE : genome.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: this file handles the class of the individual genomes
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
ALLOWED_DNA_VALUES = ["A", "T", "C", "G", "N"]
WILDCARD_READINGS = ["N"]


class RawGenome:
    """
    This class represents a raw genome - right after it's read from the fasta file
    """
    def __init__(self, identifier:str, sequence:str):
        self._identifier = identifier
        self._sequence = sequence

    @property
    def identifier(self):
        return self._identifier

    @property
    def sequence(self):
        return self._sequence



class ReferencedGenomeStats:

    def __init__(self, value: str):
        self._total_bases = len(value)
        self._unique_kmers = 0
        self._multi_mapping_kmers = 0


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
