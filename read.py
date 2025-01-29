# FILE : read.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: this module contains the Read class which represents a fastq read and a class that maps a read to a genomes
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
from typing import List, Dict, Set

import numpy as np
from program_constants import *
from validators import validate_values_in_given_list, validate_above_value


class ReadKmerMapping:
    def __init__(self):
        self._specific_kmers = []
        self._unspecific_kmers = []
        self._specific_kmers_in_genomes: Dict[str, List[str]] = {}
        self._unspecific_kmers_in_genomes: Dict[str, List[str]] = {}

    @property
    def specific_kmers(self):
        return self._specific_kmers

    def add_specific_kmer(self, kmer: str) -> None:
        self._specific_kmers.append(kmer)

    @property
    def unspecific_kmers(self):
        return self._unspecific_kmers

    def add_unspecific_kmer(self, kmer: str) -> None:
        self._unspecific_kmers.append(kmer)

    @property
    def specific_kmers_in_genomes(self):
        return self._specific_kmers_in_genomes

    @property
    def unspecific_kmers_in_genomes(self):
        return self._unspecific_kmers_in_genomes


def reverse_complement(value: str) -> str:
    reversed_complement_list = [COMPLEMENTS_DICT[nucleotide] for nucleotide in
                                reversed(value)]
    return "".join(reversed_complement_list)


class Read:
    def __init__(self, identifier: str, value: str, quality: str) -> None:
        """
        This function initializes a Read instance.
        Three validity checks that are run are:
        1. The values of the read are all valid DNA nucleotides
        2. The quality line contains only valid characters - non have an ASCII code below 33
        3. The quality and the read value are of the same length
        :param identifier: the read identifier
        :param value: the read value
        :param quality: the read quality
        """
        # Validate that the read value consists of valid DNA nucleotides
        validate_values_in_given_list(value, ALLOWED_DNA_VALUES)

        self._identifier = identifier
        self._value = value
        self._quality: List[int] = [ord(letter) - 33 for letter in quality]

        # Check that the quality line doesn't contain any "negative" numbers
        map(lambda x: validate_above_value(x, 0, True), self._quality)
        if len(quality) != len(value):
            raise ValueError(
                "Read {} does not have matching read and quality lines".format(
                    identifier))

        self._length = len(value)
        self._read_status: READ_STATUS = UNMAPPED_READ
        self._mapped_genomes:Set[str] = set()
        self._specific_kmers = []
        self._unspecific_kmers = []
        self._reverse_complement: str = reverse_complement(value)

    @property
    def identifier(self) -> str:
        return self._identifier

    @property
    def value(self) -> str:
        return self._value

    @property
    def quality(self) -> List[int]:
        return self._quality

    @property
    def read_status(self) -> READ_STATUS:
        return self._read_status

    @read_status.setter
    def read_status(self, read_status: READ_STATUS) -> None:
        self._read_status = read_status

    @property
    def reverse_complement(self) -> str:
        return self._reverse_complement

    @property
    def mapped_genomes(self) -> Set[str]:
        return self._mapped_genomes

    def add_mapped_genome(self, genome: str) -> None:
        self._mapped_genomes.add(genome)


    def calculate_mean_quality(self, starting_index:int = 0, ending_index:int = -1) -> np.floating:
        """
        This function calculates the mean quality value for the read.
        Note: Assumes starting,ending >=0 + starting<=ending
        :param starting_index: Optional, 0 default - if supplied - calculates the mean quality from the starting index
        :param ending_index: Optional - if supplied - calculates the mean quality until the ending index (excluding)
        :return: the mean value
        """
        if ending_index == -1:
            ending_index = self._length

        return np.mean(np.array(self.quality[starting_index:ending_index]))
