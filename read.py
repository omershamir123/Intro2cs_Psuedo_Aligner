# FILE : read.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: this module contains the Read class which represents a fastq read
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
from typing import List

from program_constants import *


def reverse_complement(value: str) -> str:
    reversed_complement_list = [COMPLEMENTS_DICT[nucleotide] for nucleotide in
                                reversed(value)]
    return "".join(reversed_complement_list)


class Read:
    def __init__(self, identifier: str, value: str, quality: str) -> None:
        self._identifier = identifier
        self._value = value
        self._quality: List[int] = [ord(letter) - 33 for letter in quality]
        self._length = len(value)
        self._read_status: READ_STATUS = UNMAPPED_READ
        self._reverse_complement: str = reverse_complement(value)

    @property
    def identifier(self) -> str:
        return self._identifier

    @property
    def value(self) -> str:
        return self._value

    @property
    def quality(self) -> str:
        return self._quality

    @property
    def read_status(self) -> READ_STATUS:
        return self._read_status

    @property
    def reverse_complement(self) -> str:
        return self._reverse_complement
