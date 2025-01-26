# FILE : kmer_reference.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: A file with the KmerReference Class
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
from typing import Dict, List, Generator, Tuple

import file_parsers
import validators
from genome import ReferencedGenomeStats, ALLOWED_DNA_VALUES, RawGenome, \
    WILDCARD_READINGS


def extract_kmers_from_string(sequence: str, kmer_size: int,
                              remove_wildcard: bool) -> \
        Generator[Tuple[str, int], None, None]:
    """
    This function receives a sequence and returns all kmers within that sequence.
    If a kmer contains a WILDCARD_READING - then the kmer will be thrown out
    :param sequence: the whole DNA sequence
    :param kmer_size: the length of the kmer to be extracted
    :param remove_wildcard: should we remove the kmers containing a wildcard
    :return: Tuple with the first parameter being the kmer, the second being its position in the initial sequence
    """
    if kmer_size > len(sequence):
        raise StopIteration
    else:
        current_kmer = sequence[:kmer_size]
        if not (remove_wildcard and any(
                nucleotide in WILDCARD_READINGS for nucleotide in
                current_kmer)):
            yield current_kmer, 0
    for i in range(1, len(sequence) - kmer_size + 1):
        current_kmer = current_kmer[1:].join(sequence[i + kmer_size - 1])
        if not (remove_wildcard and any(
                nucleotide in WILDCARD_READINGS for nucleotide in
                current_kmer)):
            yield current_kmer, i


class KmerReference:

    def __init__(self, kmer_size: int):
        self._kmer_db: Dict[str, Dict[str, List[int]]] = {}
        self._genomes_db: Dict[str, ReferencedGenomeStats] = {}
        self._kmer_size = kmer_size

    def add_kmers_to_db(self, genome: RawGenome) -> None:
        for kmer in extract_kmers_from_string(genome.sequence, self._kmer_size,
                                              remove_wildcard=True):
            if kmer[0] not in self._kmer_db:
                self._kmer_db[kmer[0]] = {}
            if genome.identifier not in self._kmer_db[kmer[0]]:
                self._kmer_db[kmer[0]][genome.identifier] = [kmer[1]]
            else:
                self._kmer_db[kmer[0]][genome.identifier].append(kmer[1])


    def build_kmer_reference(self, genome_file: str) -> bool:
        """
        This function builds is a kmer-reference builder.
        It receives the path to a fasta file and uses the file parser and kmer extractor
        to build the kmer DB reference
        :param genome_file: the path to the fasta file
        :return: True if the build was successful, False otherwise
        """
        try:
            for genome in file_parsers.parse_fasta_file(genome_file):
                validators.validate_values_in_given_list(genome.sequence,
                                                         ALLOWED_DNA_VALUES)
                if genome.identifier in self._genomes_db:
                    # TODO - add case when genome appears twice in fasta
                    pass
                else:
                    self._genomes_db[genome.identifier] = ReferencedGenomeStats(
                        genome.identifier)
                    self.add_kmers_to_db(genome)
        except Exception as e:
            return False
