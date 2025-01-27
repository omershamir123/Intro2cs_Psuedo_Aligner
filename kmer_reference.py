# FILE : kmer_reference.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: A file with the KmerReference Class
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
import json
from typing import Dict, List, Generator, Tuple

import file_handlers
import validators
from genome import ReferencedGenome
from program_constants import ALLOWED_DNA_VALUES, WILDCARD_READINGS


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
    # if kmer_size > len(sequence):
    #     raise StopIteration
    # else:
    if kmer_size <= len(sequence):
        current_kmer = sequence[:kmer_size]
        if not (remove_wildcard and any(
                nucleotide in WILDCARD_READINGS for nucleotide in
                current_kmer)):
            yield current_kmer, 0
    for i in range(1, len(sequence) - kmer_size + 1):
        current_kmer = current_kmer[1:] + sequence[i + kmer_size - 1]
        if not (remove_wildcard and any(
                nucleotide in WILDCARD_READINGS for nucleotide in
                current_kmer)):
            yield current_kmer, i


class KmerReference:

    def __init__(self, kmer_size: int):
        self._kmer_db: Dict[str, Dict[str, List[int]]] = {}
        self._genomes_db: Dict[str, ReferencedGenome] = {}
        self._kmer_size = kmer_size

    @property
    def kmer_db(self) -> Dict[str, Dict[str, List[int]]]:
        return self._kmer_db

    @property
    def genomes_db(self) -> Dict[str, ReferencedGenome]:
        return self._genomes_db

    def add_kmers_to_db(self, genome: ReferencedGenome) -> None:
        for kmer in extract_kmers_from_string(genome.sequence, self._kmer_size,
                                              remove_wildcard=True):
            current_kmer, current_position = kmer
            if current_kmer not in self._kmer_db:
                self._kmer_db[current_kmer] = {}
            if genome.identifier not in self._kmer_db[current_kmer]:
                # check if it's the first genome added - if so, it's unique
                if len(self._kmer_db[current_kmer].keys()) == 0:
                    self._genomes_db[genome.identifier].unique_kmers += 1
                # if there is another genome already associated with this kmer,
                # decrease the unique kmer count and increase the multi_mapping
                elif len(self._kmer_db[current_kmer].keys()) == 1:
                    other_genome_identifier = next(
                        iter(self._kmer_db[current_kmer]))
                    self._genomes_db[other_genome_identifier].unique_kmers -= 1
                    self._genomes_db[
                        other_genome_identifier].multi_mapping_kmers += 1
                # if there is more than one genome associated with this kmer
                # add the multi mapping kmer to this genome
                if len(self._kmer_db[current_kmer].keys()) >= 1:
                    self._genomes_db[
                        genome.identifier].multi_mapping_kmers += 1
                self._kmer_db[current_kmer][genome.identifier] = [
                    current_position]

            else:
                self._kmer_db[current_kmer][genome.identifier].append(
                    current_position)

    def build_kmer_reference(self, genome_file: str) -> bool:
        """
        This function builds is a kmer-reference builder.
        It receives the path to a fasta file and uses the file parser and kmer extractor
        to build the kmer DB reference
        :param genome_file: the path to the fasta file
        :return: True if the build was successful, False otherwise
        """
        try:
            for genome in file_handlers.parse_fasta_file(genome_file):
                validators.validate_values_in_given_list(genome.sequence,
                                                         ALLOWED_DNA_VALUES)
                if genome.identifier in self._genomes_db:
                    return False
                else:
                    self._genomes_db[genome.identifier] = genome
                    self.add_kmers_to_db(genome)
        except Exception as e:
            return False

        return True


    def genome_db_to_json(self) -> Dict[str,Dict[str, List[int]]]:
        return {k:v.genome_ref_to_dict() for k, v in self.genomes_db.items()}

    def to_json(self):
        return json.dumps({"Kmers":self._kmer_db, "Summary": self.genome_db_to_json()})