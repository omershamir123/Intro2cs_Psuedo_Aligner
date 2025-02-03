# FILE : kmer_reference.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: A file with the KmerReference Class
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:

import json
from typing import Dict, List, Generator, Tuple

from genome import ReferencedGenome
from program_constants import WILDCARD_READINGS, UNIQUE_KMER, MULTI_MAP_KMER


def extract_kmers_from_string(sequence: str, kmer_size: int,
                              remove_wildcard: bool = True) -> \
        Generator[Tuple[str, int], None, None]:
    """
    This function receives a sequence and returns all kmers within that sequence.
    If a kmer contains a WILDCARD_READING - then the kmer will be thrown out
    :param sequence: the whole DNA sequence
    :param kmer_size: the length of the kmer to be extracted
    :param remove_wildcard: should we remove the kmers containing a wildcard
    :return: Tuple with the first parameter being the kmer, the second being its position in the initial sequence
    """
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
        self._similarity_results: dict = {}

    @property
    def kmer_db(self) -> Dict[str, Dict[str, List[int]]]:
        return self._kmer_db

    @property
    def genomes_db(self) -> Dict[str, ReferencedGenome]:
        return self._genomes_db

    @property
    def kmer_size(self) -> int:
        return self._kmer_size

    def check_reference_was_filtered(self)->bool:
        return self._similarity_results != {}

    def add_kmers_to_db(self, genome: ReferencedGenome) -> None:
        """
        This function populates the KmerDB and updates the genomes DB based on the given genomes kmers
        It checks for each kmer in the genome whether it is its first appearance in the reference and adds it
        If it's the first time the specific kmer appears in the current genome, it updates the unique and ambiguous kmers counts
        for all genomes with that kmer within their sequence
        :param genome: The current genome to map
        :return: None
        """
        for kmer in extract_kmers_from_string(genome.sequence, self._kmer_size):
            current_kmer, current_position = kmer
            genome.add_kmer_to_genome_mapping(current_kmer)
            if current_kmer not in self._kmer_db:
                self._kmer_db[current_kmer] = {}
            if genome.identifier not in self._kmer_db[current_kmer]:
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
        from file_handlers import parse_fasta_file
        try:
            for genome in parse_fasta_file(genome_file):
                if genome.identifier in self._genomes_db:
                    raise ValueError(
                        "There are two genomes in the fasta file with the header {}".format(
                            genome.identifier))
                else:
                    self._genomes_db[genome.identifier] = genome
                    self.add_kmers_to_db(genome)
        except Exception as e:
            print(e)
            return False
        return True

    def calculate_kmers_type(self) -> None:
        for genome in self._genomes_db.values():
            genome.unique_kmers = 0
            genome.multi_mapping_kmers = 0
            for kmer in genome.kmers_set:
                if len(self._kmer_db[kmer]) == 1:
                    genome.unique_kmers += len(
                        self._kmer_db[kmer][genome.identifier])
                else:
                    genome.multi_mapping_kmers += len(
                        self._kmer_db[kmer][genome.identifier])

    def genome_db_to_dict(self) -> Dict[str, Dict[str, int]]:
        return {k: v.genome_ref_to_dict() for k, v in self.genomes_db.items()}

    def to_json(self) -> str:
        return json.dumps(
            {"Kmers": self._kmer_db, "Summary": self.genome_db_to_dict()})

    def filter_reference_based_on_similarity(self,
                                             similarity_threshold: float) -> \
            Dict[
                str, Tuple[str, float]]:
        """
        This function filters out genomes from the reference that are similar to others
        If two genomes are similar, the function also makes sure the reference is updated accordingly
        :param similarity_threshold: the similarity threshold to filter out genomes
        :return: a dictionary - the keys are the removed genome identifiers and the values are the genome that they are similar to
        """
        filtered_genomes: Dict[str, Tuple[str, float]] = {}
        # We initialize the genome_list separately because while iterating over the
        # genomes and filtering similar ones, we will update the unique + multi_map_kmers
        genome_list = [genome for genome in self.genomes_db]
        genome_list = sorted(genome_list, key=lambda genome: (
            self._genomes_db[genome].unique_kmers,
            self._genomes_db[genome].multi_mapping_kmers,
            self._genomes_db[genome].total_bases,
            self._genomes_db[genome].index_in_fasta))

        for i in range(len(genome_list)):
            current_genome = genome_list[i]
            current_genome_kmers = self._genomes_db[
                current_genome].kmers_set
            for j in range(i + 1, len(genome_list)):
                other_genome = genome_list[j]
                other_genome_kmers = self._genomes_db[
                    other_genome].kmers_set
                # case to avoid edge case of genome with 0 kmers
                if len(current_genome_kmers) != 0 and len(
                        other_genome_kmers) != 0:
                    similarity_score = len(
                        current_genome_kmers & other_genome_kmers) / min(
                        len(current_genome_kmers), len(other_genome_kmers))
                    if similarity_score > similarity_threshold:
                        self.remove_genome_by_similarity(current_genome)
                        filtered_genomes[current_genome] = (
                            other_genome, similarity_score)
                        break

        return filtered_genomes

    def remove_genome_by_similarity(self, removed_genome_id: str) -> None:
        """
        This function removed an entire genome from the kmer_db reference
        and updating the unique kmers for all related genomes
        :param removed_genome_id: the genome to remove
        :return: None
        """
        removed_genome = self._genomes_db[removed_genome_id]

        for kmer in removed_genome.kmers_set:
            if len(self._kmer_db[kmer]) == 1:
                self._kmer_db.pop(kmer)
            else:
                self._kmer_db[kmer].pop(removed_genome_id)

    def similarity_json(self, filtered_genomes: Dict[str, Tuple[str, float]]):
        similar_dict = {}
        for genome in self._genomes_db.values():
            genome_values = {"kept": "yes", "unique_kmers": genome.unique_kmers,
                             "total_kmers": genome.unique_kmers + genome.multi_mapping_kmers,
                             "genome_length": genome.total_bases,
                             "similar_to": "NA", "similarity_score": "NA"}
            if genome.identifier in filtered_genomes:
                genome_values["similar_to"] = \
                    filtered_genomes[genome.identifier][0]
                genome_values["similarity_score"] = \
                    filtered_genomes[genome.identifier][1]
                genome_values["kept"] = "no"
            similar_dict[genome.identifier] = genome_values
        for filtered_genome in filtered_genomes:
            self.genomes_db.pop(filtered_genome)
        return {"Similarity": similar_dict}

    def filter_genomes_logic(self, similarity_threshold: float) -> None:
        filtered_genomes = self.filter_reference_based_on_similarity(
            similarity_threshold)
        self._similarity_results = self.similarity_json(filtered_genomes)

    def print_similarity_results(self) -> str:
        return json.dumps(self._similarity_results)
