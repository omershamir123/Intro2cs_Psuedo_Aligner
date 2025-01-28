# FILE : pseudo_aligner.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: This file includes the pseudo aligner class that handles the pseudo_aligner algorithm
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
from typing import Dict, List, Tuple, Optional
import file_handlers
from kmer_reference import KmerReference, extract_kmers_from_string
from program_constants import UNMAPPED_READ, UNIQUE_READ, READ_STATUS, \
    AMBIGUOUS_READ
from read import Read, ReadKmerMapping


class PseudoAligner:
    def __init__(self):
        self._reads: Dict[str, Read] = {}
        self._unique_mapped_reads = 0
        self._ambiguous_mapped_reads = 0
        self._unmapped_reads = 0
        self._filtered_quality_reads = 0
        self._filtered_quality_kmers = 0
        self._filtered_hr_kmers = 0

    def align_algorithm(self, fastq_file_path: str,
                        kmer_reference: KmerReference, unique_threshold: int,
                        ambiguous_threshold: int) -> bool:
        """
        This function runs the pseudo_align algorithm as given in the instructions.
        For each read in the fastq_file:
        1. Process all kmers and classify each one as specific or non_specific
        2. Map the read to a specific genome based on the specific kmers
        3. Validate the mapping based on the unspecific kmers
        :param fastq_file_path: the path to the fastq file
        :param kmer_reference: a KmerReference object
        :param unique_threshold: threshold to distinguish between the highest specific genome to its second highest in the read
        :param ambiguous_threshold: threshold to distinguish two ambiguously mapped genomes
        :return: None
        """
        kmer_size = len(next(iter(kmer_reference.kmer_db))[0])
        try:
            for read in file_handlers.parse_fastq_file(fastq_file_path):
                if read.identifier in self._reads:
                    raise ValueError("Duplicate reads detected")
                self._reads[read.identifier] = read
                current_read_mapping = self.extract_and_map_kmers_from_read(
                    read,
                    kmer_reference,
                    kmer_size)
                if len(current_read_mapping.specific_kmers) == 0:
                    self.map_read(read, UNMAPPED_READ)
                else:
                    self._map_read_using_specific_kmers(current_read_mapping,
                                                        read, kmer_reference,
                                                        unique_threshold)
                    if read.read_status == UNIQUE_READ:
                        self._validate_uniquness_using_unspecific(
                            current_read_mapping, read, kmer_reference,
                            ambiguous_threshold)

        except Exception as e:
            print(e)
            return False
        return True

    def map_read(self, read, status: READ_STATUS,
                 kmer_reference: KmerReference = None,
                 genome_identifier: str = None) -> None:
        """
        This function receives a read and sets its status as well as the effect of
        that status on all relevant genomes
        :param read: the read to mark
        :param status: the status of the read
        :param kmer_reference: the kmer reference - None if UNMAPPED_READ
        :param genome_identifier: the genome identifier - None if UNMAPPED_READ
        :return: None
        """
        read.read_status = status
        read.add_mapped_genome(genome_identifier)
        if status == UNIQUE_READ:
            kmer_reference.genomes_db[genome_identifier].unique_reads += 1
        if status == AMBIGUOUS_READ:
            kmer_reference.genomes_db[genome_identifier].ambiguous_reads += 1

    def extract_and_map_kmers_from_read(self, read: Read,
                                        kmer_reference: KmerReference,
                                        kmer_size: int) -> ReadKmerMapping:
        """
        This function extracts the kmers from the read and maps each kmer as specific, unspecific or neither
        It does that based on each kmer's appearance in the kmer_reference
        :param read: the current read
        :param kmer_reference: the KmerReference object
        :param kmer_size: the kmer_size
        :return: The current read's kmer mapping object
        """
        current_read_kmer_classification = ReadKmerMapping()
        for kmer_tuple in extract_kmers_from_string(read.value, kmer_size,
                                                    True):
            kmer = kmer_tuple[0]
            # Check if the current kmer from the read is in the given reference
            if kmer in kmer_reference.kmer_db:
                # If only one genome has this kmer - it's a specific one
                if len(kmer_reference.kmer_db[kmer]) == 1:
                    self._add_specific_kmer_to_classification(kmer,
                                                              kmer_reference,
                                                              current_read_kmer_classification)
                # Unspecific Kmer
                else:
                    self._add_unspecific_kmer_to_classification(kmer,
                                                                kmer_reference,
                                                                current_read_kmer_classification)
        return current_read_kmer_classification

    def _add_specific_kmer_to_classification(self, kmer: str,
                                             kmer_reference: KmerReference,
                                             current_read_kmer_classification: ReadKmerMapping) -> None:
        """
        This function adds a specific kmer to a readMapping object
        Since it's a specific kmer - in the kmer reference only one genome is mapped to it
        hence we can extract the genome identifier in O(1)
        :param kmer: the specific kmer to be mapped
        :param kmer_reference: the kmer reference
        :param current_read_kmer_classification: the mapping object for a single read
        :return: None
        """
        genome_identifier = \
            next(iter(kmer_reference.kmer_db[kmer]))[0]
        if genome_identifier not in current_read_kmer_classification.specific_kmers_in_genomes:
            current_read_kmer_classification.specific_kmers_in_genomes[
                genome_identifier] = [kmer]
        else:
            current_read_kmer_classification.specific_kmers_in_genomes[
                genome_identifier].append(kmer)
        current_read_kmer_classification.specific_kmers.append(kmer)

    def _add_unspecific_kmer_to_classification(self, kmer: str,
                                               kmer_reference: KmerReference,
                                               current_read_kmer_classification: ReadKmerMapping) -> None:
        """
        This function adds an uspecific kmer to a readMapping object
        Since it's an uspecific kmer - in the kmer reference several genomes are mapped to it
        hence we will iterate over the genomes and update the mapping for each genome
        :param kmer: the specific kmer to be mapped
        :param kmer_reference: the kmer reference
        :param current_read_kmer_classification: the mapping object for a single read
        :return: None
        """
        for genome in kmer_reference.kmer_db[kmer]:
            genome_identifier = genome[0]
            if genome_identifier not in current_read_kmer_classification.unspecific_kmers_in_genomes:
                current_read_kmer_classification.unspecific_kmers_in_genomes[
                    genome_identifier] = [kmer]
            else:
                current_read_kmer_classification.unspecific_kmers_in_genomes[
                    genome_identifier].append(kmer)
        current_read_kmer_classification.unspecific_kmers.append(kmer)

    def _map_read_using_specific_kmers(self,
                                       current_read_mapping: ReadKmerMapping,
                                       read: Read,
                                       kmer_reference: KmerReference,
                                       unique_threshold: int) -> None:
        """
        This is the core part of the pseudo_aligner algorithm - after mapping all kmers
        It checks the two genome with the highest specific-kmers found in the read
        If the difference between them is greater than or equal to unique_threshold, maps the read
        uniquely to the genome with most specific kmers in this read
        Otherwise, ambiguously maps the read and increments the unique counter for each genome
        that has a specific kmer for this read
        :param current_read_mapping: the mapping for the read - specific and unspecific kmers
        :param read: the current read
        :param kmer_reference: the kmer reference to access the genomes
        :param unique_threshold: the minimum difference between the two highest appearing genomes
        :return: None
        """
        if len(current_read_mapping.specific_kmers_in_genomes) == 1:
            genome_identifier = \
                next(iter(current_read_mapping.specific_kmers_in_genomes))[0]
            self.map_read(read, UNIQUE_READ, kmer_reference, genome_identifier)
        else:
            frequent, second_frequent = self.find_two_most_frequent_genomes(
                current_read_mapping.specific_kmers_in_genomes)
            frequency_difference = frequent[1] - second_frequent[1]
            if frequency_difference >= unique_threshold:
                genome_identifier = frequent[0]
                self.map_read(read, UNIQUE_READ, kmer_reference,
                              genome_identifier)
            else:
                for genome in current_read_mapping.specific_kmers_in_genomes:
                    genome_identifier = genome[0]
                    self.map_read(read, AMBIGUOUS_READ, kmer_reference,
                                  genome_identifier)

    def find_two_most_frequent_genomes(self, specific_kmers_in_genomes: Dict[
        str, List[str]]) -> Tuple[Tuple[str, int], Tuple[str, int]]:
        """
        This function receives a dictionary of genomes and the specific kmers that matched them.
        It scans through the dictionary and returns the two highest genomes with most kmers
        :param specific_kmers_in_genomes: the dictionary of genomes and their specific kmers
        :return: Tuple of (first,second) each being a Tuple of (genome_identifier,count)
        """
        first: Optional[Tuple[str, int]] = None
        second: Optional[Tuple[str, int]] = None
        for genome_identifier, specific_kmer_list in specific_kmers_in_genomes.items():
            if first is None or len(specific_kmer_list) >= first[1]:
                first, second = (
                    genome_identifier, len(specific_kmer_list)), first
            elif second is None or len(specific_kmer_list) >= second[1]:
                second = (genome_identifier, len(specific_kmer_list))

        return first, second

    def _validate_uniquness_using_unspecific(self,
                                             current_read_mapping: ReadKmerMapping,
                                             read: Read,
                                             kmer_reference: KmerReference,
                                             ambiguous_threshold: int) -> None:

        specific_kmers_dict = current_read_mapping.specific_kmers_in_genomes
        unspecific_kmers_dict = current_read_mapping.unspecific_kmers_in_genomes

        total_counter_dict = {key: specific_kmers_dict.get(key, []).extend(
            unspecific_kmers_dict.get(key, [])) for key in
                              set(specific_kmers_dict) | set(
                                  unspecific_kmers_dict)}

        max_count = max(total_counter_dict.values())
        # Our read was uniquely mapped - hence the mapped_genomes to it is only one
        unique_genome_identifier = next(iter(read.mapped_genomes))
        map_count = total_counter_dict[unique_genome_identifier]
        if max_count - map_count > ambiguous_threshold:
            # ambiguously map all genomes
            for genome_identifier, total_count in total_counter_dict.items():
                if total_count >= map_count:
                    self.map_read(read, AMBIGUOUS_READ, kmer_reference, genome_identifier)

