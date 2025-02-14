# FILE : pseudo_aligner.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: This file includes the pseudo aligner class that handles the pseudo_aligner algorithm
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:

import numpy as np

import json
from typing import Dict, List, Tuple, Optional
from kmer_reference import KmerReference, extract_kmers_from_string
from program_constants import UNMAPPED_READ, UNIQUE_READ, READ_STATUS, \
    AMBIGUOUS_READ
from read import Read, ReadKmerMapping


class PseudoAlignerOutput:
    """
    This class is an object that holds data regarding a specific pseudo align algorithm run
    """

    def __init__(self, kmer_reference: KmerReference,
                 coverage_included: bool = False,
                 genome_list_str: Optional[str] = None):
        self._reads: Dict[str, Read] = {}
        self._kmer_reference: KmerReference = kmer_reference
        self._unique_mapped_reads = 0
        self._ambiguous_mapped_reads = 0
        self._unmapped_reads = 0
        self._filtered_quality_reads = -1
        self._filtered_quality_kmers = -1
        self._filtered_hr_kmers = -1
        self._genome_list: List[str] = extract_genomes_list(genome_list_str,
                                                            kmer_reference) if coverage_included else []

    @property
    def reads(self) -> Dict[str, Read]:
        return self._reads

    @property
    def kmer_reference(self) -> KmerReference:
        return self._kmer_reference

    @property
    def unique_mapped_reads(self) -> int:
        return self._unique_mapped_reads

    @property
    def ambiguous_mapped_reads(self) -> int:
        return self._ambiguous_mapped_reads

    @property
    def unmapped_reads(self) -> int:
        return self._unmapped_reads

    @property
    def filtered_quality_reads(self) -> int:
        return self._filtered_quality_reads

    @filtered_quality_reads.setter
    def filtered_quality_reads(self, value: int):
        self._filtered_quality_reads = value

    @property
    def filtered_quality_kmers(self) -> int:
        return self._filtered_quality_kmers

    @filtered_quality_kmers.setter
    def filtered_quality_kmers(self, value: int):
        self._filtered_quality_kmers = value

    @property
    def filtered_hr_kmers(self) -> int:
        return self._filtered_hr_kmers

    @filtered_hr_kmers.setter
    def filtered_hr_kmers(self, value: int):
        self._filtered_hr_kmers = value

    @property
    def genome_list(self) -> List[str]:
        return self._genome_list

    def add_read(self, read: Read) -> None:
        if not read.identifier in self._reads:
            self._reads[read.identifier] = read

    def update_reads_stats(self, read_identifier: str) -> None:
        """
        This function checks the status of a given read and updates the total count
        of ambiguous, unique, unmapped reads accordingly
        :param read_identifier: the ID of the read to look for in the dictionary
        :return: None
        """
        if read_identifier in self._reads:
            status = self._reads[read_identifier].read_status
            if status == UNMAPPED_READ:
                self._unmapped_reads += 1
            elif status == AMBIGUOUS_READ:
                self._ambiguous_mapped_reads += 1
            elif status == UNIQUE_READ:
                self._unique_mapped_reads += 1

    def convert_to_aln_object(self, is_reversed: bool) -> "AlnFileDataObject":
        return AlnFileDataObject(self, is_reversed=is_reversed)


class AlnFileDataObject:
    """
    This class holds the relevant data from the result of a pseudo align algorithm run.
    It is the object that is being put inside an ALN file
    """

    def __init__(self, aligner_output: PseudoAlignerOutput, is_reversed: bool):
        self.genomes_db = aligner_output.kmer_reference.genomes_db
        self._unique_mapped_reads = aligner_output.unique_mapped_reads
        self._ambiguous_mapped_reads = aligner_output.ambiguous_mapped_reads
        self._unmapped_reads = aligner_output.unmapped_reads
        self._filtered_quality_reads = aligner_output.filtered_quality_reads
        self._filtered_quality_kmers = aligner_output.filtered_quality_kmers
        self._filtered_hr_kmers = aligner_output.filtered_hr_kmers
        self._reverse_complement_included = is_reversed
        self._genome_list = aligner_output.genome_list
        self._coverage_results: dict = {}

    def check_coverage_was_applied(self) -> bool:
        return self._coverage_results != {}

    def genomes_mapped_to_dict(self) -> Dict[str, Dict[str, int]]:
        if not self._reverse_complement_included:
            return {k: v.genome_mapped_to_dict() for k, v in
                    self.genomes_db.items()}
        else:
            return {k + "_F": v.genome_mapped_to_dict(is_reversed=False) for
                    k, v in
                    self.genomes_db.items()} | {
                k + "_R": v.genome_mapped_to_dict(is_reversed=True) for k, v in
                self.genomes_db.items()}

    def read_stats_to_dict(self) -> Dict[str, int]:
        summary = {"unique_mapped_reads": self._unique_mapped_reads,
                   "ambiguous_mapped_reads": self._ambiguous_mapped_reads,
                   "unmapped_reads": self._unmapped_reads}
        if self._filtered_quality_reads >= 0:
            summary["filtered_quality_reads"] = self._filtered_quality_reads
        if self._filtered_quality_kmers >= 0:
            summary["filtered_quality_kmers"] = self._filtered_quality_kmers
        if self._filtered_hr_kmers >= 0:
            summary["filtered_hr_kmers"] = self._filtered_hr_kmers
        return summary

    def to_json(self):
        return json.dumps({"Statistics": self.read_stats_to_dict(),
                           "Summary": self.genomes_mapped_to_dict()})

    def coverage_statistics_summary(self, min_coverage: int):
        return {
            genome: self.genomes_db[genome].genome_coverage_stats(min_coverage)
            for genome in
            self._genome_list}

    def full_coverage_stats(self):
        return {genome: self.genomes_db[genome].genome_full_coverage_stats() for
                genome in
                self._genome_list}

    def coverage_logic(self, apply_full_coverage: bool, min_coverage: int):
        if not apply_full_coverage:
            self._coverage_results = {
                "Coverage": self.coverage_statistics_summary(min_coverage)}
        else:
            self._coverage_results = {
                "Coverage": self.coverage_statistics_summary(min_coverage),
                "Details": self.full_coverage_stats()}

    def print_coverage_results(self):
        return json.dumps(self._coverage_results)


def should_filter_read(read: Read, min_read_quality: int) -> bool:
    """
    This function checks for a specific read if it should be filtered based on low quality
    :param read: the read to check
    :param min_read_quality: the minimum quality threshold
    :return: True if read should be filtered, False otherwise
    """
    if read.calculate_mean_quality() < min_read_quality:
        return True
    return False


def try_mapping_using_specific_kmers(current_read_mapping,
                                     current_unique_threshold: int) -> Tuple[
    READ_STATUS, int]:
    """
    This function is used in the reverse complement extension and checks without mapping
    a read what would be its status after scanning through its specific kmers
    :param current_unique_threshold:
    :param current_read_mapping:
    :return: A tuple with the read status and an integer that contains:
    For UNIQUE_READ: the number of specific kmers for the genome mapped to it
    For AMBIGUOUS_READ: the max number of total_kmers for the genomes mapped to it
    For UNMAPPED_READ: -1
    """
    if len(current_read_mapping.specific_kmers) == 0:
        return UNMAPPED_READ, -1
    if len(current_read_mapping.specific_kmers_in_genomes) == 1:
        genome_identifier = \
            next(iter(current_read_mapping.specific_kmers_in_genomes))
        return UNIQUE_READ, len(
            current_read_mapping.specific_kmers_in_genomes[genome_identifier])

    # By this point - there are at least two genomes that have specific kmers from this read
    # hence, frequent, second_frequent won't return None
    frequent, second_frequent = find_two_most_frequent_genomes(
        current_read_mapping.specific_kmers_in_genomes)
    frequency_difference = frequent[1] - second_frequent[1]
    if frequency_difference >= current_unique_threshold:
        return UNIQUE_READ, frequent[1]
    else:
        # the mapping is ambiguous, hence we need to calculate the max total kmers of a single genome for this read
        total_kmers_in_genomes_dict = {
            genome: current_read_mapping.specific_kmers_in_genomes.get(genome,
                                                                    []) + current_read_mapping.unspecific_kmers_in_genomes.get(
                genome, []) for genome in
            current_read_mapping.specific_kmers_in_genomes.keys() | current_read_mapping.unspecific_kmers_in_genomes.keys()}

        return AMBIGUOUS_READ, max(len(kmer_list) for kmer_list in total_kmers_in_genomes_dict.values())


def should_read_be_reversed(
        read_status_if_mapped_forward: Tuple[READ_STATUS, int],
        read_status_if_mapped_reversed: Tuple[READ_STATUS, int]) -> bool:
    """
    This function is used in the reverse complement extension and receives the
    mapping results of both orientations.
    It determines via a set of cases the best orientation to prefer
    :param read_status_if_mapped_forward: Tuple containing the forward read status and it matching number
    :param read_status_if_mapped_reversed: Tuple containing the reverse read status and it matching number
    :return: True if reverse is preferable, False if forward is preferable
    """
    forward_status = read_status_if_mapped_forward[0]
    forward_number = read_status_if_mapped_forward[1]
    reverse_status = read_status_if_mapped_reversed[0]
    reverse_number = read_status_if_mapped_reversed[1]
    if forward_status == UNIQUE_READ and reverse_status == UNIQUE_READ:
        return reverse_number > forward_number
    if forward_status == UNIQUE_READ:
        return False
    if reverse_status == UNIQUE_READ:
        return True

    if forward_status == AMBIGUOUS_READ and reverse_status == AMBIGUOUS_READ:
        return reverse_number > forward_number
    # The reverse isn't unique, it's unmapped - hence return the forward is preferable
    if reverse_status == AMBIGUOUS_READ:
        return True
    # By now - the reverse must be unmapped - hence return the forward orientation
    return False


def determine_best_mapping_for_read(read: Read,
                                    aligner_output: PseudoAlignerOutput,
                                    kmer_reference: KmerReference,
                                    current_unique_threshold: int,
                                    **kwargs) -> ReadKmerMapping:
    """
    This function calculates the KmerMapping object of the current read.
    If the reverse_complement flag is inserted - it determines the better orientaion as well
    It does so by calculating the amount of specific kmers in each orientation
    If the reverse_complement flag isn't inserted, it returns the default forward orientation read kmer mapping

    :param current_unique_threshold: the unique threshold for the alignment
    :param read: the current read
    :param aligner_output: the aligner_output object
    :param kmer_reference: the kmer_reference object
    :param kwargs: more optional arguments given by the user - one being reverse_complement
    :return: the kmer mapping of this read, if needed, the reverse orientation mapping
    """
    filtered_quality_kmers_before_read_check = aligner_output.filtered_quality_kmers
    filtered_hr_kmers_before_read_check = aligner_output.filtered_hr_kmers
    check_reverse_complement = kwargs.get("reverse_complement")
    forward_read_mapping = extract_and_map_kmers_from_read(
        read,
        kmer_reference,
        aligner_output, **kwargs)
    filtered_hr_kmers_after_forward_check = aligner_output.filtered_hr_kmers
    filtered_quality_kmers_after_forward_check = aligner_output.filtered_quality_kmers
    current_read_mapping = forward_read_mapping
    if check_reverse_complement:
        read_status_if_mapped_forward = try_mapping_using_specific_kmers(
            current_read_mapping, current_unique_threshold)
        read.is_reversed = True
        reversed_read_mapping = extract_and_map_kmers_from_read(
            read,
            kmer_reference,
            aligner_output, **kwargs)
        read_status_if_mapped_reversed = try_mapping_using_specific_kmers(
            reversed_read_mapping, current_unique_threshold)
        read_to_be_reversed = should_read_be_reversed(
            read_status_if_mapped_forward, read_status_if_mapped_reversed)
        # In this case, and this case only - prefer the reverse mapping
        if read_to_be_reversed:
            current_read_mapping = reversed_read_mapping
            hr_kmers_diff = aligner_output.filtered_hr_kmers - filtered_hr_kmers_after_forward_check
            aligner_output.filtered_hr_kmers = filtered_hr_kmers_before_read_check + hr_kmers_diff

            quality_kmers_diff = aligner_output.filtered_quality_kmers - filtered_quality_kmers_after_forward_check
            aligner_output.filtered_quality_kmers = filtered_quality_kmers_before_read_check + quality_kmers_diff

        else:
            read.is_reversed = False
            aligner_output.filtered_hr_kmers = filtered_hr_kmers_after_forward_check
            aligner_output.filtered_quality_kmers = filtered_quality_kmers_after_forward_check

    return current_read_mapping


def initialize_genome_coverage(genome_list: List[str],
                               kmer_reference: KmerReference) -> bool:
    """
    This function is called once the coverage extension is activated, and it initializes the coverage
    array for each genome in the genome_list
    :param genome_list: the list of genomes to apply the coverage
    :param kmer_reference: the kmer reference which contains the genomes db
    :return: True if the initialization is successful, False otherwise (a genome in the genome list isn't in the db)
    """
    for genome_identifier in genome_list:
        if genome_identifier not in kmer_reference.genomes_db:
            print(
                "There was an error in the coverage genomes list. The genome {} isn't in the reference".format(
                    genome_identifier))
            return False
        kmer_reference.genomes_db[
            genome_identifier].initialize_coverage_arrays()
    return True


def update_genome_coverage(read: Read, current_read_mapping: ReadKmerMapping,
                           kmer_reference: KmerReference,
                           genome_list: List[str]) -> None:
    """
    This function is called once the coverage extension is activated, and it updates the genome coverage
    for all genomes in the genome_list for the current read.
    It goes over all the mapped genomes of the read - If they are in genome_list - then if does the following:
    Goes over all the specific kmers + unspecific_kmers in the ReadMap - if they are present in the genome
    It does so by looking for the genome in the kmer_reference.kmer_db
    :param genome_list: the list of genomes for the coverage
    :param read: the current read
    :param current_read_mapping: the mapping of the current read
    :param kmer_reference: the kmer reference which contains the kmer_db
    :return: None
    """

    for genome_identifier in read.mapped_genomes:
        if genome_identifier in genome_list:
            covered_bases = np.zeros(
                kmer_reference.genomes_db[genome_identifier].total_bases,
                dtype=int)
            kmers_in_genome_present_in_read = \
                current_read_mapping.specific_kmers_in_genomes[
                    genome_identifier] + \
                current_read_mapping.unspecific_kmers_in_genomes.get(
                    genome_identifier, [])

            for kmer in kmers_in_genome_present_in_read:
                for starting_position in kmer_reference.kmer_db[kmer][
                    genome_identifier]:
                    covered_bases[
                    starting_position:starting_position + kmer_reference.kmer_size] = 1

            if read.read_status == UNIQUE_READ:
                kmer_reference.genomes_db[
                    genome_identifier].unique_coverage_positions += covered_bases

            if read.read_status == AMBIGUOUS_READ:
                kmer_reference.genomes_db[
                    genome_identifier].unique_coverage_positions += covered_bases


def extract_genomes_list(genome_list_str: Optional[str],
                         kmer_reference: KmerReference) -> List[str]:
    if genome_list_str is None:
        genome_list = list(kmer_reference.genomes_db.keys())
    else:
        genome_list = str.split(genome_list_str, ",")
    return genome_list


def initialize_filtering_stats(aligner_output: PseudoAlignerOutput, **kwargs):
    """
    This function initializes the filtering stats in the aligner_output
    :param aligner_output: the aligner output
    :param kwargs: the arguments given by the user that include the filtering flag
    :return: None
    """
    min_read_quality = kwargs.get("min_read_quality")
    if min_read_quality is not None:
        aligner_output.filtered_quality_reads = 0
    min_kmer_quality = kwargs.get("min_kmer_quality")
    if min_kmer_quality is not None:
        aligner_output.filtered_quality_kmers = 0
    max_genomes = kwargs.get("max_genomes")
    if max_genomes is not None:
        aligner_output.filtered_hr_kmers = 0


def align_algorithm(fastq_file_path: str,
                    kmer_reference: KmerReference,
                    current_unique_threshold: int,
                    ambiguous_threshold: int, **kwargs) -> Optional[
    PseudoAlignerOutput]:
    """
    This function runs the pseudo_align algorithm as given in the instructions.
    For each read in the fastq_file:
    1. Process all kmers and classify each one as specific or non_specific
    2. Map the read to a specific genome based on the specific kmers
    3. Validate the mapping based on the unspecific kmers
    :param fastq_file_path: the path to the fastq file
    :param kmer_reference: a KmerReference object
    :param current_unique_threshold: threshold to distinguish between the highest specific genome to its second highest in the read
    :param ambiguous_threshold: threshold to distinguish two ambiguously mapped genomes
    :return: None
    """
    # mypy might be unhappy but in the argsparser, the value is store_true
    check_coverage = kwargs.get("coverage")
    genome_list_str = kwargs.get("genomes")
    # mypy might be unhappy but there is an action in the argsparse - store_true
    aligner_output = PseudoAlignerOutput(kmer_reference, check_coverage,
                                         genome_list_str)
    initialize_filtering_stats(aligner_output, **kwargs)
    # mypy might be unhappy but in the argsparser, there is a default value
    min_read_quality = kwargs.get("min_read_quality")
    filter_read_by_quality = min_read_quality is not None
    check_coverage = kwargs.get("coverage")
    genome_list = aligner_output.genome_list
    if check_coverage:
        if not initialize_genome_coverage(genome_list, kmer_reference):
            return None
    from file_handlers import parse_fastq_file
    try:
        for read in parse_fastq_file(fastq_file_path):
            if read.identifier in aligner_output.reads:
                raise ValueError("Duplicate reads detected")
            # Check whether the quality of the entire read is sufficient
            # same here for mypy...
            if filter_read_by_quality and should_filter_read(read,
                                                             min_read_quality):
                aligner_output.filtered_quality_reads += 1
                continue
            aligner_output.add_read(read)

            current_read_mapping = determine_best_mapping_for_read(read,
                                                                   aligner_output,
                                                                   kmer_reference,
                                                                   current_unique_threshold,
                                                                   **kwargs)

            if len(current_read_mapping.specific_kmers) == 0:
                map_read(read, UNMAPPED_READ, kmer_reference)
            else:
                _map_read_using_specific_kmers(current_read_mapping,
                                               read, kmer_reference,
                                               current_unique_threshold)
                if read.read_status == UNIQUE_READ and ambiguous_threshold >= 0:
                    validate_uniqueness_using_unspecific(
                        current_read_mapping, read, kmer_reference,
                        ambiguous_threshold)
            aligner_output.update_reads_stats(read.identifier)
            if check_coverage:
                update_genome_coverage(read, current_read_mapping,
                                       kmer_reference, genome_list)
    except Exception as e:
        print(e)
        return None

    return aligner_output


def map_read(read, status: READ_STATUS,
             kmer_reference: KmerReference,
             genome_identifier: str = "") -> None:
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
    if genome_identifier:
        read.add_mapped_genome(genome_identifier)
    if status == UNIQUE_READ:
        if read.is_reversed:
            kmer_reference.genomes_db[
                genome_identifier].unique_reads_reverse += 1
        else:
            kmer_reference.genomes_db[genome_identifier].unique_reads += 1
    if status == AMBIGUOUS_READ:
        if read.is_reversed:
            kmer_reference.genomes_db[
                genome_identifier].ambiguous_reads_reverse += 1
        else:
            kmer_reference.genomes_db[genome_identifier].ambiguous_reads += 1


def extract_and_map_kmers_from_read(read: Read,
                                    kmer_reference: KmerReference,
                                    aligner_output: PseudoAlignerOutput,
                                    **kwargs) -> ReadKmerMapping:
    """
    This function extracts the kmers from the read and maps each kmer as specific, unspecific or neither
    It does that based on each kmer's appearance in the kmer_reference
    :param aligner_output: the pseudo aligner output object
    :param read: the current read
    :param kmer_reference: the KmerReference object
    :return: The current read's kmer mapping object
    """
    kmer_size = kmer_reference.kmer_size
    min_kmer_quality = kwargs.get("min_kmer_quality")
    max_genomes = kwargs.get("max_genomes")
    current_read_kmer_classification = ReadKmerMapping()
    for kmer_tuple in extract_kmers_from_string(read.value, kmer_size,
                                                True):
        kmer, kmer_position = kmer_tuple
        # Check if the quality of the kmer from the read is sufficient
        if min_kmer_quality is not None:
            min_kmer_quality = int(min_kmer_quality)
            mean_quality_of_kmer = read.calculate_mean_quality(kmer_position,
                                                               kmer_position + kmer_size)
            if mean_quality_of_kmer < min_kmer_quality:
                aligner_output.filtered_quality_kmers += 1
                continue
        # Check if the current kmer from the read is in the given reference
        if kmer in kmer_reference.kmer_db:
            # Check whether the kmer is mapped to more than max_genomes
            if max_genomes is not None:
                max_genomes = int(max_genomes)
                if len(kmer_reference.kmer_db[kmer]) > max_genomes:
                    aligner_output.filtered_hr_kmers += 1
                    continue
            # If only one genome has this kmer - it's a specific one
            if len(kmer_reference.kmer_db[kmer]) == 1:
                _add_specific_kmer_to_classification(kmer,
                                                     kmer_reference,
                                                     current_read_kmer_classification)
            # Unspecific Kmer
            else:
                _add_unspecific_kmer_to_classification(kmer,
                                                       kmer_reference,
                                                       current_read_kmer_classification)
    return current_read_kmer_classification


def _add_specific_kmer_to_classification(kmer: str,
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
        next(iter(kmer_reference.kmer_db[kmer]))
    if genome_identifier not in current_read_kmer_classification.specific_kmers_in_genomes:
        current_read_kmer_classification.specific_kmers_in_genomes[
            genome_identifier] = [kmer]
    else:
        current_read_kmer_classification.specific_kmers_in_genomes[
            genome_identifier].append(kmer)
    current_read_kmer_classification.specific_kmers.append(kmer)


def _add_unspecific_kmer_to_classification(kmer: str,
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
    for genome_identifier in kmer_reference.kmer_db[kmer]:
        if genome_identifier not in current_read_kmer_classification.unspecific_kmers_in_genomes:
            current_read_kmer_classification.unspecific_kmers_in_genomes[
                genome_identifier] = [kmer]
        else:
            current_read_kmer_classification.unspecific_kmers_in_genomes[
                genome_identifier].append(kmer)
    current_read_kmer_classification.unspecific_kmers.append(kmer)


def _map_read_using_specific_kmers(current_read_mapping: ReadKmerMapping,
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
            next(iter(current_read_mapping.specific_kmers_in_genomes))
        map_read(read, UNIQUE_READ, kmer_reference, genome_identifier)
    else:
        # By this point - there are at least two genomes that have specific kmers from this read
        # hence, frequent, second_frequent won't return None
        frequent, second_frequent = find_two_most_frequent_genomes(
            current_read_mapping.specific_kmers_in_genomes)
        frequency_difference = frequent[1] - second_frequent[1]
        if frequency_difference >= unique_threshold:
            genome_identifier = frequent[0]
            map_read(read, UNIQUE_READ, kmer_reference,
                     genome_identifier)
        else:
            for genome_identifier in current_read_mapping.specific_kmers_in_genomes:
                map_read(read, AMBIGUOUS_READ, kmer_reference,
                         genome_identifier)


def find_two_most_frequent_genomes(specific_kmers_in_genomes: Dict[
    str, List[str]]) -> Tuple[Tuple[str, int], Tuple[str, int]]:
    """
    This function receives a dictionary of genomes and the specific kmers that matched them.
    It scans through the dictionary and returns the two highest genomes with most kmers
    Note: this function might return None if there are less than 2 items in specific_kmers_in_genomes
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
    # mypy might be unhappy - but this function will only be called in the case when first,second won't be None
    return first, second


def validate_uniqueness_using_unspecific(current_read_mapping: ReadKmerMapping,
                                         read: Read,
                                         kmer_reference: KmerReference,
                                         ambiguous_threshold: int) -> None:
    """
    This is a validation function that a uniquely mapped read matches up with the unspecific kmers
    It calculates the total kmer-count for all genomes, and records the highest total count
    The difference between the max and the count of the uniquely mapped genome count must be below or equal to p
    If not, classify this read as ambiguous
    :param current_read_mapping: the mapping for the read - specific and unspecific kmers
    :param read: the current read
    :param kmer_reference: the kmer reference to access the genomes
    :param ambiguous_threshold: 'p' the ambiguous threshold
    :return: None
    """
    specific_kmers_dict = current_read_mapping.specific_kmers_in_genomes
    unspecific_kmers_dict = current_read_mapping.unspecific_kmers_in_genomes

    total_counter_dict = {key: len(specific_kmers_dict.get(key)) +
                               len(unspecific_kmers_dict.get(key, [])) for key
                          in
                          set(specific_kmers_dict)}

    max_count = max(total_counter_dict.values())
    # Our read was uniquely mapped - hence the mapped_genomes to it is only one
    unique_genome_identifier = next(iter(read.mapped_genomes))
    map_count = total_counter_dict[unique_genome_identifier]
    if max_count - map_count > ambiguous_threshold:
        # This read is no longer uniquely mapped to the genome
        if read.is_reversed:
            kmer_reference.genomes_db[
                unique_genome_identifier].unique_reads_reverse -= 1
        else:
            kmer_reference.genomes_db[
                unique_genome_identifier].unique_reads -= 1
        # ambiguously map all genomes
        for genome_identifier, total_count in total_counter_dict.items():
            if total_count >= map_count:
                map_read(read, AMBIGUOUS_READ, kmer_reference,
                         genome_identifier)
