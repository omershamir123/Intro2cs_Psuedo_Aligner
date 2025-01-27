# FILE : pseudo_aligner.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: This file includes the pseudo aligner class that handles the pseudo_aligner algorithm
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
import file_handlers


class PseudoAligner:
    def __init__(self):
        self._reads = {}
        self._unique_mapped_reads = 0
        self._ambiguous_mapped_reads = 0
        self._unmapped_reads = 0
        self._filtered_quality_reads = 0
        self._filtered_quality_kmers = 0
        self._filtered_hr_kmers = 0


    def load_all_reads(self, fastq_file:str) -> None:
        for read in file_handlers.parse_fastq_file(fastq_file):
            if read.identifier in self._reads:
                raise ValueError("Duplicate reads detected")
            self._reads[read.identifier] = read

