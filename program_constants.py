# FILE : program_constants.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: This module simply contains all the constants used in this program
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:

FASTA_FILE_TYPES = [".fa"]
FASTQ_FILE_TYPES = [".fq"]
KDB_FILE_TYPES = [".kdb"]
ALN_FILE_TYPES = [".aln"]

ALLOWED_DNA_VALUES = ["A", "T", "C", "G", "N"]
WILDCARD_READINGS = ["N"]

UNIQUE_KMER = 0
MULTI_MAP_KMER = 1
KMER_TYPE = int

UNIQUE_READ = 0
AMBIGUOUS_READ = 1
UNMAPPED_READ = 2
READ_STATUS = int

COMPLEMENTS_DICT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
