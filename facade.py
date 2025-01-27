# FILE : facade.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : Final Project 2025
# DESCRIPTION: the facade of the pseudo aligner algorithm, gets the input from command args
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
import json
from argparse import Namespace
from typing import Dict, Callable, Optional

import file_handlers
import validators
from kmer_reference import KmerReference

from validators import validate_not_empty


def build_reference(kmer_size: int, genomefile: str) -> Optional[KmerReference]:
    """
    This function receives the kmer_size and the genome file path and builds a kmer-reference object
    If one of the steps in the build phase isn't successful, the returned object is None
    :param genomefile: the path to the FASTA file
    :param kmer_size: the size of each Kmer in our reference
    :return: the KmerReference object if the build was successful, otherwise None
    """
    try:
        validators.validate_above_value(kmer_size, 0, allow_equality=False)
    except (ValueError, TypeError) as e:
        print(e)
        return None
    reference = KmerReference(kmer_size)
    build_successful = reference.build_kmer_reference(genomefile)
    # If build was not successful for various reasons, the function will not continue on
    if build_successful:
        return reference
    return None


def reference_command(args: Namespace) -> None:
    """
    This function carries out the reference command - of building a reference and writing it to a kdb file
    :param args: the arguments given by the user
    :return: None
    """
    reference = build_reference(args.kmer_size, args.genomefile)
    if reference is not None:
        if validate_not_empty(args.referencefile):
            # TODO Think about the case when the write was unsuccessful - is there a need to return something?
            file_handlers.write_to_kdb_file(args.referencefile, reference)
        else:
            print("There was a problem with the given kdb file")

def dumpref_command(args: Namespace) -> None:
    """
    This function carries out the dumpref command
    :param args:
    :return:
    """
    # TODO make sure there is a definitive was to know the parameters (there are no extra for example)
    if validate_not_empty(args.genomefile) and validate_not_empty(args.kmer_size):
        reference = build_reference(args.kmer_size, args.genomefile)
    else:
        reference = file_handlers.decompress_kdb_file(args.referencefile)
    if reference is not None:
        print(reference.to_json())


def align_command(args: Namespace) -> None:
    # TODO make sure there is a definitive was to know the parameters (there are no extra for example)
    if validate_not_empty(args.genomefile) and validate_not_empty(args.kmer_size):
        reference = build_reference(args.kmer_size, args.genomefile)
    else:
        reference = file_handlers.decompress_kdb_file(args.referencefile)
    if reference is not None:
        pseudo


# TODO Callable return type
def initialize_function_calls() -> Dict[str, Callable[[Namespace], None]]:
    """
    This function initializes the function dictionary for each task in the program
    :return: Dictionary of task names to functions
    """
    function_calls: Dict[str, Callable[[Namespace], None]] = {}
    function_calls["reference"] = reference_command
    function_calls["dumpref"] = dumpref_command
    return function_calls


def start_program(args: Namespace) -> None:
    """
    This function starts the program with the given parameters from the user
    :param args: the arguments given by the user
    :return: None
    """
    function_calls = initialize_function_calls()
    if args.task not in function_calls:
        print("Invalid task")
        return
    function_calls[args.task](args)
