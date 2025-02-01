# FILE : facade.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : Final Project 2025
# DESCRIPTION: the facade of the pseudo aligner algorithm, gets the input from command args
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
from argparse import Namespace
from typing import Dict, Callable, Optional

import file_handlers
import pseudo_aligner
import validators
from kmer_reference import KmerReference
from program_constants import KDB_FILE_TYPES, ALN_FILE_TYPES
from pseudo_aligner import AlnFileDataObject

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
    if not args.genomefile or not args.kmer_size or not args.referencefile:
        print(
            "For the reference command - a genome file, a reference file and a kmer-size must be provided")
        return

    reference = build_reference(args.kmer_size, args.genomefile)
    if reference is not None:
        # By now we have built the entire reference from all genomes.
        # Now, if filer-similar flag is True we check the similarity score between each genome
        # and update the reference file accordingly
        if args.filter_similar:
            print(reference.filter_genomes_logic(args.similarity_threshold))
        file_handlers.write_to_pickle_file(args.referencefile, reference,
                                           KDB_FILE_TYPES)


def extract_reference(args: Namespace) -> Optional[KmerReference]:
    """
    This function checks the arguements and extracts the reference file either by building it
    or by running the build reference function.
    If genome file, a kmer-reference file and a kmer-size are provided, it's an invalid input
    :param args: the arguments given by the user
    :return: if the parameters are valid - a kmer-reference object else None
    """
    if args.genomefile and args.kmer_size:
        if args.referencefile:
            print(
                "Both a genome file and a reference file were provided - invalid input for the {} command".format(
                    args.task))
            return None
        reference = build_reference(args.kmer_size, args.genomefile)
    else:
        if not args.referencefile:
            print(
                "No reference file given - invalid input for the {} command".format(
                    args.task))
            return None
        # mypy might be unhappy about this but there is an extra check later on that
        # reference is indeed a KmerReference object
        reference = file_handlers.decompress_pickle_file(args.referencefile,
                                                         KDB_FILE_TYPES)
    if isinstance(reference, KmerReference):
        return reference
    return None


def dumpref_command(args: Namespace) -> None:
    """
    This function carries out the dumpref command
    :param args: the arguments given by the user
    :return: None
    """
    reference = extract_reference(args)
    if reference is not None:
        if args.filter_similar:
            print(reference.filter_genomes_logic(args.similarity_threshold))
        print(reference.to_json())


def align_command(args: Namespace) -> None:
    """
    This function handles the align command - runs the pseudo align algorithm and saves the output to an aln file
    The function also runs validity checks that the user has entered a valid combination of arguments
    For example the user must enter a value for the reads file and the align file
    :param args: arguments given by the user
    :return: None
    """
    if not args.alignfile or not args.reads:
        print(
            "For the alignment command - an align file path and a reads file must be provided")
        return
    reference = extract_reference(args)
    if reference is not None:
        try:
            validators.validate_above_value(args.unique_threshold, 0, True)
        except ValueError:
            print("Invalid unique threshold value provided")
            return
        align_output = pseudo_aligner.align_algorithm(args.reads, reference,
                                                      args.unique_threshold,
                                                      args.ambiguous_threhold, **vars(args))
        if align_output is not None:
            align_file_object = align_output.convert_to_aln_object(args.reverse_complement)
            if args.coverage:
                try:
                    validators.validate_above_value(args.min_coverage, 0, True)
                except ValueError:
                    print("Invalid min coverage provided")
                    return
                print(align_file_object.to_coverage_json(args.full_coverage,
                                                         args.min_coverage))
            file_handlers.write_to_pickle_file(args.alignfile, align_file_object,
                                               ALN_FILE_TYPES)


def dumpalign_command(args: Namespace) -> None:
    align_file_object = None
    if args.alignfile:
        if any(argument is not None for argument in
               [args.reads, args.referencefile, args.genomefile,
                args.kmer_size]):
            print(
                "For the dumpalign command - if -a is a given flag, no other flags can be inputted")
            return
        align_file_object = file_handlers.decompress_pickle_file(args.alignfile,
                                                                 ALN_FILE_TYPES)
    else:
        if not args.reads:
            print(
                "No reads file given - invalid input for the dumpalign command")
            return
        reference = extract_reference(args)
        if reference is not None:
            try:
                validators.validate_above_value(args.unique_threshold, 0, True)
            except ValueError:
                print("Invalid unique threshold value provided")
                return
            align_output = pseudo_aligner.align_algorithm(args.reads, reference,
                                                          args.unique_threshold,
                                                          args.ambiguous_threhold, **vars(args))
            if align_output is not None:
                align_file_object = align_output.convert_to_aln_object(args.reverse_complement)
            else:
                align_file_object = None
    if align_file_object is not None and isinstance(align_file_object, AlnFileDataObject):
        if args.coverage:
            try:
                validators.validate_above_value(args.min_coverage, 0, True)
            except ValueError:
                print("Invalid min coverage provided")
                return

            print(align_file_object.to_coverage_json(args.full_coverage, args.min_coverage))
        print(align_file_object.to_json())


def initialize_function_calls() -> Dict[str, Callable[[Namespace], None]]:
    """
    This function initializes the function dictionary for each task in the program
    :return: Dictionary of task names to functions
    """
    function_calls: Dict[str, Callable[[Namespace], None]] = {
        "reference": reference_command, "dumpref": dumpref_command,
        "align": align_command, "dumpalign": dumpalign_command}
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
