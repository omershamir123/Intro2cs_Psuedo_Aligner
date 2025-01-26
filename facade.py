# FILE : facade.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : Final Project 2025
# DESCRIPTION: the facade of the pseudo aligner algorithm, gets the input from command args
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:

from argparse import Namespace
from typing import Dict, Callable

import validators


def build_reference(args: Namespace) -> bool:
    try:
        validators.validate_above_value(args.kmersize, 0, allow_equality=False)
    except ValueError:
        return False




def initialize_function_calls() -> Dict[str, Callable[[Namespace], bool]]:
    """
    This function initializes the function dictionary for each task in the program
    :return: Dictionary of task names to functions
    """
    function_calls: Dict[str, Callable[[Namespace], bool]] = {}
    function_calls["reference"] = build_reference
    return function_calls

def start_program(args: Namespace) -> None:
    """
    This function starts the program with the given parameters from the user
    :param args: the arguments given by the user
    :return: None
    """
    function_calls = initialize_function_calls()
    if args.command not in function_calls:
        return
    function_calls[args.command](args)
