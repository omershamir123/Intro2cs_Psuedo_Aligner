# FILE : test_validators.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: This module contains all tests of the validators module
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:

import pytest

from genome import ReferencedGenome
from program_constants import KDB_FILE_TYPES, FASTQ_FILE_TYPES
from validators import validate_not_empty, validate_file_type, \
    validate_above_value, validate_values_in_given_list


def test_validate_not_empty():
    """
    This function tests that various params aren't None
    :return:
    """
    assert validate_not_empty(None) == False
    assert validate_not_empty(3) == True
    assert validate_not_empty(ReferencedGenome("gdfgdf", "A", 3)) == True


def test_validate_file_type():
    """
    Thisfunction tests given file paths and their matching file types and makes sure
    it's the right file type
    :return:
    """
    try:
        validate_file_type("example.kdb", KDB_FILE_TYPES)
        validate_file_type("example.fa.kdb", KDB_FILE_TYPES)
    except TypeError:
        assert False, "An error has occurred while checking the validity of VALID kdb file names"

    with pytest.raises(TypeError):
        validate_file_type("example.fa", KDB_FILE_TYPES)
        validate_file_type("example.kdb", FASTQ_FILE_TYPES)
        validate_file_type("example.fa", KDB_FILE_TYPES + FASTQ_FILE_TYPES)


def test_validate_above_value():
    """
    This function tests the validation of an object being above another value
    :return:
    """
    try:
        validate_above_value(5,4,False)
        validate_above_value(5,5,True)
        validate_above_value(5.0000000,5,True)
        validate_above_value(3.87867, 3.8 ,False)
    except ValueError as e:
        # all above lines are valid, hence no exception should be thrown
        assert False, "validate above value has failed with true inputs - {}".format(e)

    with pytest.raises(ValueError):
        validate_above_value(5,6,True)
        validate_above_value(5,5.0,False)


def test_values_in_given_list():
    """
    This function checks the validation that all items in a given iterable are in the given list
    :return:
    """
    try:
        validate_values_in_given_list("ABCDEFG",["A", "C", "B", "G", "F", "G", "D", "E", "N"])
        validate_values_in_given_list("AAAAAAAAAAAAAAA", ["A"])
        validate_values_in_given_list("", ["A"])
    except ValueError as e:
        assert False, "validate_values_int_list should not have failed for the input: {}".format(e)

    with pytest.raises(ValueError):
        validate_values_in_given_list("ABCD", [])
        validate_values_in_given_list("ABCD", ["a", "A", "N","C", "D"])

