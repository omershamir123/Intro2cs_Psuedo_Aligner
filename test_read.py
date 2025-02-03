# FILE : test_read.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: This module tests the entire read module functionality
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
import numpy as np
import pytest

import read
from read import reverse_complement


def test_read_initialization():
    """
    This function checks the initialization of Read, whether it's because of wrong values, no quality, invalid quality
    :return:
    """
    with pytest.raises(ValueError):
        # no quality
        read.Read("test_read", "ATCGN", "")
        # wrong letter in sequence
        read.Read("test_read", "ATCGN@", "")
        # wrong quality length
        read.Read("test_read", "ATCGN", "@@@@")
        # invalid char in quality (below 33)
        read.Read("test_read", "ATCGN", "@@@@"+chr(32))

    read1 = read.Read("test_read", "ATCGN", "@@@@A")
    assert np.array_equal(read1.quality, np.array([31,31,31,31,32]))
    assert read1.reverse_complement == "NCGAT"
    assert np.array_equal(read1.reversed_quality, np.array([32,31,31,31,31]))


def test_reverse_complement():
    """
    This function checks the calculation of the reverse complement
    :return:
    """
    assert reverse_complement("ATCGT") == "ACGAT"
    assert reverse_complement("ATCGN") == "NCGAT"
    assert reverse_complement("GG") == "CC"

def test_read_mean_quality():
    """
    This function checks the calculation of the mean quality (uses np.mean())
    :return:
    """
    read1 = read.Read("test_read", "ATCGN", "@@@@A")
    assert np.array_equal(read1.quality, np.array([31,31,31,31,32]))
    assert read1.calculate_mean_quality() == 31.2
    assert read1.calculate_mean_quality(1,3) == float(31)
    assert read1.calculate_mean_quality(2) == float(31.33333333333333333333)
    assert read1.calculate_mean_quality(4) == float(32)
    assert read1.calculate_mean_quality(ending_index=4) == float(31)