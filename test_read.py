# FILE : test_read.py
# WRITER : Omer Shamir , omer.shamir , 322593120
# EXERCISE : intro2cs Final Project 2025
# DESCRIPTION: This module tests the entire read module functionality
# STUDENTS I DISCUSSED THE EXERCISE WITH:
# WEB PAGES I USED:
# NOTES:
import pytest

import read
from read import reverse_complement


def test_read_initialization():
    with pytest.raises(ValueError):
        read.Read("test_read", "ATCGN", "")
        read.Read("test_read", "ATCGN@", "")
        read.Read("test_read", "ATCGN", "@@@@")
        read.Read("test_read", "ATCGN", "@@@@"+chr(32))

    read1 = read.Read("test_read", "ATCGN", "@@@@A")
    assert read1.quality == [31,31,31,31,32]
    assert read1.reverse_complement == "NCGAT"
    assert read1.reversed_quality == [32,31,31,31,31]


def test_reverse_complement():
    assert reverse_complement("ATCGT") == "ACGAT"
    assert reverse_complement("ATCGN") == "NCGAT"
    assert reverse_complement("GG") == "CC"

def test_read_mean_quality():
    read1 = read.Read("test_read", "ATCGN", "@@@@A")
    assert read1.quality == [31,31,31,31,32]
    assert read1.calculate_mean_quality() == 31.2
    assert read1.calculate_mean_quality(1,3) == float(31)
    assert read1.calculate_mean_quality(2) == float(31.33333333333333333333)
    assert read1.calculate_mean_quality(4) == float(32)
    assert read1.calculate_mean_quality(ending_index=4) == float(31)