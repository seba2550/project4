# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    expected_align_mat = np.array([[  0, -np.inf, -np.inf, -np.inf],
 [-np.inf,   5, -10, -11,],
 [-np.inf, -11,   4,  -8,],
 [-np.inf, -10,  -1,   5,],
 [-np.inf, -11,  -6,   4,]])

    expected_gapA_mat = np.array([[-np.inf, -10 , -10, -10],
 [-np.inf,   -21, -6, -7,],
 [-np.inf, -21,   -17,  -7,],
 [-np.inf, -21,  -18,   -12,],
 [-np.inf, -21,  -19,   -17,]])

    expected_gapB_mat = np.array([[-np.inf, -np.inf , -np.inf, -np.inf],
 [-10,   -21, -21, -21,],
 [-10, -6,   -17,  -18,],
 [-10, -7,  -7,   -18,],
 [-10, -8,  -8,   -6,]])


    

    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    sub_mat = "./substitution_matrices/BLOSUM62.mat"
    nw = NeedlemanWunsch(sub_mat, -10, -1)
    score, align_seq1, align_seq2 = nw.align(seq1, seq2)
    assert np.array_equal(expected_align_mat, nw._align_matrix)
    assert np.array_equal(expected_gapA_mat, nw._gapA_matrix)
    assert np.array_equal(expected_gapB_mat, nw._gapB_matrix)

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    sub_mat = "./substitution_matrices/BLOSUM62.mat"

    expected_score = 17 # Given by the README.md file in the Github Repo for this assignment
    expected_seq3_align = "MAVHQLIRRP"
    expected_seq4_align = "M---QLIRHP"
    
    nw = NeedlemanWunsch(sub_mat, -10, -1)

    nw_align_res = nw.align(seq3, seq4)

    assert expected_score == nw_align_res[0]
    assert expected_seq3_align == nw_align_res[1]
    assert expected_seq4_align == nw_align_res[2]






