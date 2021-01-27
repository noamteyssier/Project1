#!/usr/bin/env python3

import pytest
import numpy as np
import sys
import os

sys.path.append(
	os.path.abspath(
		os.path.join(
			os.path.dirname(__file__),
			os.path.pardir
		)
	)
)

from align import algs

@pytest.fixture
def some_relevant_data():
	return np.ones(10)

def test_fasta_io():

	expected_header = "d1flp__ 1.1.1.1.2"
	expected_seq = "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"

	io = algs.io()
	header, seq = io.read_fasta("../sequences/prot-0004.fa")

	assert (header == expected_header) & (seq == expected_seq)

def test_scoring_matrix_io():

	# BLOSUM50 Hardcorded
	expected_mat = np.array([
		[5,-2,-1,-2,-1,-1,-1,0,-2,-1,-2,-1,-1,-3,-1,1,0,-3,-2,0,-2,-1,-1,-5],
		[-2,7,-1,-2,-4,1,0,-3,0,-4,-3,3,-2,-3,-3,-1,-1,-3,-1,-3,-1,0,-1,-5],
		[-1,-1,7,2,-2,0,0,0,1,-3,-4,0,-2,-4,-2,1,0,-4,-2,-3,4,0,-1,-5],
		[-2,-2,2,8,-4,0,2,-1,-1,-4,-4,-1,-4,-5,-1,0,-1,-5,-3,-4,5,1,-1,-5],
		[-1,-4,-2,-4,13,-3,-3,-3,-3,-2,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1,-3,-3,-2,-5],
		[-1,1,0,0,-3,7,2,-2,1,-3,-2,2,0,-4,-1,0,-1,-1,-1,-3,0,4,-1,-5],
		[-1,0,0,2,-3,2,6,-3,0,-4,-3,1,-2,-3,-1,-1,-1,-3,-2,-3,1,5,-1,-5],
		[0,-3,0,-1,-3,-2,-3,8,-2,-4,-4,-2,-3,-4,-2,0,-2,-3,-3,-4,-1,-2,-2,-5],
		[-2,0,1,-1,-3,1,0,-2,10,-4,-3,0,-1,-1,-2,-1,-2,-3,2,-4,0,0,-1,-5],
		[-1,-4,-3,-4,-2,-3,-4,-4,-4,5,2,-3,2,0,-3,-3,-1,-3,-1,4,-4,-3,-1,-5],
		[-2,-3,-4,-4,-2,-2,-3,-4,-3,2,5,-3,3,1,-4,-3,-1,-2,-1,1,-4,-3,-1,-5],
		[-1,3,0,-1,-3,2,1,-2,0,-3,-3,6,-2,-4,-1,0,-1,-3,-2,-3,0,1,-1,-5],
		[-1,-2,-2,-4,-2,0,-2,-3,-1,2,3,-2,7,0,-3,-2,-1,-1,0,1,-3,-1,-1,-5],
		[-3,-3,-4,-5,-2,-4,-3,-4,-1,0,1,-4,0,8,-4,-3,-2,1,4,-1,-4,-4,-2,-5],
		[-1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,-1,-1,-4,-3,-3,-2,-1,-2,-5],
		[1,-1,1,0,-1,0,-1,0,-1,-3,-3,0,-2,-3,-1,5,2,-4,-2,-2,0,0,-1,-5],
		[0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,2,5,-3,-2,0,0,-1,0,-5],
		[-3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1,1,-4,-4,-3,15,2,-3,-5,-2,-3,-5],
		[-2,-1,-2,-3,-3,-1,-2,-3,2,-1,-1,-2,0,4,-3,-2,-2,2,8,-1,-3,-2,-1,-5],
		[0,-3,-3,-4,-1,-3,-3,-4,-4,4,1,-3,1,-1,-3,-2,0,-3,-1,5,-4,-3,-1,-5],
		[-2,-1,4,5,-3,0,1,-1,0,-4,-4,0,-3,-4,-2,0,0,-5,-3,-4,5,2,-1,-5],
		[-1,0,0,1,-3,4,5,-2,0,-3,-3,1,-1,-4,-1,0,-1,-2,-2,-3,2,5,-1,-5],
		[-1,-1,-1,-1,-2,-1,-1,-2,-1,-1,-1,-1,-1,-2,-2,-1,0,-3,-1,-1,-1,-1,-1,-5],
		[-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,1]
	])

	io = algs.io()
	idx_lookup, submat = io.read_substitution_matrix("../scoring_matrices/BLOSUM50.mat")

	assert (submat - expected_mat).sum() == 0

def test_identical():
	expected_seq = "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"

	io = algs.io()
	header, seq = io.read_fasta("../sequences/prot-0004.fa")
	idx_lookup, submat = io.read_substitution_matrix("../scoring_matrices/BLOSUM50.mat")

	for method in [algs.SmithWaterman, algs.NeedlemanWunsch]:
		m_obj = method(
			seq, seq, submat, idx_lookup
			)

		score, q1, q2 = m_obj.align()

		assert q1 == q2
		assert expected_seq == q1

def test_alignment_score():

	io = algs.io()
	h1, s1 = io.read_fasta("../sequences/prot-0004.fa")
	h2, s2 = io.read_fasta("../sequences/prot-0008.fa")
	idx_lookup, submat = io.read_substitution_matrix("../scoring_matrices/BLOSUM50.mat")

	sw = algs.NeedlemanWunsch(
		s1, s2, submat, idx_lookup,
		opening_penalty = 10, extension_penalty=5
	)

	score, q1, q2 = sw.align()
	assert (score + 17 == 0)

# runs tests
def main():

	test_fasta_io()

	test_scoring_matrix_io()

	test_identical()

	test_alignment_score()

if __name__ == "__main__":
	main()
