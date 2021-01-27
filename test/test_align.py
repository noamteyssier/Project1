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

	sr = algs.SeqReader()
	header, seq = sr.read_fasta("../data/prot-0004.fa")

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


	pw = algs.PairwiseAligner(
		"../data/prot-0004.fa",
		"../data/prot-0004.fa",
		"../scoring_matrices/BLOSUM50.mat"
		)

	assert (pw.submat - expected_mat).sum() == 0

def test_identical():
	expected_seq = "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"

	for method in [algs.SmithWaterman, algs.NeedlemanWunsch]:
		m_obj = method(
			"../data/prot-0004.fa",
			"../data/prot-0004.fa",
			"../scoring_matrices/BLOSUM50.mat"
			)

		score, q1, q2 = m_obj.align()

		assert q1 == q2
		assert expected_seq == q1

def test_alignment_score():
	sw = algs.NeedlemanWunsch(
		"../data/prot-0004.fa",
		"../data/prot-0008.fa",
		"../scoring_matrices/BLOSUM50.mat",
		gamma = 5, opening_penalty = 10
	)
	score, q1, q2 = sw.align()
	print(q1)
	print(q2)
	print(score)

	assert True

# runs tests
def main():

	test_fasta_io()

	test_scoring_matrix_io()

	test_identical()

	test_alignment_score()

if __name__ == "__main__":
	main()
