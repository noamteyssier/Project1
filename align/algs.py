#!/usr/bin/env python3

"""
Align
====================================
Tools for Sequence / Substitution IO and Sequence Alignment (NeedlemanWunsch + SmithWaterman)
"""

import numpy as np
import argparse
import os

class io:

	"""
	Functions for reading in sequences from Fasta files and substitution matrices
	"""

	def _tokenize(self, s):

		"""
		Tokenize white space in a string and return generator

		Parameters
		----------
		s
			string to split on white space and yield tokens

		"""

		for i in s.strip().split(" "):
			if i == "":
				continue
			yield i

	def read_fasta(self, file_path):
		"""
		Reads in a fasta file and returns the first sequence as (header, sequence) tuple

		Parameters
		---------
		file_path

			path to a \*.fa format file.

		Returns
		-------
		tuple
			(header, sequence)

		"""
		assert file_path[-3:] == ".fa"

		header = ""
		sequence = ""

		with open(file_path, 'r') as f:

			while True:

				try:
					# read in line
					line = next(f)

					# check if header
					if ">" == line[0]:

						if len(header) == 0:
							header = line[1:].strip()

						# second header found, breaking loop
						else:
							break

					# append sequence if multi-line
					else:
						sequence += line.strip().upper()

				# EOF
				except StopIteration:
					break

		# return header and sequence as tuple
		return (header, sequence)

	def read_substitution_matrix(self, file_path):

		"""

		Reads in substitution matrix and populates idx_lookup and submat

		Parameters
		---------
		file_path

			path to a \*.mat formatted matrix.

		Returns
		-------
		tuple
			(character lookup, score matrix)

		"""

		assert file_path[-4:] == ".mat"

		idx_lookup = {}

		row_idx = 0
		with open(file_path, "r") as f:
			while True:
				try:
					line = next(f)

					# ignores comments
					if "#" == line[0]:
						continue
					elif "." == line[0]:
						continue

					# first uncommented line expected as header
					elif len(idx_lookup) == 0:
						idx_lookup = {
							char : index for index,char in enumerate(self._tokenize(line))
							}

						# instantiate substitution matrix
						submat = np.zeros(
							(len(idx_lookup), len(idx_lookup)),
							dtype=int
							)

					# otherwise fill substitution matrix
					else:
						for col_idx, val in enumerate(self._tokenize(line)):
							submat[row_idx][col_idx] = int(val)

						row_idx += 1

				# EOF
				except StopIteration:
					break


		return (idx_lookup, submat)


class PairwiseAligner:

	"""
	Parent class for Alignment Methods

	Parameters
	----------

	seq1
		string of reference protein sequence

	seq2
	 	string of query protein sequence

	submat
		a subsitution matrix as a 2D numpy array

	idx_lookup
		a dictionary to map protein characters to their indices in submat

	opening_penalty
		the cost of opening a gap

	extension_penalty
		the cost of extending a previously opened gap

	"""

	def __init__(self, seq1, seq2, submat, idx_lookup, opening_penalty = 11, extension_penalty = 3):

		# sequences
		self.seq1 = seq1
		self.seq2 = seq2

		# substitution matrix
		self.submat = submat

		# character lookup table
		self.idx_lookup = idx_lookup

		# gap penalties
		self.opening_penalty = opening_penalty
		self.extension_penalty = extension_penalty

		# expected sequence sizes
		self.n = len(self.seq1) + 1
		self.m = len(self.seq2) + 1

		# scoring matrix
		self.score_mat = self._InitMatrix()

		# trace matrix
		## -1 : stop
		##  0 : diag
		##  1 : left
		##  2 : up
		self.trace_mat = self._InitMatrix(dtype=int)

		# gap matrices for affine gap calculations
		self.gap_l = self._InitMatrix(dtype=int)
		self.gap_u = self._InitMatrix(dtype=int)

	def _InitMatrix(self, dtype=float):

		"""
		Initialize matrices of an expected size

		Returns
		-------
		2darray
			An empty matrix of size (mxn) where m and n are the length of the given sequences
		"""

		shape = (self.n, self.m)
		return np.zeros(shape, dtype=dtype)

	def _GapLoss(self, idx, jdx, movement):
		"""
		Implements Gap Loss

		G = opening_penalty + (num_gaps * extension_penalty)

		Parameters
		---------
		idx
			row position
		jdx
			col position
		movement
			the prefix direction (1 : left prefix, 2 : up prefix)

		Returns
		-------
		numeric
			returns the gap loss at a given position with a given movement
		"""

		# no opening penalty required (non-affine gap)
		if self.opening_penalty == 0:
			return self.extension_penalty

		# opening penalty required (affine gap)
		else:

			# left gap state matrix
			if movement == 1:

				# case where gap has not been already opened
				if self.gap_l[(idx, jdx-1)] == 0:
					return self.opening_penalty #+ self.extension_penalty

				# gap has already been opened
				else:
					return self.extension_penalty


			# upward gap state matrix
			else:

				# case where gap has not been already opened
				if self.gap_u[(idx-1,jdx)] == 0:

					return self.opening_penalty + self.extension_penalty

				# gap has already been opened
				else:
					return self.extension_penalty

	def _SubScore(self, idx, jdx):

		"""
		Returns a substitution score from sub matrix for two given indices

		Because of 1gap padded strings idx and jdx will be treated as 1-indexed

		Parameters
		---------
		idx
			row position
		jdx
			col position

		Returns
		-------
		numeric
			returns the substitution cost at a given position
		"""

		c1 = self.seq1[idx-1]
		c2 = self.seq2[jdx-1]

		if c1 not in self.idx_lookup:
			c1 = "*"
		if c2 not in self.idx_lookup:
			c2 = "*"

		sub_idx = self.idx_lookup[self.seq1[idx-1]]
		sub_jdx = self.idx_lookup[self.seq2[jdx-1]]
		return self.submat[(sub_idx, sub_jdx)]

	def _ScoreDiag(self, idx, jdx):

		"""
		Calculates score with upper-left diagonal prefix

		Parameters
		---------
		idx
			row position
		jdx
			col position

		Returns
		-------
		numeric
			returns the score of diagonal-prefix movement
		"""

		return self.score_mat[(idx-1, jdx-1)] + self._SubScore(idx, jdx)

	def _ScoreLeft(self, idx, jdx):

		"""
		Calculates score with left prefix

		Parameters
		---------
		idx
			row position
		jdx
			col position

		Returns
		-------
		numeric
			returns the score of left-prefix movement
		"""

		return self.score_mat[(idx, jdx-1)] - self._GapLoss(idx, jdx, movement=1)

	def _ScoreUp(self, idx, jdx):

		"""
		Calculates score with upper prefix

		Parameters
		---------
		idx
			row position
		jdx
			col position

		Returns
		-------
		numeric
			returns the score of up-prefix movement
		"""

		return self.score_mat[(idx-1, jdx)] - self._GapLoss(idx, jdx, movement=2)

	def _UpdateMatrices(self, idx, jdx, chosen_score, chosen_direction):

		"""
		Updates Score/Trace/Gap Matrices with chosen score/direction

		Parameters
		---------
		idx
			row position
		jdx
			col position
		chosen_score
			the maximum score of [diagonal, left, up]
		chosen_direction
			the direction of the maximum score
		"""

		self.score_mat[(idx, jdx)] = chosen_score
		self.trace_mat[(idx, jdx)] = chosen_direction
		self._UpdateGaps(idx, jdx, chosen_direction)

	def _UpdateGaps(self, idx, jdx, chosen_direction):

		"""
		Updates required gap matrix with opening/extension decision

		Parameters
		---------
		idx
			row position
		jdx
			col position
		chosen_direction
			the direction of the maximum score
		"""

		if chosen_direction == 1:
			self.gap_l[(idx, jdx)] = 1

		elif chosen_direction == 2:
			self.gap_u[(idx, jdx)] = 1

	def _ChooseDirection(self, score_list):

		"""
		Returns the direction of the highest score.

		If there are ties with the diagonal the diagonal will win.
		Otherwise the choice is random between the two prefixes

		Parameters
		---------
		score_list
			1darray of of scores from [diag, left, up]

		Returns
		-------
		int
			direction to prefix
		"""

		m = np.max(score_list)
		mask = (score_list >= m)

		# multiple possible choices
		if mask.sum() > 1:

			# if diagonal is possible follow diagonal
			if mask[0]:
				return 0

			# else if both directions are possible return a random choice
			elif (mask[1] == mask[2]) & mask[1]:
				return np.random.choice([1,2])


		# return maximum index
		return np.where(mask)[0]

	def align(self):

		"""
		Align 2 protein sequences

		Returns
		-------
		tuple
			(score, aligned_seq1, aligned_seq2)
		"""

		return self.__align__()


class SmithWaterman(PairwiseAligner):

	"""
	Protein Sequencing Alignment using Smith-Waterman Local Alignment

	Parameters
	----------

	seq1
		string of reference protein sequence

	seq2
	 	string of query protein sequence

	submat
		a subsitution matrix as a 2D numpy array

	idx_lookup
		a dictionary to map protein characters to their indices in submat

	opening_penalty
		the cost of opening a gap

	extension_penalty
		the cost of extending a previously opened gap
	"""

	def __align__(self):

		"""
		Implementation of alignment for SmithWaterman
		"""

		# row iter
		for idx in np.arange(self.score_mat.shape[0]):

			# col iter
			for jdx in np.arange(self.score_mat.shape[1]):

				# top left corner (double gap)
				if (idx == 0) and (jdx == 0):
					self.score_mat[(idx, jdx)] = 0
					self.trace_mat[(idx, jdx)] = -1

				# top row
				elif (idx == 0) and (jdx > 0):
					score_left = self._ScoreLeft(idx, jdx)

					if score_left < 0:
						score_left = 0

					self.score_mat[(idx,jdx)] = score_left
					self.trace_mat[(idx,jdx)] = 1

				# left column
				elif (idx > 0) and (jdx == 0):
					score_up = self._ScoreUp(idx, jdx)

					if score_up < 0:
						score_up = 0

					self.score_mat[(idx, jdx)] = score_up
					self.trace_mat[(idx,jdx)] = 2

				# main loop
				else:
					score_list = np.array([
						self._ScoreDiag(idx, jdx),
						self._ScoreLeft(idx, jdx),
						self._ScoreUp(idx, jdx)
						])

					# sets minimum bound on score
					score_list[score_list < 0] = 0

					chosen_direction = self._ChooseDirection(score_list)
					chosen_score = score_list[chosen_direction]

					self._UpdateMatrices(
						idx, jdx, chosen_score, chosen_direction
						)



		score, q1, q2 = self.__trace__()
		return (score, q1, q2)

	def __trace__(self):

		"""
		Implementation of traceback for SmithWaterman
		"""

		# Find the highest scoring points on the matrix
		score = np.max(self.score_mat)

		# Select the index furthest down the diagonal
		pos_idx, pos_jdx = np.where(self.score_mat >= score)
		chosen_pos = np.argmax(pos_idx + pos_jdx)
		idx, jdx = (pos_idx[chosen_pos], pos_jdx[chosen_pos])

		q1 = ""
		q2 = ""

		while True:

			if (idx == 0) and (jdx == 0):
				break

			elif self.score_mat[(idx,jdx)] == 0:
				break

			current_direction = self.trace_mat[(idx, jdx)]

			# traverse diagonally
			if current_direction == 0:
				q1 += self.seq1[idx-1]
				q2 += self.seq2[jdx-1]
				idx -= 1
				jdx -= 1

			# traverse left
			elif current_direction == 1:
				q1 += "-"
				q2 += self.seq2[jdx-1]
				jdx -= 1

			# traverse up
			elif current_direction == 2:
				q1 += self.seq1[idx-1]
				q2 += "-"
				idx -= 1

		q1 = q1[::-1]
		q2 = q2[::-1]
		return (score, q1, q2)


class NeedlemanWunsch(PairwiseAligner):

	"""
	Protein Sequencing Alignment using Needleman-Wunsch Global Alignment

	Parameters
	----------

	seq1
		string of reference protein sequence

	seq2
	 	string of query protein sequence

	submat
		a subsitution matrix as a 2D numpy array

	idx_lookup
		a dictionary to map protein characters to their indices in submat

	opening_penalty
		the cost of opening a gap

	extension_penalty
		the cost of extending a previously opened gap
	"""

	def __align__(self):

		"""
		Implementation of alignment for NeedlemanWunsch
		"""

		# row iter
		for idx in np.arange(self.score_mat.shape[0]):

			# col iter
			for jdx in np.arange(self.score_mat.shape[1]):

				# top left corner (double gap)
				if (idx == 0) and (jdx == 0):
					self.score_mat[(idx, jdx)] = 0
					self.trace_mat[(idx, jdx)] = -1

				# top row
				elif (idx == 0) and (jdx > 0):
					self.score_mat[(idx,jdx)] = self._ScoreLeft(idx, jdx)
					self.trace_mat[(idx,jdx)] = 1

				# left column
				elif (idx > 0) and (jdx == 0):
					self.score_mat[(idx, jdx)] = self._ScoreUp(idx, jdx)
					self.trace_mat[(idx,jdx)] = 2

				# main loop
				else:

					score_list = np.array([
						self._ScoreDiag(idx, jdx),
						self._ScoreLeft(idx, jdx),
						self._ScoreUp(idx, jdx)
						])

					chosen_direction = self._ChooseDirection(score_list)
					chosen_score = score_list[chosen_direction]

					self._UpdateMatrices(
						idx, jdx, chosen_score, chosen_direction
						)

		score, q1, q2 = self.__trace__()
		return (score, q1, q2)

	def __trace__(self):

		"""
		Implementation of traceback for NeedlemanWunsch
		"""

		idx, jdx = [i - 1 for i in self.trace_mat.shape]

		q1 = ""
		q2 = ""


		score = self.score_mat[tuple(np.array(self.score_mat.shape) - 1)]


		while True:

			if (idx == 0) and (jdx == 0):
				break

			current_direction = self.trace_mat[(idx, jdx)]

			# traverse diagonally
			if current_direction == 0:
				q1 += self.seq1[idx-1]
				q2 += self.seq2[jdx-1]
				idx -= 1
				jdx -= 1

			# traverse left
			elif current_direction == 1:
				q1 += "-"
				q2 += self.seq2[jdx-1]
				jdx -= 1

			# traverse up
			elif current_direction == 2:
				q1 += self.seq1[idx-1]
				q2 += "-"
				idx -= 1

		q1 = q1[::-1]
		q2 = q2[::-1]
		return (score, q1, q2)



def get_args():

	p = argparse.ArgumentParser()

	p.add_argument(
		"-i", "--input_1", required=True,
		help = "fasta sequence to align"
		)

	p.add_argument(
		"-I", "--input_2", required=True,
		help = "fasta sequence to align"
		)

	p.add_argument(
		"-m", '--method', required=True, type = str,
		help = "Global or Local Alignment [gG[lobal] : global, lL[ocal] : local]"
		)

	p.add_argument(
		"-s", "--sub", required=False,
		default=os.path.join(
			os.path.dirname(os.path.realpath(__file__)), "../scoring_matrices/BLOSUM50.mat"
			),
		help = "path of substitution matrix to use"
		)

	p.add_argument(
		"-g", "--gap_open", required=False, default=11, type=int,
		help = "cost of opening a gap"
		)

	p.add_argument(
		"-e", "--gap_extension", required=False, default=3, type=int,
		help = "cost of extending an existing gap"
		)



	args = p.parse_args()
	return args

def main():

	args = get_args()

	io_obj = io()
	h1, s1 = io_obj.read_fasta(args.input_1)
	h2, s2 = io_obj.read_fasta(args.input_2)
	idx_lookup, submat = io_obj.read_substitution_matrix(args.sub)

	method = SmithWaterman if args.method[0].upper() == "G" else NeedlemanWunsch

	met_obj = method(
		s1, s2, submat, idx_lookup,
		opening_penalty=args.gap_open,
		extension_penalty=args.gap_extension
		)
	score, q1, q2 = met_obj.align()

	print()
	print(q1)
	print(q2)
	print()
	print("Score : {}".format(score))


if __name__ == '__main__':
	main()
