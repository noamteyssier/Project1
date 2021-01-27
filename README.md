# Project 1 - Sequence Alignment
## Due 01/27/2021

![BuildStatus](https://github.com/ucsf-bmi-203-2021/HW1/workflows/HW1/badge.svg?event=push)

### main
To run the code from the command line use the following command :

```
python -m align
```

This requires two input sequences and a decision of the alignment algorithm (global/local)
```
# example alignment (local)
python -m align \
  -i sequences/prot-0004.fa \
  -I sequences/prot-0008.fa \
  -m l

# example alignment (global)
python -m align \
  -i sequences/prot-0004.fa \
  -I sequences/prot-0008.fa \
  -m g

# example local alignment with different matrix
python -m align \
  -i sequences/prot-0004.fa \
  -I sequences/prot-0008.fa \
  -s scoring_matrices/BLOSUM62.mat \
  -m l

# example global alignment with given matrix and non-default gap loss (opening 12 / extension 5)
python -m align \
  -i sequences/prot-0004.fa \
  -I sequences/prot-0008.fa \
  -s scoring_matrices/BLOSUM62.mat \
  -m g -g 12 -e 5

```

### testing
Testing is as simple as running
```
python -m pytest
```
from the root directory of this project.

### Repo Organization

All source code is found in `align/`

Analysis and question answers found in `notebooks/`

Sequencing data in `sequences/`

Scoring matrices found in `scoring_matrices/`
