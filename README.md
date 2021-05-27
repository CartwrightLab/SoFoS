# SoFoS
Rescale and smooth genetic polymorphism data to match a common sample size.

[![CircleCI](https://circleci.com/gh/CartwrightLab/SoFoS.svg?style=svg)](https://circleci.com/gh/CartwrightLab/SoFoS)
[![codecov](https://codecov.io/gh/CartwrightLab/SoFoS/branch/master/graph/badge.svg)](https://codecov.io/gh/CartwrightLab/SoFoS)

## Usage

`sofos [OPTION]... [FILE] > [OUTPUT]`

With no FILE or when FILE is -, read standard input.

Input format is a VCF/BCF file with genotypes.

| Argument            |Description                   |
|---------------------|------------------------------|
|-a number -b number  |shape parameters of beta prior|
|-n integer           |number of gene copies in posterior resample|
|-f -u                |generated (f)olded or (u)nfolded distributions|
|-t -r                |use AA (t)ag or (r)eference allele as ancestral|
|-e number            |probability of ancestral allele misassignment|
|-p [or] -pp          |use GP tag to estimate allele frequencies|
|-z number            |add extra invariant sites to manage ascertainment bias|
|-P number			  |average ploidy of samples (used with -z)|
|-q -v                |(q)uiet progress info or be (v)erbose|
|-h                   |print usage information|

Defaults: `sofos -u -a 1.0 -b 1.0 -n 10 -P 2`

### Notes

 - Unless otherwise stated -f enables -r and -u enables -t.
 - -p specifies that GP contains probabilities in the range 0 and 1.
 - -pp specifies that GP contains phred-scaled probabilities.
 - -e is only used for generating unfolded spectra.
