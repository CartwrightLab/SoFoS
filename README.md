# SoFoS
Rescale genetic polymorphism data to match a common sample size.


## Usage

`sofos [OPTION]... [FILE] > [OUTPUT]`

With no FILE or when FILE is -, read standard input.

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

Defaults: `sofos -f -a 1.0 -b 1.0 -n 9 -d 2 -e 0`

### Notes

 - Unless otherwise stated -f enables -r and -u enables -t.
 - -p specifies that GP contains probabilities in the range 0 and 1.
 - -pp specifies that GP contains phred-scaled probabilities.
 - -e is only used for generating unfolded spectra.
