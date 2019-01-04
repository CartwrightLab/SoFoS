# SoFoS
Rescale genetic polymorphism data to match a common sample size.


## Usage

`sofos [OPTION]... [FILE] > [OUTPUT]`

With no FILE or when FILE is -, read standard input.

  -a number -b number  shape parameters of beta prior
  -n integer           number of gene copies in posterior resample
  -f -u                generated (f)olded or (u)nfolded distributions
  -t -r                use AA (t)ag or (r)eference allele as ancestral
  -z number            add extra invariant sites to manage ascertainment bias
  -q -v                (q)uiet progress info or be (v)erbose
  -h                   print usage information

Default: sofos -f -a 1.0 -b 1.0 -n 9
Note: Unless otherwise stated -f enables -r and -u enables -t.
