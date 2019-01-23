/*****************************************************************************
Copyright (c) 2019 Reed A. Cartwright, PhD <reed@cartwrig.ht>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*****************************************************************************/

#ifdef SOFOS_UNIT_TESTS
#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#endif

#include "sofos.hpp"

#include <string>
#include <cstring>
#include <stdexcept>
#include <iostream>

#include <unistd.h>

// print a usage message for sofos
void print_usage(const char* exe, std::ostream& os) {
    const char* p = std::strrchr(exe, '/');
    if(p != nullptr && p[0] != '\0' && p[1] != '\0') {
        exe = p+1;
    }
    os << "Usage: " << exe << " [OPTION]... [FILE] > [OUTPUT]\n"
    "Rescale genetic polymorphism data to match a common sample size.\n"
    "\n"
    "With no FILE or when FILE is -, read standard input.\n"
    "\n"
    "  -a number -b number  shape parameters of beta prior\n"
    "  -n integer           number of gene copies in posterior resample\n"
    "  -f -u                generated (f)olded or (u)nfolded distributions\n"
    "  -t -r                use AA (t)ag or (r)eference allele as ancestral\n"
    "  -e number            probability of ancestral allele misassignment\n"
    "  -p [or] -pp          use GP tag to estimate allele frequencies\n"
    "  -z number            add extra invariant sites to manage ascertainment bias\n"
    "  -P number            average ploidy of samples (used with -z)\n"
    "  -q -v                (q)uiet progress info or be (v)erbose\n"
    "  -h                   print usage information\n"
    "\n"
    "Default: " << exe << " -f -a 1.0 -b 1.0 -n 9\n"
    "Notes: Unless otherwise stated -f enables -r and -u enables -t.\n"
    "       -p specifies that GP contains probabilities in the range 0 and 1.\n"
    "       -pp specifies that GP contains phred-scaled probabilities.\n"
    "       -e is only used when generating unfolded spectra.\n"
    "\n"
    "Copyright (c) 2019 Reed A. Cartwright, PhD <reed@cartwrig.ht>\n"
    "\n";
}

// Main program entry point
#ifndef CATCH_CONFIG_MAIN
int main(int argc, char *argv[]) {
    try {
        // default parameters
        double alpha = 1.0;
        double beta = 1.0;
        double zero = 0.0;
        int size = 9;
        bool folded = true;
        int refalt = -1;
        double ploidy = 2;
        double error_rate = 0.0;
        int use_gp = 0;

        // Process program options via getopt
        char c = 0;
        while((c = getopt(argc, argv, "a:b:n:e:hufrtpz:P:qv")) != -1) {
            switch(c) {
            case 'a':
                alpha = std::stod(optarg);
                break;
            case 'b':
                beta = std::stod(optarg);
                break;
            case 'n':
                size =  std::stoi(optarg);
                break;
            case 'e':
                error_rate = std::stod(optarg);
                break;
            case 'z':
                zero = std::stod(optarg);
                break;
            case 'P':
                ploidy = std::stod(optarg);
                break;
            case 'u':
                folded = false;
                break;
            case 'f':
                folded = true;
                break;
            case 'r':
                refalt = 1;
                break;
            case 't':
                refalt = 0;
                break;
            case 'p':
                ++use_gp;
                break;
            case 'q':
                g_sofos_quiet = true;
                break;
            case 'v':
                g_sofos_quiet = false;
                break;
            case 'h':
                print_usage(argv[0], std::cout);
                return EXIT_SUCCESS;
            case '?':
                print_usage(argv[0], std::cerr);
                return EXIT_FAILURE;                
            };
        }
        // read input from a file or fall back to stdin
        const char *path = nullptr;
        if(optind == argc) {
            path = "-";
        } else if(optind+1 == argc) {
            path = argv[optind];
        } else {
            throw std::invalid_argument("more than one input file was specified.");
        }
        // Setup flags for sofos_main
        unsigned int flags = SOFOS_FLAG_DEFAULT;
        flags |= (folded) ? SOFOS_FLAG_FOLDED : 0; 
        flags |= (refalt == 1 || (refalt != 0 && folded)) ? SOFOS_FLAG_REFALT : 0;
        flags |= (use_gp > 0) ? SOFOS_FLAG_USE_GP : 0;
        flags |= (use_gp > 1) ? SOFOS_FLAG_PHRED_GP : 0;

        return sofos_main(path, alpha, beta, size, error_rate, zero, ploidy, flags);

    } catch(std::exception &e) {
        // If an exception is thrown, print it to stderr.
        std::cerr << "ERROR: " << e.what() << std::endl;
    }
    return EXIT_FAILURE;
}
#endif
