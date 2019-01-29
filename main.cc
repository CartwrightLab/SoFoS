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

#include <unistd.h>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>

// print a usage message for sofos
void print_usage(const char* exe, std::ostream& os) {
    const char* p = std::strrchr(exe, '/');
    if(p != nullptr && p[0] != '\0' && p[1] != '\0') {
        exe = p + 1;
    }
    os << "Usage: " << exe
       << " [OPTION]... [FILE] > [OUTPUT]\n"
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
          "Default: "
       << exe
       << " -f -a 1.0 -b 1.0 -n 9\n"
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
int main(int argc, char* argv[]) {
    try {
        // default parameters
        sofos_params_t params;

        int refalt = -1;

        int use_gp = 0;

        // Process program options via getopt
        char c = 0;
        while((c = getopt(argc, argv, "a:b:n:e:hufrtpz:P:qv")) != -1) {
            switch(c) {
            case 'a':
                params.alpha = std::stod(optarg);
                break;
            case 'b':
                params.beta = std::stod(optarg);
                break;
            case 'n':
                params.size = std::stoi(optarg);
                break;
            case 'e':
                params.error_rate = std::stod(optarg);
                break;
            case 'z':
                params.zero_count = std::stod(optarg);
                break;
            case 'P':
                params.ploidy = std::stod(optarg);
                break;
            case 'u':
                params.flag_folded = false;
                break;
            case 'f':
                params.flag_folded = true;
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
        // Setup flags for sofos_main
        params.flag_refalt = (refalt == 1 || (refalt != 0 && params.flag_folded));
        params.flag_use_gp = (use_gp > 0);
        params.flag_phred_gp = (use_gp > 1);

        Sofos sofos{params};

        // read input from a file or fall back to stdin
        std::vector<const char*> paths(argv + optind, argv + argc);
        if(paths.empty()) {
            paths.push_back("-");
        }

        output_header(std::cout, sofos, paths);

        for(auto&& path : paths) {
            sofos.RescaleBcf(path);
        }
        sofos.FinishHistogram();

        output_body(std::cout, sofos);

        return EXIT_SUCCESS;

    } catch(std::exception& e) {
        // If an exception is thrown, print it to stderr.
        std::cerr << "ERROR: " << e.what() << std::endl;
    }
    return EXIT_FAILURE;
}
#endif
