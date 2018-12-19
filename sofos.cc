/*****************************************************************************
Copyright (c) 2018 Reed A. Cartwright, PhD <reed@cartwrig.ht>

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

#include <memory>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <chrono>
#include <numeric>

#include <unistd.h>

#include <htslib/vcf.h>
#include <htslib/vcfutils.h>

// print a usage message for sofos
void print_usage(const char* exe, std::ostream& os) {
    const char* p = strrchr(exe, '/');
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
    "  -z number            add extra invariant sites to manage ascertainment bias\n"
    "  -q -v                (q)uiet progress info or be (v)erbose\n"
    "  -h                   print usage information\n"
    "\n"
    "Default: " << exe << " -f -a 1.0 -b 1.0 -n 9\n"
    "Note: Unless otherwise stated -f enables -r and -u enables -t.\n"
    "\n"
    "Copyright (c) 2018 Reed A. Cartwright, PhD <reed@cartwrig.ht>\n"
    "\n";
}

// Globals
bool g_sofos_quiet=false;

// Flags used to control sofos_main
constexpr unsigned int SOFOS_FLAG_DEFAULT=0;
constexpr unsigned int SOFOS_FLAG_UNFOLDED=1;
constexpr unsigned int SOFOS_FLAG_REFALT=2;

// Function and Class Declarations
int sofos_main(const char *path, double alpha, double beta, int size,
    double zero, unsigned int flags=SOFOS_FLAG_DEFAULT);
std::pair<std::string, std::string> timestamp();
bool is_ref_missing(bcf1_t* record);
bool is_allele_missing(const char* a);
void update_counts(double a, double b, double weight, std::vector<double> *counts);
void update_bins(double a, double b, double weight, std::vector<double> *counts);
void fold_histogram(std::vector<double> *counts);

// The *_free_t classes are used enable RAII on pointers created by htslib.
struct buffer_free_t {
    void operator()(void* ptr) const {
        free(ptr);
    }
};
struct file_free_t{
    void operator()(void* ptr) const {
        vcf_close(reinterpret_cast<vcfFile*>(ptr));
    }
};
struct header_free_t{
    void operator()(void* ptr) const {
        bcf_hdr_destroy(reinterpret_cast<bcf_hdr_t*>(ptr));
    }
};
struct bcf_free_t {
    void operator()(void* ptr) const {
        bcf_destroy(reinterpret_cast<bcf1_t*>(ptr));
    }
};

// Templates and functions for handling buffers used by htslib
template<typename T>
using buffer_t = std::unique_ptr<T[],buffer_free_t>;

template<typename T>
inline
buffer_t<T> make_buffer(std::size_t sz) {
    void *p = std::malloc(sizeof(T)*sz);
    if(p == nullptr) {
        throw std::bad_alloc{};
    }
    return buffer_t<T>{ reinterpret_cast<T*>(p) };
}

int get_string(const bcf_hdr_t *header, bcf1_t *record,
    const char *tag, buffer_t<char>* buffer, int *capacity);

// Main program entry point
int main(int argc, char *argv[]) {
    try {

        // default parameters
        double alpha = 1.0;
        double beta = 1.0;
        double zero = 0.0;
        int size = 9;
        bool folded = true;
        int refalt = -1;

        // Process program options via getopt
        char c = 0;
        while((c = getopt(argc, argv, "a:b:n:hufrtz:qv")) != -1) {
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
            case 'z':
                zero = std::stod(optarg);
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
        flags |= (!folded) ? SOFOS_FLAG_UNFOLDED : 0; 
        flags |= (refalt == 1 || (refalt != 0 && folded)) ? SOFOS_FLAG_REFALT : 0;
        return sofos_main(path, alpha, beta, size, zero, flags);

    } catch(std::exception &e) {
        // If an exception is thrown, print it to stderr.
        std::cerr << "ERROR: " << e.what() << std::endl;
    }
    return EXIT_FAILURE;
}

// the main processing function
int sofos_main(const char *path, double alpha, double beta, int size, double zero, unsigned int flags) {
    assert(path != nullptr);

    // open the input file or throw on error
    vcfFile *input_ptr = vcf_open(path,"r");
    if(input_ptr == nullptr) {
        throw std::runtime_error(std::string{"unable to open input file: '"} + path +"'.");
    }

    // Output header including program runtime information and options
    auto stamp = timestamp();
    std::cout << "#SoFoS v2.0\n";
    std::cout << "#date=" << stamp.first << "\n";
    std::cout << "#epoch=" << stamp.second << "\n";

    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::cout << "#path=" << path << "\n";
    std::cout << "#alpha=" << alpha << "\n";
    std::cout << "#beta=" << beta << "\n";
    std::cout << "#size=" << size << "\n";
    std::cout << "#folded=" << ((flags & SOFOS_FLAG_UNFOLDED) ? 0 : 1) << "\n";
    std::cout << "#refalt=" << ((flags & SOFOS_FLAG_REFALT) ? 1 : 0) << "\n";

    // Setup for reading from vcf file
    std::unique_ptr<vcfFile,file_free_t> input{input_ptr};
    std::unique_ptr<bcf_hdr_t,header_free_t> header{bcf_hdr_read(input.get())};
    if(!header) {
        throw std::invalid_argument("unable to read header from input.");
    }
    std::unique_ptr<bcf1_t,bcf_free_t> record{bcf_init()};
    if(!record) {
        throw std::invalid_argument("unable to initialize empty vcf record.");
    }
    size_t nsamples = bcf_hdr_nsamples(header.get());
    if(nsamples <= 0) {
        throw std::invalid_argument("input vcf has no samples.");
    }
    int char_capacity = 64;
    auto char_buffer = make_buffer<char>(char_capacity);

    // A sample of N copies can have [0,N] derived alleles.
    std::vector<double> posterior(size+1, 0.0);
    std::vector<double> bins(size+1, 0.0);

    // Setup for reading AC tags
    std::vector<int> ac_buffer;
    int ac_which = BCF_UN_FMT;
    // if AN and AC tags are defined, use them
    if(bcf_hdr_id2int(header.get(), BCF_DT_ID, "AN") >= 0 &&
       bcf_hdr_id2int(header.get(), BCF_DT_ID, "AC") >= 0 ) {
        ac_which |= BCF_UN_INFO;
    }

    // begin timer
    auto start = std::chrono::steady_clock::now();
    auto last = start;

    // begin reading from the input file
    size_t nsites = 0;
    while(bcf_read(input.get(), header.get(), record.get()) == 0) {
        // identify which allele is ancestral
        // if REFALT is specified, assume 0 is the ancestral allele.
        int anc = 0;
        if(!(flags & SOFOS_FLAG_REFALT)) {
            // If REFALT is not specified, use the AA tag
            int n = get_string(header.get(), record.get(), "AA", &char_buffer, &char_capacity);
            if(n <= 0) {
                // missing, unknown, or empty tag
                continue;
            }
            // Compare the AA tag to the refalt allele to find the index of the ancestral allele
            // If this fails, anc = n_allele
            bcf_unpack(record.get(), BCF_UN_STR);
            for(anc=0;anc<record->n_allele;++anc) {
                if(strcmp(char_buffer.get(), record->d.allele[anc]) == 0) {
                    break;
                }
            }
        }

        // Calculate the counts of each allele from the AC/AN tags or directly from GT.
        // Skip line if it fails
        ac_buffer.assign(record->n_allele,0.0);
        if(bcf_calc_ac(header.get(), record.get(), ac_buffer.data(),
                ac_which) == 0) {
            continue;
        }
        int n_total = std::accumulate(ac_buffer.begin(), ac_buffer.end(), 0);
        if(n_total == 0) {
            continue;
        }
        // Number of ancestral copies and derived copies.
        int n_anc = (anc < ac_buffer.size()) ? ac_buffer[anc] : 0;
        int n_der = n_total - n_anc;

        // Update posterior using observed counts, which may be 0 and 0.
        if(flags & SOFOS_FLAG_UNFOLDED) {
            update_counts(alpha + n_der, beta + n_anc, 1.0, &posterior);
        } else {
            // When calculating the folded spectrum, we can't assume that either
            // allele is ancestral. Thus we will update the posterior weighted
            // by the possibility of each occurrence.
            // The probability that an allele is ancestral is its frequency.

            double f = (1.0*n_der)/n_total;
            // der is ancestral
            update_counts(alpha + n_anc, beta + n_der, f, &posterior);
            // anc is ancestral
            update_counts(alpha + n_der, beta + n_anc, 1.0-f, &posterior);
        }

        // Update the observed bins based on observed counts, skipping 0 and 0.
        // This maps the observed frequency onto a bin based on size.
        if(n_der+n_anc > 0) {
            update_bins(n_der, n_anc, 1.0, &bins);
        }

        // Increment the number of observed sites.
        nsites += 1;

        // Output a progress message every minute or so.
        if(nsites % 1000 == 0 && !g_sofos_quiet) {
            auto now = std::chrono::steady_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - last);
            if(elapsed.count() >= 60) {
                elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start);
                last = now;
                std::cerr << "## [" << elapsed.count() << "s elapsed] "
                          << nsites << " sites processed --- at "
                          << bcf_seqname(header.get(), record.get()) << ":" << record->pos+1
                          << std::endl;
            }
        }
    }
    // Output a progress message at end.
    if(!g_sofos_quiet) {
        auto now = std::chrono::steady_clock::now();
        auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - start);
        std::cerr << "## [" << elapsed.count() << "s elapsed] "
                  << nsites << " sites processed"
                  << std::endl;
    }

    // Add zero-count sites (for controlling ascertainment bias)
    update_counts(alpha, beta+nsamples, zero, &posterior);
    update_bins(0, nsamples, zero, &bins);

    // calculate prior assuming no data
    std::vector<double> prior(posterior.size(), 0.0);
    update_counts(alpha, beta, nsites+zero, &prior);

    // if we are outputting a folded histogram, update the vectors
    if(!(flags & SOFOS_FLAG_UNFOLDED)) {
        fold_histogram(&prior);
        fold_histogram(&bins);
        fold_histogram(&posterior);
    }

    // spacer
    std::cout << "#\n";

    // output resulting scale
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::cout << "Number,Prior,Observed,Posterior\n";
    for(int i=0;i<posterior.size();++i) {
        std::cout << i
                  << "," << prior[i]
                  << "," << bins[i]
                  << "," << posterior[i]
                  << "\n";
    }
    std::cout << std::flush;
    
    return EXIT_SUCCESS;
}

// htslib may call realloc on our pointer. When using a managed buffer,
// we need to check to see if it needs to be updated.
inline
int get_string(const bcf_hdr_t *header, bcf1_t *record,
    const char *tag, buffer_t<char>* buffer, int *capacity)
{
    char *p = buffer->get();
    int n = bcf_get_info_string(header, record, tag, &p, capacity);
    if(n == -4) {
        throw std::bad_alloc{};
    } else if(p != buffer->get()) {
        // update pointer
        buffer->release();
        buffer->reset(p);
    }
    return n;
}

// an allele is missing if its value is '.', 'N', or 'n'.
inline
bool is_allele_missing(const char* a) {
    if(a == nullptr) {
        return true;
    }
    if(a[0] == '\0') {
        return true;
    }
    if((a[0] == '.' || a[0] == 'N' || a[0] == 'n') && a[1] == '\0') {
        return true;
    }
    return false;
}

// determine if the reference allele is missing
inline
bool is_ref_missing(bcf1_t* record) {
    if(record->n_allele == 0) {
        return true;
    }
    bcf_unpack(record, BCF_UN_STR);
    return is_allele_missing(record->d.allele[0]);
}

// Generate a timestamp. Returns an ISO formatted timestring
// and seconds since epoch.
std::pair<std::string, std::string> timestamp() {
    using namespace std;
    using namespace std::chrono;
    std::string buffer(127, '\0');
    auto now = system_clock::now();
    auto now_t = system_clock::to_time_t(now);
    size_t sz = strftime(&buffer[0], 127, "%FT%T%z",
                         localtime(&now_t));
    buffer.resize(sz);
    auto epoch = std::chrono::duration_cast<std::chrono::milliseconds>(
                     now.time_since_epoch());
    return {buffer, to_string(epoch.count())};
}

// This is the core algorithm of SoFoS.
// It calculates the histogram of Beta-Binomial(N, a, b) and
// adds that histogram to counts based on weight.
//
//                  N!      Gamma[a+k]   Gamma[b+N-k]    Gamma[a+b]
// f(k;N,a,b) =  -------- * ---------- * ------------ * ------------
//               k!(N-k)!    Gamma[a]      Gamma[b]     Gamma[a+b+N]
//
// To generated the posterior, a = alpha+n_der and b = beta+n_anc
//
void update_counts(double a, double b, double weight, std::vector<double> *counts) {
    assert(counts != nullptr);
    assert(a > 0.0 && b > 0.0);
    using std::lgamma;
    double ab = a+b;
    int size = counts->size()-1;
    double scale = lgamma(a)+lgamma(b)+lgamma(ab+size)-lgamma(ab)-lgamma(size+1);
    for(int k=0;k<=size;++k) {
        // rescale assuming beta-binomial model
        double d = lgamma(a+k)+lgamma(b+size-k);
        d -= lgamma(k+1) + lgamma(size-k+1);
        d -= scale;
        (*counts)[k] += weight*exp(d);
    }
}

// maps a/(a+b) to a bin in the range of [0,N]/N
// rounds to nearest bin. Thus the bins on the edge may
// have lower counts than those in the middle.
void update_bins(double a, double b, double weight, std::vector<double> *bins) {
    assert(bins != nullptr);
    assert(a >= 0.0 && b >= 0.0);
    assert(a+b > 0);
    // bins are from [0,n]
    int n = bins->size()-1;
    double f = a/(a+b)*n;
    int k = static_cast<int>(f+0.5);
    if(k > n) {
        k = n;
    }
    (*bins)[k] += weight;
}

// folds the histogram so that the second half is added to the first half.
// handles both odd and even vector sizes.
inline
void fold_histogram(std::vector<double> *counts) {
    assert(counts != nullptr);
    for(int k=0;k<counts->size()/2;++k) {
        (*counts)[k] += (*counts)[counts->size()-k-1];
    }
    counts->resize((counts->size()+1)/2);
}
