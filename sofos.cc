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
    "Copyright (c) 2018 Reed A. Cartwright, PhD <reed@cartwrig.ht>\n"
    "\n";
}

// Globals
bool g_sofos_quiet=false;

// Flags used to control sofos_main
constexpr unsigned int SOFOS_FLAG_DEFAULT=0;
constexpr unsigned int SOFOS_FLAG_FOLDED=1;
constexpr unsigned int SOFOS_FLAG_REFALT=2;
constexpr unsigned int SOFOS_FLAG_USE_GP=4;
constexpr unsigned int SOFOS_FLAG_PHRED_GP=8;

// Function and Class Declarations
int sofos_main(const char *path, double alpha, double beta, int size,
    double error_rate,
    double zero, double ploidy,
    unsigned int flags=SOFOS_FLAG_DEFAULT);

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
struct buffer_t {
    std::unique_ptr<T[],buffer_free_t> data;
    int capacity;
};

template<typename T>
inline
buffer_t<T> make_buffer(int sz) {
    void *p = std::malloc(sizeof(T)*sz);
    if(p == nullptr) {
        throw std::bad_alloc{};
    }
    return {std::unique_ptr<T[],buffer_free_t>{reinterpret_cast<T*>(p)}, sz};
}

int get_info_string(const bcf_hdr_t *header, bcf1_t *record,
    const char *tag, buffer_t<char>* buffer);
int get_info_int32(const bcf_hdr_t *header, bcf1_t *record,
    const char *tag, buffer_t<int32_t>* buffer);

int get_format_float(const bcf_hdr_t *header, bcf1_t *record,
    const char *tag, buffer_t<float>* buffer);
int get_genotypes(const bcf_hdr_t *header, bcf1_t *record,
    buffer_t<int32_t>* buffer);

double quality_to_p01(int x);
double quality_to_p01(float x);
double phred_to_p01(float x);
double phred_to_p01(int x);

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

// Utility class for handling the combinatorial number system used
// in format fields with size=G.
// https://en.wikipedia.org/wiki/Combinatorial_number_system
class Combinadic {
public:
    Combinadic(int ploidy) {
        Reset(ploidy);
    }

    // Reset to the lowest k-combination
    void Reset(int ploidy) {
        data_ = 0;
        for(int j = 0; j < ploidy; ++j) {
            data_ = (data_ << 1) | 1;
        }
    }

    // Permute the k-combination to the next highest value
    void Next() {
        // Gosper's Hack
        unsigned long long u = data_ & -data_;
        unsigned long long v = u + data_;
        data_ = v + (((v^data_)/u)>>2);
    }

    // Access an integer representing a k-combination
    unsigned long long Get() const { return data_; }

protected:
    unsigned long long data_;
};

// the main processing function
int sofos_main(const char *path, double alpha, double beta, int size, double error_rate,
        double zero, double ploidy, unsigned int flags) {
    assert(path != nullptr);
    assert(alpha > 0.0 && beta > 0.0);
    assert(size > 0);
    assert(0.0 <= error_rate && error_rate <= 1.0);
    assert(zero >= 0.0);
    assert(ploidy >= 0.0);

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
    if(error_rate > 0.0) {
        std::cout << "#error_rate=" << error_rate << "\n";
    }
    std::cout << "#folded=" << ((flags & SOFOS_FLAG_FOLDED) ? 1 : 0) << "\n";
    std::cout << "#refalt=" << ((flags & SOFOS_FLAG_REFALT) ? 1 : 0) << "\n";
    if(zero > 0.0) {
        std::cout << "#zero=" << zero << "\n";
        std::cout << "#ploidy=" << ploidy << "\n";
    }
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
    auto char_buffer = make_buffer<char>(64);
    auto int_buffer = make_buffer<int32_t>(1);
    auto gp_buffer = make_buffer<float>(3*nsamples);
    auto gt_buffer = make_buffer<int32_t>(2*nsamples);

    // A sample of N copies can have [0,N] derived alleles.
    std::vector<double> posterior(size+1, 0.0);
    std::vector<double> bins(size+1, 0.0);

    // Setup for reading AC tags
    std::vector<int> ac_buffer;
    std::vector<double> af_buffer;
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
        // Make a copy of error_rate because we may update it.
        double local_error_rate = error_rate;
        // identify which allele is ancestral
        // if REFALT is specified, assume 0 is the ancestral allele.
        int anc = 0;
        if(!(flags & SOFOS_FLAG_REFALT)) {
            // If REFALT is not specified, use the AA tag
            int n = get_info_string(header.get(), record.get(), "AA", &char_buffer);
            if(n <= 0) {
                // missing, unknown, or empty tag
                continue;
            }
            // If error_rate is missing, try to look for an AAQ tag
            if(local_error_rate == 0.0) {
                n = get_info_int32(header.get(), record.get(), "AAQ", &int_buffer);
                if(n > 0) {
                    local_error_rate = phred_to_p01(int_buffer.data[0]);
                }
            }
            // Compare the AA tag to the refalt allele to find the index of the ancestral allele
            // If this fails, anc = n_allele
            bcf_unpack(record.get(), BCF_UN_STR);
            for(anc=0;anc<record->n_allele;++anc) {
                if(strcmp(char_buffer.data.get(), record->d.allele[anc]) == 0) {
                    break;
                }
            }
        }
        double n_der,n_anc,n_total;

        if(!(flags & SOFOS_FLAG_USE_GP)) {
            // Calculate the counts of each allele from the AC/AN tags or directly from GT.
            // Skip line if it fails
            ac_buffer.assign(record->n_allele,0);
            if(bcf_calc_ac(header.get(), record.get(), ac_buffer.data(),
                    ac_which) <= 0) {
                continue;
            }
            int n = std::accumulate(ac_buffer.begin(), ac_buffer.end(), 0);
            if(n == 0) {
                continue;
            }
            // Number of ancestral copies and derived copies.
            n_anc = (anc < ac_buffer.size()) ? ac_buffer[anc] : 0;
            n_der = n - n_anc;
            n_total = n;
        } else {
            // Calculate counts from GP
            int n = get_format_float(header.get(), record.get(), "GP", &gp_buffer);
            if(n <= 0) {
                // missing, unknown, or empty tag
                continue;
            }
            int gp_width = n/nsamples;
            // convert GP buffer if needed
            if(flags & SOFOS_FLAG_PHRED_GP) {
                for(int i=0;i<n;++i) {
                    gp_buffer.data[i] = phred_to_p01(gp_buffer.data[i]);
                }
            }

            // Calculate GT (for ploidy)
            n = get_genotypes(header.get(), record.get(), &gt_buffer);
            if(n <= 0) {
                // missing, unknown, or empty tag
                continue;
            }
            int gt_width = n/nsamples;

            // Allocate vector to estimate allele counts from GP
            af_buffer.assign(record->n_allele,0.0);
            for(int i=0; i<nsamples; i++) {
                // identify ploidy for this sample
                int *p_gt = gt_buffer.data.get() + i*gt_width;
                int sample_ploidy=0;
                for(;sample_ploidy < gt_width; ++sample_ploidy) {
                    if(p_gt[sample_ploidy] == bcf_int32_vector_end) {
                        break;
                    }
                }
                // Setup a sequence of k-permutations
                Combinadic gp_alleles{sample_ploidy};
                int gp_alleles_sz = record->n_allele + sample_ploidy - 1;

                // Step thorough each element of the gp
                float *p_gp = gp_buffer.data.get() + i*gp_width;
                for(int j=0; j < gp_width; ++j, gp_alleles.Next()) {
                    // check for missing and end values
                    if(bcf_float_is_vector_end(p_gp[j])) {
                        break;
                    }
                    if(bcf_float_is_missing(p_gp[j])) {
                        continue;
                    }
                    // step through the alleles for this GP
                    auto alleles = gp_alleles.Get();
                    int pos = 0;
                    for(int h=0; h < gp_alleles_sz; ++h) {
                        // if the lowest bit is 1, we have an allele
                        // if the lowest bit is 0, we move on to the next allele
                        if((alleles&0x1) == 1) {
                            af_buffer[pos] += p_gp[j];
                        } else {
                            pos += 1;
                        }
                        // shit to the next bit
                        alleles = alleles >> 1;
                    }
                }
            }
            n_total = std::accumulate(af_buffer.begin(), af_buffer.end(), 0.0);
            if(n_total == 0.0) {
                continue;
            }
            n_anc = (anc < af_buffer.size()) ? af_buffer[anc] : 0;
            n_der = n_total - n_anc;
        }

        // Update posterior using observed counts, which may be 0 and 0.
        if(!(flags & SOFOS_FLAG_FOLDED)) {
            if(local_error_rate == 0.0) {
                update_counts(alpha + n_der, beta + n_anc, 1.0, &posterior);
                // Update the observed bins based on observed counts, skipping 0 and 0.
                // This maps the observed frequency onto a bin based on size.
                update_bins(n_der, n_anc, 1.0, &bins);
            } else {
                // If there is some level of uncertainty about the ancestral allele,
                // weight the two possibilities.
                // anc is ancestral
                update_counts(alpha + n_der, beta + n_anc, 1.0-local_error_rate, &posterior);
                // der is ancestral
                update_counts(alpha + n_anc, beta + n_der, local_error_rate, &posterior);

                update_bins(n_der, n_anc, 1.0-error_rate, &bins);
                update_bins(n_anc, n_der, error_rate, &bins);
            }
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

            update_bins(n_anc, n_der, f, &bins);
            update_bins(n_der, n_anc, 1.0-f, &bins);
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
    if(zero > 0.0) {
        update_counts(alpha, beta+ploidy*nsamples, zero, &posterior);
        update_bins(0, ploidy*nsamples, zero, &bins);
    }

    // calculate prior assuming no data
    std::vector<double> prior(posterior.size(), 0.0);
    update_counts(alpha, beta, nsites+zero, &prior);

    // if we are outputting a folded histogram, update the vectors
    if(flags & SOFOS_FLAG_FOLDED) {
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
int get_info_string(const bcf_hdr_t *header, bcf1_t *record,
    const char *tag, buffer_t<char>* buffer)
{
    char *p = buffer->data.get();
    int n = bcf_get_info_string(header, record, tag, &p, &buffer->capacity);
    if(n == -4) {
        throw std::bad_alloc{};
    } else if(p != buffer->data.get()) {
        // update pointer
        buffer->data.release();
        buffer->data.reset(p);
    }
    return n;
}

inline
int get_info_int32(const bcf_hdr_t *header, bcf1_t *record,
    const char *tag, buffer_t<int32_t>* buffer)
{
    int32_t *p = buffer->data.get();
    int n = bcf_get_info_int32(header, record, tag, &p, &buffer->capacity);
    if(n == -4) {
        throw std::bad_alloc{};
    } else if(p != buffer->data.get()) {
        // update pointer
        buffer->data.release();
        buffer->data.reset(p);
    }
    return n;
}

inline
int get_format_float(const bcf_hdr_t *header, bcf1_t *record,
    const char *tag, buffer_t<float>* buffer)
{
    float *p = buffer->data.get();
    int n = bcf_get_format_float(header, record, tag, &p, &buffer->capacity);
    if(n == -4) {
        throw std::bad_alloc{};
    } else if(p != buffer->data.get()) {
        // update pointer
        buffer->data.release();
        buffer->data.reset(p);
    }
    return n;
}

inline
int get_genotypes(const bcf_hdr_t *header, bcf1_t *record,
    buffer_t<int32_t>* buffer)
{
    int32_t *p = buffer->data.get();
    int n = bcf_get_genotypes(header, record, &p, &buffer->capacity);
    if(n == -4) {
        throw std::bad_alloc{};
    } else if(p != buffer->data.get()) {
        // update pointer
        buffer->data.release();
        buffer->data.reset(p);
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
    if(a+b == 0.0) {
        return;
    }
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

inline
double quality_to_p01(float x) {
    return -expm1(-x*(M_LN10/10.0));
}

inline
double phred_to_p01(float x) {
    return exp(-x*(M_LN10/10.0));
}

inline
double quality_to_p01(int x) {
    assert(x >= 0);
    static float data[96] = {
        0.000000000, 0.205671765, 0.369042656, 0.498812766, 0.601892829, 0.683772234,
        0.748811357, 0.800473769, 0.841510681, 0.874107459, 0.900000000, 0.920567177,
        0.936904266, 0.949881277, 0.960189283, 0.968377223, 0.974881136, 0.980047377,
        0.984151068, 0.987410746, 0.990000000, 0.992056718, 0.993690427, 0.994988128,
        0.996018928, 0.996837722, 0.997488114, 0.998004738, 0.998415107, 0.998741075,
        0.999000000, 0.999205672, 0.999369043, 0.999498813, 0.999601893, 0.999683772,
        0.999748811, 0.999800474, 0.999841511, 0.999874107, 0.999900000, 0.999920567,
        0.999936904, 0.999949881, 0.999960189, 0.999968377, 0.999974881, 0.999980047,
        0.999984151, 0.999987411, 0.999990000, 0.999992057, 0.999993690, 0.999994988,
        0.999996019, 0.999996838, 0.999997488, 0.999998005, 0.999998415, 0.999998741,
        0.999999000, 0.999999206, 0.999999369, 0.999999499, 0.999999602, 0.999999684,
        0.999999749, 0.999999800, 0.999999842, 0.999999874, 0.999999900, 0.999999921,
        0.999999937, 0.999999950, 0.999999960, 0.999999968, 0.999999975, 0.999999980,
        0.999999984, 0.999999987, 0.999999990, 0.999999992, 0.999999994, 0.999999995,
        0.999999996, 0.999999997, 0.999999997, 0.999999998, 0.999999998, 0.999999999,
        0.999999999, 0.999999999, 0.999999999, 0.999999999, 1.000000000, 1.000000000
    };
    return (x < 96) ? data[x] : 1.0;
}

inline
double phred_to_p01(int x) {
    assert(x >= 0);
    static float data[96] = {
        1.000000000, 0.794328235, 0.630957344, 0.501187234, 0.398107171, 0.316227766,
        0.251188643, 0.199526231, 0.158489319, 0.125892541, 0.100000000, 0.079432823,
        0.063095734, 0.050118723, 0.039810717, 0.031622777, 0.025118864, 0.019952623,
        0.015848932, 0.012589254, 0.010000000, 0.007943282, 0.006309573, 0.005011872,
        0.003981072, 0.003162278, 0.002511886, 0.001995262, 0.001584893, 0.001258925,
        0.001000000, 0.000794328, 0.000630957, 0.000501187, 0.000398107, 0.000316228,
        0.000251189, 0.000199526, 0.000158489, 0.000125893, 0.000100000, 0.000079433,
        0.000063096, 0.000050119, 0.000039811, 0.000031623, 0.000025119, 0.000019953,
        0.000015849, 0.000012589, 0.000010000, 0.000007943, 0.000006310, 0.000005012,
        0.000003981, 0.000003162, 0.000002512, 0.000001995, 0.000001585, 0.000001259,
        0.000001000, 0.000000794, 0.000000631, 0.000000501, 0.000000398, 0.000000316,
        0.000000251, 0.000000200, 0.000000158, 0.000000126, 0.000000100, 0.000000079,
        0.000000063, 0.000000050, 0.000000040, 0.000000032, 0.000000025, 0.000000020,
        0.000000016, 0.000000013, 0.000000010, 0.000000008, 0.000000006, 0.000000005,
        0.000000004, 0.000000003, 0.000000003, 0.000000002, 0.000000002, 0.000000001,
        0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000000, 0.000000000
    };
    return (x < 96) ? data[x] : 0.0;
}
