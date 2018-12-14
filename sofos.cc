/*****************************************************************************
Copyright (c) 2010,2018 Reed A. Cartwright, PhD <reed@cartwrig.ht>

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

#include <unistd.h>

#include <htslib/vcf.h>

// bcf_get_format_* uses realloc internally, so this buffer
// will be managed by malloc and free
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

inline
int get_genotypes(const bcf_hdr_t *header, bcf1_t *record,
    buffer_t<int32_t>* buffer, int *capacity)
{
    int *p = buffer->get();
    int n = bcf_get_genotypes(header, record, &p, capacity);
    if(n == -4) {
        throw std::bad_alloc{};
    } else if(p != buffer->get()) {
        // update pointer
        buffer->release();
        buffer->reset(p);
    }
    return n;
}

inline
bool is_ref_missing(const bcf1_t* record) {
    if(record->n_allele == 0) {
        return true;
    }
    const char* ref = record->d.allele[0];
    size_t len = strlen(ref);
    if(len == 0) {
        return true;
    }
    if(len == 1 && (ref[0] == '.' || ref[0] == 'N' || ref[0] == 'n')) {
        return true;
    }
    return false;
}

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

void update_counts(double a, double b, std::vector<double> *counts) {
    assert(counts != nullptr);
    using std::lgamma;
    double ab = a+b;
    int size = counts->size()-1;
    double scale = lgamma(a)+lgamma(b)+lgamma(ab+size)-lgamma(ab)-lgamma(size+1);
    for(int k=0;k<=size;++k) {
        // rescale assuming beta-binomial model
        double d = lgamma(a+k)+lgamma(b+size-k);
        d -= lgamma(k+1) + lgamma(size-k+1);
        d -= scale;
        (*counts)[k] += exp(d);
    }
}

int sofos_main(const char *path, double alpha, double beta, int size) {
    assert(path != nullptr);
    vcfFile *input_ptr = vcf_open(path,"r");
    if(input_ptr == nullptr) {
        throw std::runtime_error(std::string{"unable to open input file: '"} + path +"'.");
    }

    auto stamp = timestamp();
    std::cout << "#SoFoS v2.0\n";
    std::cout << "#date=" << stamp.first << "\n";
    std::cout << "#epoch=" << stamp.second << "\n";

    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::cout << "#path=" << path << "\n";
    std::cout << "#alpha=" << alpha << "\n";
    std::cout << "#beta=" << beta << "\n";
    std::cout << "#size=" << size << "\n";

    // spacer
    std::cout << "#\n";

    std::unique_ptr<vcfFile,file_free_t> input{input_ptr};
    std::unique_ptr<bcf_hdr_t,header_free_t> header{bcf_hdr_read(input.get())};
    std::unique_ptr<bcf1_t,bcf_free_t> record{bcf_init()};
    
    size_t nsamples = bcf_hdr_nsamples(header.get());
    int capacity = nsamples*4;
    auto buffer = make_buffer<int32_t>(capacity);

    // A sample of N copies can have [0,N] non-ref alleles.
    std::vector<double> counts(size+1, 0.0);

    size_t nsites = 0;

    while(bcf_read(input.get(), header.get(), record.get()) == 0) {
        bcf_unpack(record.get(), BCF_UN_STR);

        //std::cerr << record.get()->rid << ":" << record.get()->pos << std::endl;

        if(is_ref_missing(record.get())) {
            continue;
        }
        int n_ref = 0;
        int n_alt = 0;

        int n = get_genotypes(header.get(), record.get(), &buffer, &capacity);
        int max_ploidy = n/nsamples;
        for(int sample = 0; sample < nsamples; ++sample) {
            int32_t *gt = buffer.get()+sample*max_ploidy;
            for(int i=0; i < max_ploidy; ++i) {
                if(gt[i] == bcf_int32_vector_end) {
                    break;
                }
                if(bcf_gt_is_missing(gt[i])) {
                    continue;
                }
                int allele = bcf_gt_allele(gt[i]);
                n_ref += (allele == 0);
                n_alt += (allele != 0);
            }
        }
        int n_total = n_ref+n_alt;
        // Question: Should sites with no data be included in posterior?
        if(n_total == 0) {
            continue;
        }
        update_counts(alpha + n_alt, beta + n_ref, &counts);
        nsites += 1;
    }
    // calculate prior
    std::vector<double> prior(counts.size(), 0.0);
    update_counts(alpha, beta, &prior);

    // output resulting scale
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::cout << "Number,Posterior,Prior\n";
    for(int i=0;i<=size;++i) {
        std::cout << i << "," << counts[i] << "," << nsites*prior[i] << "\n";
    }
    std::cout << std::flush;
    
    return EXIT_SUCCESS;
}

void print_usage(const char* exe, std::ostream& os) {
    os << "Usage: " << exe << " -a alpha_param -b beta_param -n rescale_size [filepath] > output.csv\n";
}

int main(int argc, char *argv[]) {
    try {
        char c = 0;
        double alpha = 1.0;
        double beta = 1.0;
        int size = 9;
        while((c = getopt(argc, argv, "a:b:n:h")) != -1) {
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
            case 'h':
                print_usage(argv[0], std::cout);
                return EXIT_SUCCESS;
            case '?':
                print_usage(argv[0], std::cerr);
                return EXIT_FAILURE;                
            };
        }
        const char *path = nullptr;
        if(optind == argc) {
            path = "-";
        } else if(optind+1 == argc) {
            path = argv[optind];
        } else {
            throw std::invalid_argument("more than one input file was specified.");
        }
        return sofos_main(path, alpha, beta, size);

    } catch(std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
    }
    return EXIT_FAILURE;
}