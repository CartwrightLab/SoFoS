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

#include <unistd.h>

#include "sofos.hpp"
#include "vcf.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <chrono>
#include <numeric>

bool g_sofos_quiet = false;

std::pair<std::string, std::string> timestamp();

double quality_to_p01(int x);
double quality_to_p01(float x);
double phred_to_p01(float x);
double phred_to_p01(int x);

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

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("Combinadic generates the combinatorial number system") {
    SECTION("Combinadic with ploidy 2") {
        Combinadic h{2};
        REQUIRE(h.Get() == 0x3);
        h.Next();
        REQUIRE(h.Get() == 0x5);
        h.Next();
        REQUIRE(h.Get() == 0x6);        
    }
    SECTION("Combinadic with ploidy 3") {
        Combinadic h{3};
        REQUIRE(h.Get() == 0x7);
        h.Next();
        REQUIRE(h.Get() == 0xb);
        h.Next();
        REQUIRE(h.Get() == 0xd);
        h.Next();
        REQUIRE(h.Get() == 0xe);
    }
}
#endif

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
    BcfReader reader{path};

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
    auto header_ptr = reader.header();

    size_t nsamples = bcf_hdr_nsamples(header_ptr);
    if(nsamples <= 0) {
        throw std::invalid_argument("input vcf has no samples.");
    }
    auto char_buffer = make_buffer<char>(64);
    auto int_buffer = make_buffer<int32_t>(1);
    auto gp_buffer = make_buffer<float>(3*nsamples);
    auto gt_buffer = make_buffer<int32_t>(2*nsamples);

    // Construct Histogram
    SofosHistogram histogram{size, alpha, beta};

    // Setup for reading AC tags
    std::vector<int> ac_buffer;
    std::vector<double> af_buffer;
    int ac_which = BCF_UN_FMT;
    // if AN and AC tags are defined, use them
    if(bcf_hdr_id2int(header_ptr, BCF_DT_ID, "AN") >= 0 &&
       bcf_hdr_id2int(header_ptr, BCF_DT_ID, "AC") >= 0 ) {
        ac_which |= BCF_UN_INFO;
    }

    // begin timer
    auto start = std::chrono::steady_clock::now();
    auto last = start;

    // begin reading from the input file
    size_t nsites = 0;
    reader([&](bcf1_t *record, const bcf_hdr_t *header){
        // Make a copy of error_rate because we may update it.
        double local_error_rate = error_rate;
        // identify which allele is ancestral
        // if REFALT is specified, assume 0 is the ancestral allele.
        int anc = 0;
        if(!(flags & SOFOS_FLAG_REFALT)) {
            // If REFALT is not specified, use the AA tag
            int n = get_info_string(header, record, "AA", &char_buffer);
            if(n <= 0) {
                // missing, unknown, or empty tag
                return;
            }
            // If error_rate is missing, try to look for an AAQ tag
            if(local_error_rate == 0.0) {
                n = get_info_int32(header, record, "AAQ", &int_buffer);
                if(n > 0) {
                    local_error_rate = phred_to_p01(int_buffer.data[0]);
                }
            }
            // Compare the AA tag to the refalt allele to find the index of the ancestral allele
            // If this fails, anc = n_allele
            bcf_unpack(record, BCF_UN_STR);
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
            if(bcf_calc_ac(header, record, ac_buffer.data(),
                    ac_which) <= 0) {
                return;
            }
            int n = std::accumulate(ac_buffer.begin(), ac_buffer.end(), 0);
            if(n == 0) {
                return;
            }
            // Number of ancestral copies and derived copies.
            n_anc = (anc < ac_buffer.size()) ? ac_buffer[anc] : 0;
            n_der = n - n_anc;
            n_total = n;
        } else {
            // Calculate counts from GP
            int n = get_format_float(header, record, "GP", &gp_buffer);
            if(n <= 0) {
                // missing, unknown, or empty tag
                return;
            }
            int gp_width = n/nsamples;
            // convert GP buffer if needed
            if(flags & SOFOS_FLAG_PHRED_GP) {
                for(int i=0;i<n;++i) {
                    gp_buffer.data[i] = phred_to_p01(gp_buffer.data[i]);
                }
            }

            // Calculate GT (for ploidy)
            n = get_genotypes(header, record, &gt_buffer);
            if(n <= 0) {
                // missing, unknown, or empty tag
                return;
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
                return;
            }
            n_anc = (anc < af_buffer.size()) ? af_buffer[anc] : 0;
            n_der = n_total - n_anc;
        }

        if(flags & SOFOS_FLAG_FOLDED) {
            // When calculating the folded spectrum, we can't assume that either
            // allele is ancestral. Thus we will update the posterior weighted
            // by the possibility of each occurrence.
            // The probability that an der is ancestral is its frequency.
            local_error_rate = (1.0*n_der)/n_total;
        }

        // Update Posterior and Observed using observed counts
        if(!(flags & SOFOS_FLAG_FOLDED)) {
            if(local_error_rate == 0.0) {
                histogram.AddCounts(n_der, n_anc, 1.0);
            } else {
                // If there is some level of uncertainty about the ancestral allele,
                // weight the two possibilities.
                histogram.AddCounts(n_der, n_anc, 1.0-local_error_rate);
                histogram.AddCounts(n_anc, n_der, local_error_rate);
            }
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
                          << bcf_seqname(header, record) << ":" << record->pos+1
                          << std::endl;
            }
        }
    });

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
        histogram.AddCounts(0.0, ploidy*nsamples, zero);
    }

    // calculate prior assuming no data
    histogram.CalculatePrior();

    // fold histogram if necessary
    if(flags & SOFOS_FLAG_FOLDED) {
        histogram.Fold();
    }

    // spacer
    std::cout << "#\n";

    // output resulting scale
    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
    std::cout << "Number,Prior,Observed,Posterior\n";
    for(int i=0;i< histogram.posterior().size();++i) {
        std::cout << i
                  << "," << histogram.prior()[i]
                  << "," << histogram.observed()[i]
                  << "," << histogram.posterior()[i]
                  << "\n";
    }
    std::cout << std::flush;
    
    return EXIT_SUCCESS;
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
    int k = lround(f);
    if(k > n) {
        k = n;
    }
    (*bins)[k] += weight;
}

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("update_bins maps a/(a+b) to a bin in the range of [0,N]/N") {
    using namespace Catch::literals;

    std::vector<double> bins(4,0.0);

    SECTION("a+b == 0.0") {
        update_bins(0,0,2,&bins);
        REQUIRE(bins.size() == 4);
        CHECK(bins[0] == 0.0);
        CHECK(bins[1] == 0.0);
        CHECK(bins[2] == 0.0);
        CHECK(bins[3] == 0.0);
    }
    SECTION("a == 1 && b == 0") {
        update_bins(1,0,2,&bins);
        REQUIRE(bins.size() == 4);
        CHECK(bins[0] == 0.0);
        CHECK(bins[1] == 0.0);
        CHECK(bins[2] == 0.0);
        CHECK(bins[3] == 2.0);
    }
    SECTION("a == 0 && b == 1") {
        update_bins(0,1,2,&bins);
        REQUIRE(bins.size() == 4);
        CHECK(bins[0] == 2.0);
        CHECK(bins[1] == 0.0);
        CHECK(bins[2] == 0.0);
        CHECK(bins[3] == 0.0);
    }
    SECTION("after many calls") {
        update_bins(0,1,2,&bins);
        update_bins(1,0,2,&bins);
        update_bins(1,2,2,&bins);
        update_bins(1.1,2,2,&bins);
        REQUIRE(bins.size() == 4);
        CHECK(bins[0] == 2.0);
        CHECK(bins[1] == 4.0);
        CHECK(bins[2] == 0.0);
        CHECK(bins[3] == 2.0);
    }    
}
#endif

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

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("fold_histogram folds second half of a vector onto the first") {
    using namespace Catch::literals;
    SECTION("fold histogram with even number of elements") {
        std::vector<double> v = {1,2,3,10,20,30};
        fold_histogram(&v);
        REQUIRE(v.size() == 3);
        CHECK(v[0] == 31.0_a);
        CHECK(v[1] == 22.0_a);
        CHECK(v[2] == 13.0_a);
    }
    SECTION("fold histogram with odd number of elements") {
        std::vector<double> v = {1,2,3,10,20};
        fold_histogram(&v);
        REQUIRE(v.size() == 3);
        CHECK(v[0] == 21.0_a);
        CHECK(v[1] == 12.0_a);
        CHECK(v[2] == 3.0_a);
    }
}
#endif

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


#ifdef SOFOS_UNIT_TESTS
TEST_CASE("quality_to_p01 converts a quality value to a [0,1] value.") {
    using namespace Catch::literals;
    SECTION("quality_to_p01(int)") {
        CHECK(quality_to_p01(0)  == 0.0_a);
        CHECK(quality_to_p01(10) == 0.9_a);
        CHECK(quality_to_p01(20) == 0.99_a);
        CHECK(quality_to_p01(30) == 0.999_a);
        CHECK(quality_to_p01(40) == 0.9999_a);
        CHECK(quality_to_p01(255) == 1.0_a);
    }
    SECTION("quality_to_p01(double)") {
        CHECK(quality_to_p01(0.0f)  == 0.0_a);
        CHECK(quality_to_p01(10.0f) == 0.9_a);
        CHECK(quality_to_p01(20.0f) == 0.99_a);
        CHECK(quality_to_p01(30.0f) == 0.999_a);
        CHECK(quality_to_p01(40.0f) == 0.9999_a);
        CHECK(quality_to_p01(255.0f) == 1.0_a);
    }
}

TEST_CASE("phred_to_p01 converts a phred value to a [0,1] value.") {
    using namespace Catch::literals;
    auto zero = Approx(0.0).margin(std::numeric_limits<float>::epsilon()*100);
    SECTION("phred_to_p01(int)") {
        CHECK(phred_to_p01(0)  == 1.0_a);
        CHECK(phred_to_p01(10) == 0.1_a);
        CHECK(phred_to_p01(20) == 0.01_a);
        CHECK(phred_to_p01(30) == 0.001_a);
        CHECK(phred_to_p01(40) == 0.0001_a);
        CHECK(phred_to_p01(255) == zero);
    }
    SECTION("phred_to_p01(double)") {
        CHECK(phred_to_p01(0.0f)  == 1.0_a);
        CHECK(phred_to_p01(10.0f) == 0.1_a);
        CHECK(phred_to_p01(20.0f) == 0.01_a);
        CHECK(phred_to_p01(30.0f) == 0.001_a);
        CHECK(phred_to_p01(40.0f) == 0.0001_a);
        CHECK(phred_to_p01(255.0f) == zero);
    }    
}
#endif




