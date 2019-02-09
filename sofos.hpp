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
#ifndef SOFOS_SOFOS_HPP
#define SOFOS_SOFOS_HPP

#include "vcf.hpp"

#include <array>
#include <cassert>
#include <cmath>
#include <numeric>
#include <ostream>
#include <vector>

// Globals
extern bool g_sofos_quiet;

struct sofos_params_t {
    double alpha{1.0};
    double beta{1.0};
    int size{10};
    double error_rate{0.0};
    double zero_count{0.0};
    double ploidy{2.0};

    bool flag_folded{false};
    bool flag_refalt{false};
    bool flag_use_gp{false};
    bool flag_phred_gp{false};
};

void update_counts(double a, double b, double weight, std::vector<double> *counts);
void update_bins(double a, double b, double weight, std::vector<double> *counts);
void fold_histogram(std::vector<double> *counts);

class SofosHistogram {
   public:
    SofosHistogram(int size, double alpha, double beta);

    void AddCounts(double num_derived_copies, double num_ancestral_copies, double weight);

    void CalculatePrior();

    void Fold();

    void ZeroCounts();

    size_t num_rows() const {
        assert(prior_.size() == posterior_.size());
        assert(observed_.size() == posterior_.size());

        return posterior_.size();
    }

    std::array<double, 3> row(size_t i) const {
        assert(i < num_rows());
        return {prior_[i], observed_[i], posterior_[i]};
    }
    const std::vector<double> &col(size_t i) const {
        assert(i < 3);
        if(i == 0) {
            return prior_;
        } else if(i == 1) {
            return observed_;
        } else {
            return posterior_;
        }
    }

   private:
    int size_;
    double alpha_;
    double beta_;

    std::vector<double> prior_;
    std::vector<double> observed_;
    std::vector<double> posterior_;
};

class Sofos {
   public:
    explicit Sofos(const sofos_params_t &params) : params_{params}, histogram_{params.size, params.alpha, params.beta} {
        assert(0.0 <= params.error_rate && params.error_rate <= 1.0);
        assert(params.zero_count >= 0.0);
        assert(params.ploidy >= 0.0);
    }

    void RescaleBcf(const char *path);

    void ResetHistogram() { histogram_.ZeroCounts(); }

    void FinishHistogram();

    const SofosHistogram histogram() const { return histogram_; }

    // This is a template in case SofosHistogram::AddCounts becomes an overloaded function
    template <typename T>
    void AddCounts(const std::vector<T> &counts, int anc, double error_rate);

   private:
    sofos_params_t params_;

    SofosHistogram histogram_;

   public:
    friend void output_header(std::ostream &os, const Sofos &sofos, const std::vector<const char *> &paths);
    friend void output_body(std::ostream &os, const Sofos &sofos);
};

inline SofosHistogram::SofosHistogram(int size, double alpha, double beta) : size_{size}, alpha_{alpha}, beta_{beta} {
    assert(size > 0);
    assert(alpha > 0.0 && beta > 0.0);

    ZeroCounts();
}

inline void SofosHistogram::ZeroCounts() {
    // A sample of N copies can have [0,N] derived alleles.
    prior_.assign(size_ + 1, 0.0);
    observed_.assign(size_ + 1, 0.0);
    posterior_.assign(size_ + 1, 0.0);
}

inline void SofosHistogram::AddCounts(double num_derived_copies, double num_ancestral_copies, double weight) {
    update_counts(num_derived_copies + alpha_, num_ancestral_copies + beta_, weight, &posterior_);
    update_bins(num_derived_copies, num_ancestral_copies, weight, &observed_);
}

inline void SofosHistogram::CalculatePrior() {
    prior_.assign(size_ + 1, 0.0);
    double weight = std::accumulate(posterior_.begin(), posterior_.end(), 0.0);
    update_counts(alpha_, beta_, weight, &prior_);
}

inline void SofosHistogram::Fold() {
    fold_histogram(&prior_);
    fold_histogram(&observed_);
    fold_histogram(&posterior_);
}

template <typename T>
void Sofos::AddCounts(const std::vector<T> &counts, int anc, double error_rate) {
    auto n_total = std::accumulate(counts.begin(), counts.end(), T{0});
    if(n_total == T{0}) {
        return;
    }
    // Number of ancestral copies and derived copies.
    auto n_anc = (anc < counts.size()) ? counts[anc] : T{0};
    auto n_der = n_total - n_anc;

    if(params_.flag_folded) {
        // When calculating the folded spectrum, we can't assume that either
        // allele is ancestral. Thus we will update the posterior weighted
        // by the possibility of each occurrence.
        // The probability that an der is ancestral is its frequency.
        error_rate = static_cast<double>(n_der) / n_total;
    }

    // Update Posterior and Observed using observed counts
    if(error_rate == 0.0) {
        histogram_.AddCounts(n_der, n_anc, 1.0);
    } else {
        // If there is some level of uncertainty about the ancestral allele,
        // weight the two possibilities.
        histogram_.AddCounts(n_der, n_anc, 1.0 - error_rate);
        histogram_.AddCounts(n_anc, n_der, error_rate);
    }
}

// calculate prior and fold if necessary
// should be done after the all data has been processed
inline void Sofos::FinishHistogram() {
    // calculate prior assuming no data
    histogram_.CalculatePrior();

    // fold histogram if necessary
    if(params_.flag_folded) {
        histogram_.Fold();
    }
}

// folds the histogram so that the second half is added to the first half.
// handles both odd and even vector sizes.
inline void fold_histogram(std::vector<double> *counts) {
    assert(counts != nullptr);
    for(int k = 0; k < counts->size() / 2; ++k) {
        (*counts)[k] += (*counts)[counts->size() - k - 1];
    }
    counts->resize((counts->size() + 1) / 2);
}

// calculate allele counts from AC+AN or GT values.
bool calculate_ac(bcf1_t *record, const bcf_hdr_t *header, int ac_which, std::vector<int> *ac_buffer);

// a functor to calculate allele counts from GP values
struct calculate_af_t {
    explicit calculate_af_t(int nsamples) {
        gp_buffer = make_buffer<float>(3 * nsamples);
        gt_buffer = make_buffer<int32_t>(2 * nsamples);
    }
    bool operator()(bcf1_t *record, const bcf_hdr_t *header, bool phred_scaled, std::vector<double> *af_buffer);

    // reusable buffers
    buffer_t<float> gp_buffer;
    buffer_t<int32_t> gt_buffer;
};

struct calculate_aa_t {
    calculate_aa_t() : char_buffer_(make_buffer<char>(64)) {}

    bool operator()(bcf1_t *record, const bcf_hdr_t *header, int *anc_allele);

    buffer_t<char> char_buffer_;
};

struct calculate_aaq_t {
    calculate_aaq_t() : int_buffer_(make_buffer<int32_t>(1)) {}

    bool operator()(bcf1_t *record, const bcf_hdr_t *header, double *error_rate);

    buffer_t<int32_t> int_buffer_;
};

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
    explicit Combinadic(int ploidy) : data_{0} { Reset(ploidy); }

    // Reset to the lowest k-combination
    void Reset(int ploidy) {
        data_ = 0;
        for(int j = 0; j < ploidy; ++j) {
            data_ = (data_ << 1) | 0x1u;
        }
    }

    // Permute the k-combination to the next highest value
    void Next() {
        // Gosper's Hack
        uint64_t u = data_ & -data_;
        uint64_t v = u + data_;
        data_ = v + (((v ^ data_) / u) >> 2);
    }

    // Access an integer representing a k-combination
    uint64_t Get() const { return data_; }

   protected:
    uint64_t data_;
};

inline double quality_to_p01(float x) { return -expm1(-x * (M_LN10 / 10.0)); }

inline double phred_to_p01(float x) { return exp(-x * (M_LN10 / 10.0)); }

inline double quality_to_p01(int x) {
    assert(x >= 0);
    static float data[96] = {
        0.000000000, 0.205671765, 0.369042656, 0.498812766, 0.601892829, 0.683772234, 0.748811357, 0.800473769,
        0.841510681, 0.874107459, 0.900000000, 0.920567177, 0.936904266, 0.949881277, 0.960189283, 0.968377223,
        0.974881136, 0.980047377, 0.984151068, 0.987410746, 0.990000000, 0.992056718, 0.993690427, 0.994988128,
        0.996018928, 0.996837722, 0.997488114, 0.998004738, 0.998415107, 0.998741075, 0.999000000, 0.999205672,
        0.999369043, 0.999498813, 0.999601893, 0.999683772, 0.999748811, 0.999800474, 0.999841511, 0.999874107,
        0.999900000, 0.999920567, 0.999936904, 0.999949881, 0.999960189, 0.999968377, 0.999974881, 0.999980047,
        0.999984151, 0.999987411, 0.999990000, 0.999992057, 0.999993690, 0.999994988, 0.999996019, 0.999996838,
        0.999997488, 0.999998005, 0.999998415, 0.999998741, 0.999999000, 0.999999206, 0.999999369, 0.999999499,
        0.999999602, 0.999999684, 0.999999749, 0.999999800, 0.999999842, 0.999999874, 0.999999900, 0.999999921,
        0.999999937, 0.999999950, 0.999999960, 0.999999968, 0.999999975, 0.999999980, 0.999999984, 0.999999987,
        0.999999990, 0.999999992, 0.999999994, 0.999999995, 0.999999996, 0.999999997, 0.999999997, 0.999999998,
        0.999999998, 0.999999999, 0.999999999, 0.999999999, 0.999999999, 0.999999999, 1.000000000, 1.000000000};
    return (x < 96) ? data[x] : 1.0;  // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
}

inline double phred_to_p01(int x) {
    assert(x >= 0);
    static float data[96] = {
        1.000000000, 0.794328235, 0.630957344, 0.501187234, 0.398107171, 0.316227766, 0.251188643, 0.199526231,
        0.158489319, 0.125892541, 0.100000000, 0.079432823, 0.063095734, 0.050118723, 0.039810717, 0.031622777,
        0.025118864, 0.019952623, 0.015848932, 0.012589254, 0.010000000, 0.007943282, 0.006309573, 0.005011872,
        0.003981072, 0.003162278, 0.002511886, 0.001995262, 0.001584893, 0.001258925, 0.001000000, 0.000794328,
        0.000630957, 0.000501187, 0.000398107, 0.000316228, 0.000251189, 0.000199526, 0.000158489, 0.000125893,
        0.000100000, 0.000079433, 0.000063096, 0.000050119, 0.000039811, 0.000031623, 0.000025119, 0.000019953,
        0.000015849, 0.000012589, 0.000010000, 0.000007943, 0.000006310, 0.000005012, 0.000003981, 0.000003162,
        0.000002512, 0.000001995, 0.000001585, 0.000001259, 0.000001000, 0.000000794, 0.000000631, 0.000000501,
        0.000000398, 0.000000316, 0.000000251, 0.000000200, 0.000000158, 0.000000126, 0.000000100, 0.000000079,
        0.000000063, 0.000000050, 0.000000040, 0.000000032, 0.000000025, 0.000000020, 0.000000016, 0.000000013,
        0.000000010, 0.000000008, 0.000000006, 0.000000005, 0.000000004, 0.000000003, 0.000000003, 0.000000002,
        0.000000002, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000001, 0.000000000, 0.000000000};
    return (x < 96) ? data[x] : 0.0;  // NOLINT(cppcoreguidelines-pro-bounds-constant-array-index)
}

#endif  // SOFOS_SOFOS_HPP
