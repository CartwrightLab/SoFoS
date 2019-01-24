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
#pragma once
#ifndef SOFOS_SOFOS_HPP
#define SOFOS_SOFOS_HPP

#include "vcf.hpp"

#include <vector>
#include <numeric>
#include <cassert>
#include <ostream>

// Globals
extern bool g_sofos_quiet;

// Flags used to control sofos_main
constexpr unsigned int SOFOS_FLAG_DEFAULT=0;
constexpr unsigned int SOFOS_FLAG_FOLDED=1;
constexpr unsigned int SOFOS_FLAG_REFALT=2;
constexpr unsigned int SOFOS_FLAG_USE_GP=4;
constexpr unsigned int SOFOS_FLAG_PHRED_GP=8;

struct sofos_params_t {
    double alpha{1.0};
    double beta{1.0};
    int    size{9};
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

    const std::vector<double>& prior() const { return prior_; }
    const std::vector<double>& observed() const { return observed_; }
    const std::vector<double>& posterior() const { return posterior_; }

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
    Sofos(const sofos_params_t& params) : params_{params},
        histogram_{params.size, params.alpha, params.beta}
    {
        assert(0.0 <= params.error_rate && params.error_rate <= 1.0);
        assert(params.zero_count >= 0.0);
        assert(params.ploidy >= 0.0);

        char_buffer_ = make_buffer<char>(64);
        int_buffer_ = make_buffer<int32_t>(1);
    }

    void RescaleFile(const char *path);

    void ResetHistogram() {
        histogram_.ZeroCounts();
    }

    void FinishHistogram();

    const SofosHistogram histogram() const { return histogram_; }

protected:
    std::pair<int,double> GetAncestor(bcf1_t *record, const bcf_hdr_t *header);

    // This is a template in case SofosHistogram::AddCounts becomes an overloaded function
    template<typename T>
    void AddCounts(const std::vector<T> &counts, int anc, double error_rate);

private:
    sofos_params_t params_;

    SofosHistogram histogram_;

    buffer_t<char> char_buffer_;
    buffer_t<int32_t> int_buffer_;

public:
    friend void output_header(std::ostream &os, const Sofos &sofos,
        const std::vector<const char *> &paths);
    friend void output_body(std::ostream &os, const Sofos &sofos);
};

inline
SofosHistogram::SofosHistogram(int size, double alpha, double beta) :
    size_{size}, alpha_{alpha}, beta_{beta}
{
    assert(size > 0);
    assert(alpha > 0.0 && beta > 0.0);

    ZeroCounts();
}

inline
void SofosHistogram::ZeroCounts() {
    // A sample of N copies can have [0,N] derived alleles.
    prior_.assign(size_+1, 0.0);
    observed_.assign(size_+1,0.0);
    posterior_.assign(size_+1, 0.0);
}

inline
void SofosHistogram::AddCounts(double num_derived_copies, double num_ancestral_copies, double weight) {
    update_counts(num_derived_copies+alpha_, num_ancestral_copies+beta_, weight, &posterior_);
    update_bins(num_derived_copies, num_ancestral_copies, weight, &observed_);
}

inline
void SofosHistogram::CalculatePrior() {
    prior_.assign(size_+1, 0.0);
    double weight = std::accumulate(posterior_.begin(),posterior_.end(),0.0);
    update_counts(alpha_, beta_, weight, &prior_);
}

inline
void SofosHistogram::Fold() {
    fold_histogram(&prior_);
    fold_histogram(&observed_);
    fold_histogram(&posterior_);
}

template<typename T>
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
        error_rate = static_cast<double>(n_der)/n_total;
    }

    // Update Posterior and Observed using observed counts
    if(error_rate == 0.0) {
        histogram_.AddCounts(n_der, n_anc, 1.0);
    } else {
        // If there is some level of uncertainty about the ancestral allele,
        // weight the two possibilities.
        histogram_.AddCounts(n_der, n_anc, 1.0-error_rate);
        histogram_.AddCounts(n_anc, n_der, error_rate);
    }
}

#endif // SOFOS_SOFOS_HPP