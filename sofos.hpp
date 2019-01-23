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

#include <vector>
#include <numeric>

// Globals
extern bool g_sofos_quiet;

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

inline
SofosHistogram::SofosHistogram(int size, double alpha, double beta) :
    size_{size}, alpha_{alpha}, beta_{beta}
{
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

#endif // SOFOS_SOFOS_HPP
