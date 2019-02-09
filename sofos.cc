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

#include "sofos.hpp"
#include "vcf.hpp"

#include <unistd.h>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <numeric>
#include <string>
#include <vector>

bool g_sofos_quiet = false;

// Utility class to output a message every X number of seconds
class ProgressMessage {
   public:
    ProgressMessage(std::ostream &os, int step) : output_(os), step_{step}, next_{step} {
        start_ = std::chrono::steady_clock::now();
    }

    template <typename call_back_t>
    void operator()(call_back_t call_back, bool force = false);

   protected:
    std::ostream &output_;
    int step_;
    int next_;

    decltype(std::chrono::steady_clock::now()) start_;
};

template <typename call_back_t>
void ProgressMessage::operator()(call_back_t call_back, bool force) {
    auto now = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = now - start_;
    if(elapsed.count() >= next_ || force) {
        next_ += step_;
        output_ << "## [" << std::fixed << std::setprecision(4) << elapsed.count() << "s elapsed] ";
        output_ << std::scientific << std::setprecision(6);
        call_back(output_);
        output_ << std::endl;
    }
}

// Process a VCF file and add its contents to the histogram
void Sofos::RescaleBcf(const char *path) {
    assert(path != nullptr);

    // open the input file or throw on error
    BcfReader reader{path};

    // Setup for reading from vcf file
    std::vector<int> ac_buffer;
    std::vector<double> af_buffer;

    size_t nsamples = bcf_hdr_nsamples(reader.header());
    if(nsamples <= 0) {
        throw std::invalid_argument("input variant file has no samples.");
    }
    // resize the gp and gt buffers
    calculate_af_t calculate_af(nsamples);
    calculate_aa_t calculate_aa{};
    calculate_aaq_t calculate_aaq{};

    // Setup for reading AC tags
    int ac_which = BCF_UN_FMT;
    // if AN and AC tags are defined, use them
    if(bcf_hdr_id2int(reader.header(), BCF_DT_ID, "AN") >= 0 && bcf_hdr_id2int(reader.header(), BCF_DT_ID, "AC") >= 0) {
        ac_which |= BCF_UN_INFO;
    }

    // begin timer that activates every 60s
    ProgressMessage progress{std::cerr, 60};

    // begin reading from the input file
    size_t nsites = 0;
    reader([&](bcf1_t *record, const bcf_hdr_t *header) {
        assert(record != nullptr);
        assert(header != nullptr);

        // identify which allele is ancestral and its error rate
        int anc = 0;
        double error_rate = params_.error_rate;
        if(params_.flag_refalt == false) {
            if(calculate_aa(record, header, &anc) == false) {
                // skip line on failure
                return;
            }
            if(error_rate == 0.0) {
                // check for AAQ tag
                calculate_aaq(record, header, &error_rate);
            }
        }
        if(!params_.flag_use_gp) {
            // Calculate the counts of each allele from the AC/AN tags or directly from GT.
            if(!calculate_ac(record, header, ac_which, &ac_buffer)) {
                // Skip line if it fails
                return;
            }
            // Add this site to the histogram
            AddCounts(ac_buffer, anc, error_rate);
        } else {
            // Calculate the counts of each allele from the GP tags
            if(!calculate_af(record, header, params_.flag_phred_gp, &af_buffer)) {
                // Skip line if it fails
                return;
            }
            // Add this site to the histogram
            AddCounts(af_buffer, anc, error_rate);
        }

        // Increment the number of observed sites.
        nsites += 1;

        // Output a progress message every minute or so.
        if(!g_sofos_quiet) {
            progress([=](std::ostream &os) {
                os << nsites << " sites processed --- at " << bcf_seqname(header, record) << ":" << record->pos + 1;
            });
        }
    });

    // Output a progress message at end.
    if(!g_sofos_quiet) {
        progress([=](std::ostream &os) { os << nsites << " sites processed"; }, true);
    }

    // Add zero-count sites (for controlling ascertainment bias)
    if(params_.zero_count > 0.0) {
        histogram_.AddCounts(0.0, params_.ploidy * nsamples, params_.zero_count);
    }
}

// Determine the ancestral allele
bool calculate_aa_t::operator()(bcf1_t *record, const bcf_hdr_t *header, int *anc_allele) {
    assert(record != nullptr);
    assert(header != nullptr);
    assert(anc_allele != nullptr);

    int n = get_info_string(header, record, "AA", &char_buffer_);
    if(n <= 0) {
        // missing, unknown, or empty tag
        return false;
    }
    // Compare the AA tag to the refalt allele to find the index of the ancestral allele
    // If this fails, anc = n_allele
    bcf_unpack(record, BCF_UN_STR);
    int anc;
    for(anc = 0; anc < record->n_allele; ++anc) {
        if(strcmp(char_buffer_.data.get(), record->d.allele[anc]) == 0) {
            break;
        }
    }
    *anc_allele = anc;
    return true;
}

// Determine ancestral error rate
bool calculate_aaq_t::operator()(bcf1_t *record, const bcf_hdr_t *header, double *error_rate) {
    assert(record != nullptr);
    assert(header != nullptr);
    assert(error_rate != nullptr);

    // If error_rate is missing, try to look for an AAQ tag
    int n = get_info_int32(header, record, "AAQ", &int_buffer_);
    if(n <= 0) {
        return false;
    }
    *error_rate = phred_to_p01(int_buffer_.data[0]);
    return true;
}

// Calculate allele counts from AC/AN or GT values
bool calculate_ac(bcf1_t *record, const bcf_hdr_t *header, int ac_which, std::vector<int> *ac_buffer) {
    assert(record != nullptr);
    assert(header != nullptr);
    assert(ac_buffer != nullptr);

    ac_buffer->assign(record->n_allele, 0);
    if(bcf_calc_ac(header, record, ac_buffer->data(), ac_which) <= 0) {
        return false;
    }
    return true;
}

// Calculate allele counts from GP values
bool calculate_af_t::operator()(bcf1_t *record, const bcf_hdr_t *header, bool phred_scaled,
                                std::vector<double> *af_buffer) {
    assert(record != nullptr);
    assert(header != nullptr);
    assert(af_buffer != nullptr);

    size_t nsamples = bcf_hdr_nsamples(header);
    assert(nsamples > 0);

    // Calculate counts from GP
    int n = get_format_float(header, record, "GP", &gp_buffer);
    if(n <= 0) {
        // missing, unknown, or empty tag
        return false;
    }
    int gp_width = n / nsamples;

    // convert GP buffer if needed
    if(phred_scaled) {
        for(int i = 0; i < n; ++i) {
            gp_buffer.data[i] = phred_to_p01(gp_buffer.data[i]);
        }
    }

    // Calculate GT (for ploidy)
    n = get_genotypes(header, record, &gt_buffer);
    if(n <= 0) {
        // missing, unknown, or empty tag
        return false;
    }
    int gt_width = n / nsamples;

    // Allocate vector to estimate allele counts from GP
    af_buffer->assign(record->n_allele, 0.0);
    for(int i = 0; i < nsamples; i++) {
        // identify ploidy for this sample
        int *p_gt = gt_buffer.data.get() + i * gt_width;
        int sample_ploidy = 0;
        for(; sample_ploidy < gt_width; ++sample_ploidy) {
            if(p_gt[sample_ploidy] == bcf_int32_vector_end) {
                break;
            }
        }
        if(sample_ploidy == 0) {
            continue;
        }
        // Setup a sequence of k-permutations
        Combinadic gp_alleles{sample_ploidy};
        int gp_alleles_sz = record->n_allele + sample_ploidy - 1;

        // Step thorough each element of the gp
        float *p_gp = gp_buffer.data.get() + i * gp_width;
        for(int j = 0; j < gp_width; ++j, gp_alleles.Next()) {
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
            for(int h = 0; h < gp_alleles_sz; ++h) {
                // if the lowest bit is 1, we have an allele
                // if the lowest bit is 0, we move on to the next allele
                if((alleles & 0x1) == 1) {
                    (*af_buffer)[pos] += p_gp[j];
                } else {
                    pos += 1;
                }
                // shit to the next bit
                alleles = alleles >> 1;
            }
        }
    }
    return true;
}

// Generate a timestamp. Returns an ISO formatted timestring
// and seconds since epoch.
std::pair<std::string, std::string> timestamp() {
    char buffer[32];
    auto now = std::chrono::system_clock::now();
    auto now_t = std::chrono::system_clock::to_time_t(now);
    std::strftime(&buffer[0], 32, "%FT%T%z", std::localtime(&now_t));
    auto epoch = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch());
    return {buffer, std::to_string(epoch.count())};
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
    double ab = a + b;
    int size = counts->size() - 1;
    double scale = lgamma(a) + lgamma(b) + lgamma(ab + size) - lgamma(ab) - lgamma(size + 1);
    for(int k = 0; k <= size; ++k) {
        // rescale assuming beta-binomial model
        double d = lgamma(a + k) + lgamma(b + size - k);
        d -= lgamma(k + 1) + lgamma(size - k + 1);
        d -= scale;
        (*counts)[k] += weight * exp(d);
    }
}

// maps a/(a+b) to a bin in the range of [0,N]/N
// rounds to nearest bin. Thus the bins on the edge may
// have lower counts than those in the middle.
void update_bins(double a, double b, double weight, std::vector<double> *bins) {
    assert(bins != nullptr);
    assert(a >= 0.0 && b >= 0.0);
    if(a + b == 0.0) {
        return;
    }
    // bins are from [0,n]
    int n = bins->size() - 1;
    double f = a / (a + b) * n;
    int k = lround(f);
    if(k > n) {
        k = n;
    }
    (*bins)[k] += weight;
}

// Output header including program runtime information and options
void output_header(std::ostream &os, const Sofos &sofos, const std::vector<const char *> &paths) {
    auto stamp = timestamp();
    os << "#SoFoS v2.0\n";
    os << "#date=" << stamp.first << "\n";
    os << "#epoch=" << stamp.second << "\n";
    for(auto &&path : paths) {
        os << "#path=" << path << "\n";
    }

    os << std::setprecision(std::numeric_limits<double>::max_digits10);
    os << "#alpha=" << sofos.params_.alpha << "\n";
    os << "#beta=" << sofos.params_.beta << "\n";
    os << "#size=" << sofos.params_.size << "\n";
    os << "#folded=" << (sofos.params_.flag_folded ? 1 : 0) << "\n";
    os << "#refalt=" << (sofos.params_.flag_refalt ? 1 : 0) << "\n";

    if(sofos.params_.flag_use_gp) {
        os << "#use_gp=" << (sofos.params_.flag_phred_gp ? 2 : 1) << "\n";
    }

    if(sofos.params_.error_rate > 0.0) {
        os << "#error_rate=" << sofos.params_.error_rate << "\n";
    }
    if(sofos.params_.zero_count > 0.0) {
        os << "#zero=" << sofos.params_.zero_count << "\n";
        os << "#ploidy=" << sofos.params_.ploidy << "\n";
    }

    // spacer
    os << "#\n";
}

// Output body
void output_body(std::ostream &os, const Sofos &sofos) {
    // output resulting scale
    os << std::setprecision(std::numeric_limits<double>::max_digits10);
    os << "Number,Prior,Observed,Posterior\n";
    for(int i = 0; i < sofos.histogram().num_rows(); ++i) {
        auto &&row = sofos.histogram().row(i);
        os << i << "," << row[0] << "," << row[1] << "," << row[2] << "\n";
    }
    os << std::flush;
}
