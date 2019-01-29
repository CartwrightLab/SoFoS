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
#include "catch.hpp"
#include "test_utils.hpp"

#include <sstream>
#endif

#include "sofos.hpp"
#include "vcf.hpp"

#include <unistd.h>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <numeric>
#include <string>
#include <vector>

#ifndef SOFOS_UNIT_TESTS
bool g_sofos_quiet = false;
#else
bool g_sofos_quiet = true;

struct kstring_obj_t : public kstring_t {
    kstring_obj_t() : kstring_t{0, 0, nullptr} {}

    ~kstring_obj_t() {
        if(s != nullptr) {
            free(s);
        }
    }
};
#endif

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
    calculate_aa_t() : char_buffer_{make_buffer<char>(64)} {}

    bool operator()(bcf1_t *record, const bcf_hdr_t *header, int *anc_allele);

    buffer_t<char> char_buffer_;
};

struct calculate_aaq_t {
    calculate_aaq_t() : int_buffer_{make_buffer<int32_t>(1)} {}

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

// Utility class to output a message every X number of seconds
class ProgressMessage {
   public:
    ProgressMessage(std::ostream &os, int step) : output_{os}, step_{step}, next_{step} {
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

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("Sofos::RescaleBcf rescales SFS to a new sample size") {
    using namespace Catch::literals;

    char header_str[] =
        "##fileformat=VCFv4.2\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n"
        "##INFO=<ID=AAQ,Number=1,Type=Integer,Description=\"AA Quality\">\n"
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">\n"
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype Posteriors\">\n"
        "##contig=<ID=1,length=10000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n";

    auto rescale = [&](const sofos_params_t &params, const char *lines) -> SofosHistogram {
        std::string data = "data:";
        {
            // detect what version of the data string we can use
            auto p = hts_open("data:,##fileformat=VCFv4.2", "r");
            if(p != nullptr) {
                if(p->format.format == vcf) {
                    data += ",";
                }
                hts_close(p);
            }
        }
        data += header_str;
        data += lines;
        Sofos sofos{params};
        sofos.RescaleBcf(data.c_str());
        sofos.FinishHistogram();
        return sofos.histogram();
    };

    sofos_params_t params;
    params.size = 2;

    SECTION("when refalt=false") {
        params.flag_refalt = false;
        const char str[] = "1\t1\t.\tA\tC\t.\t.\tAA=C;AN=10;AC=9\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.3333333_a, 0.3333333_a, 0.3333333_a};
        std::vector<Approx> expected_observed = {1.0_a, 0.0_a, 0.0_a};
        std::vector<Approx> expected_posterior = {0.70512821_a, 0.25641026_a, 0.03846154_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
    SECTION("when refalt=false AA is missing") {
        params.flag_refalt = false;
        const char str[] = "1\t1\t.\tA\tC\t.\t.\tAN=10;AC=9\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.0_a, 0.0_a, 0.0_a};
        std::vector<Approx> expected_observed = {0.0_a, 0.0_a, 0.0_a};
        std::vector<Approx> expected_posterior = {0.0_a, 0.0_a, 0.0_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
    SECTION("when refalt=true") {
        params.flag_refalt = true;
        const char str[] = "1\t1\t.\tA\tC\t.\t.\tAA=G;AN=10;AC=1\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.3333333_a, 0.3333333_a, 0.3333333_a};
        std::vector<Approx> expected_observed = {1.0_a, 0.0_a, 0.0_a};
        std::vector<Approx> expected_posterior = {0.70512821_a, 0.25641026_a, 0.03846154_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
    SECTION("when AC and GT are missing") {
        params.flag_refalt = true;
        const char str[] = "1\t1\t.\tA\tC\t.\t.\t.\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.0_a, 0.0_a, 0.0_a};
        std::vector<Approx> expected_observed = {0.0_a, 0.0_a, 0.0_a};
        std::vector<Approx> expected_posterior = {0.0_a, 0.0_a, 0.0_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
    SECTION("when there are two sites") {
        params.flag_refalt = true;
        const char str[] =
            "1\t1\t.\tA\tC\t.\t.\tAN=10;AC=1\n"
            "1\t2\t.\tA\tC\t.\t.\tAN=10;AC=5\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.6666667_a, 0.6666667_a, 0.6666667_a};
        std::vector<Approx> expected_observed = {1.0_a, 1.0_a, 0.0_a};
        std::vector<Approx> expected_posterior = {0.9743590_a, 0.7179487_a, 0.3076923_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
    SECTION("when use_gp=true") {
        params.flag_use_gp = true;
        const char str[] = "1\t1\t.\tA\tC\t.\t.\tAA=C;AN=10;AC=0\tGT:GP\t0/0:1,0,0\t0/1:0,1,0\t1/1:0,0,1\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.3333333_a, 0.3333333_a, 0.3333333_a};
        std::vector<Approx> expected_observed = {0.0_a, 1.0_a, 0.0_a};
        std::vector<Approx> expected_posterior = {0.2777778_a, 0.4444444_a, 0.2777778_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
    SECTION("when use_gp=true and GP is missing") {
        params.flag_use_gp = true;
        const char str[] = "1\t1\t.\tA\tC\t.\t.\tAA=C;AN=10;AC=0\tGT\t0/0\t0/1\t1/1\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.0_a, 0.0_a, 0.0_a};
        std::vector<Approx> expected_observed = {0.0_a, 0.0_a, 0.0_a};
        std::vector<Approx> expected_posterior = {0.0_a, 0.0_a, 0.0_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
    SECTION("when use_gp=true and GT is missing") {
        params.flag_use_gp = true;
        const char str[] = "1\t1\t.\tA\tC\t.\t.\tAA=C;AN=10;AC=0\tGP\t1,0,0\t0,1,0\t0,0,1\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.0_a, 0.0_a, 0.0_a};
        std::vector<Approx> expected_observed = {0.0_a, 0.0_a, 0.0_a};
        std::vector<Approx> expected_posterior = {0.0_a, 0.0_a, 0.0_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
    SECTION("when there is no data") {
        params.flag_use_gp = true;
        const char str[] = "1\t1\t.\tA\tC\t.\t.\tAA=C;AN=0;AC=0\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.0_a, 0.0_a, 0.0_a};
        std::vector<Approx> expected_observed = {0.0_a, 0.0_a, 0.0_a};
        std::vector<Approx> expected_posterior = {0.0_a, 0.0_a, 0.0_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
    SECTION("when use_gp=true and gp_phred=true") {
        params.flag_use_gp = true;
        params.flag_phred_gp = true;
        const char str[] = "1\t1\t.\tA\tC\t.\t.\tAA=C;AN=10;AC=0\tGT:GP\t0/0:0,100,100\t0/1:100,0,100\t1/1:100,100,0\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.3333333_a, 0.3333333_a, 0.3333333_a};
        std::vector<Approx> expected_observed = {0.0_a, 1.0_a, 0.0_a};
        std::vector<Approx> expected_posterior = {0.2777778_a, 0.4444444_a, 0.2777778_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
    SECTION("when zero_count is specified") {
        params.flag_refalt = true;
        params.zero_count = 1.0;
        params.ploidy = 10.0 / 3.0;
        const char str[] = "1\t1\t.\tA\tC\t.\t.\tAA=G;AN=10;AC=1\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.6666667_a, 0.6666667_a, 0.6666667_a};
        std::vector<Approx> expected_observed = {2.0_a, 0.0_a, 0.0_a};
        std::vector<Approx> expected_posterior = {1.55128205_a, 0.39743590_a, 0.05128205_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
    SECTION("when error_rate is specified") {
        params.flag_refalt = false;
        params.error_rate = 0.5;
        const char str[] = "1\t1\t.\tA\tC\t.\t.\tAA=A;AN=10;AC=1\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.3333333_a, 0.3333333_a, 0.3333333_a};
        std::vector<Approx> expected_observed = {0.5_a, 0.0_a, 0.5_a};
        std::vector<Approx> expected_posterior = {0.3717949_a, 0.2564103_a, 0.3717949_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
    SECTION("when AAQ is specified") {
        params.flag_refalt = false;
        const char str[] = "1\t1\t.\tA\tC\t.\t.\tAA=A;AN=10;AC=1;AAQ=3\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.3333333_a, 0.3333333_a, 0.3333333_a};
        std::vector<Approx> expected_observed = {0.4988128_a, 0.0_a, 0.5011872_a};
        std::vector<Approx> expected_posterior = {0.3710034_a, 0.2564103_a, 0.3725864_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
    SECTION("when AAQ and error_rate are specified") {
        params.flag_refalt = false;
        params.error_rate = 0.5;
        const char str[] = "1\t1\t.\tA\tC\t.\t.\tAA=A;AN=10;AC=1;AAQ=10\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.3333333_a, 0.3333333_a, 0.3333333_a};
        std::vector<Approx> expected_observed = {0.5_a, 0.0_a, 0.5_a};
        std::vector<Approx> expected_posterior = {0.3717949_a, 0.2564103_a, 0.3717949_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
    SECTION("when folded=true") {
        params.flag_refalt = true;
        params.flag_folded = true;
        const char str[] = "1\t1\t.\tA\tC\t.\t.\tAA=G;AN=10;AC=1\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.6666667_a, 0.3333333_a};
        std::vector<Approx> expected_observed = {1.0_a, 0.0_a};
        std::vector<Approx> expected_posterior = {0.7435897_a, 0.25641026_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
}
#endif

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

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("calculate_aa_t calculates the ancestral allele id") {
    using namespace Catch::literals;

    char header_str[] =
        "##fileformat=VCFv4.2\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n"
        "##INFO=<ID=AAQ,Number=1,Type=Integer,Description=\"AA Quality\">\n"
        "##contig=<ID=1,length=10000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n";

    std::unique_ptr<bcf1_t, detail::bcf_free_t> record_{bcf_init()};
    REQUIRE((record_));

    std::unique_ptr<bcf_hdr_t, detail::header_free_t> header_{bcf_hdr_init("w")};
    REQUIRE((header_));

    auto header = header_.get();
    auto record = record_.get();

    int hret = bcf_hdr_parse(header, header_str);
    REQUIRE(hret == 0);

    auto parse_aa = [=](const char *line, int *anc) -> bool {
        kstring_obj_t kstr;
        if(kputs(line, &kstr) < 0) {
            return false;
        }
        int ret = vcf_parse(&kstr, header, record);
        if(ret != 0) {
            return false;
        }
        calculate_aa_t calculate_aa;

        return calculate_aa(record, header, anc);
    };
    int anc;
    bool ret;

    ret = parse_aa("1\t1\t.\tA\tC\t.\t.\tAA=A", &anc);
    REQUIRE(ret == true);
    CHECK(anc == 0);

    ret = parse_aa("1\t1\t.\tA\tC\t.\t.\tAA=C", &anc);
    REQUIRE(ret == true);
    CHECK(anc == 1);

    ret = parse_aa("1\t1\t.\tA\tC\t.\t.\tAA=G", &anc);
    REQUIRE(ret == true);
    CHECK(anc == 2);

    ret = parse_aa("1\t1\t.\tA\tC\t.\t.\t.", &anc);
    REQUIRE(ret == false);
}
#endif

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

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("calculate_aaq_t calculates the ancestral allele error_rate") {
    using namespace Catch::literals;

    char header_str[] =
        "##fileformat=VCFv4.2\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n"
        "##INFO=<ID=AAQ,Number=1,Type=Integer,Description=\"AA Quality\">\n"
        "##contig=<ID=1,length=10000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n";

    std::unique_ptr<bcf1_t, detail::bcf_free_t> record_{bcf_init()};
    REQUIRE((record_));

    std::unique_ptr<bcf_hdr_t, detail::header_free_t> header_{bcf_hdr_init("w")};
    REQUIRE((header_));

    auto header = header_.get();
    auto record = record_.get();

    int hret = bcf_hdr_parse(header, header_str);
    REQUIRE(hret == 0);

    auto parse_aaq = [=](const char *line, double *err) -> bool {
        kstring_obj_t kstr;
        if(kputs(line, &kstr) < 0) {
            return false;
        }
        int ret = vcf_parse(&kstr, header, record);
        if(ret != 0) {
            return false;
        }
        calculate_aaq_t calculate_aaq;

        return calculate_aaq(record, header, err);
    };
    double error_rate;
    bool ret;

    ret = parse_aaq("1\t1\t.\tA\tC\t.\t.\tAAQ=10", &error_rate);
    REQUIRE(ret == true);
    CHECK(error_rate == 0.1_a);

    ret = parse_aaq("1\t1\t.\tA\tC\t.\t.\tAAQ=100", &error_rate);
    REQUIRE(ret == true);
    CHECK(error_rate == 0.0);

    ret = parse_aaq("1\t1\t.\tA\tC\t.\t.\t.", &error_rate);
    REQUIRE(ret == false);
}
#endif

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

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("calculate_ac() calculates allele counts from GT or AC/AN values") {
    char header_str[] =
        "##fileformat=VCFv4.2\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">\n"
        "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##contig=<ID=1,length=10000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n";

    std::unique_ptr<bcf1_t, detail::bcf_free_t> record_{bcf_init()};
    REQUIRE((record_));

    std::unique_ptr<bcf_hdr_t, detail::header_free_t> header_{bcf_hdr_init("w")};
    REQUIRE((header_));

    auto header = header_.get();
    auto record = record_.get();

    int ret = bcf_hdr_parse(header, header_str);
    REQUIRE(ret == 0);

    auto parse_ac = [=](const char *line, int ac_which, int sz) -> std::vector<int> {
        kstring_obj_t kstr;
        if(kputs(line, &kstr) < 0) {
            return {};
        }
        int ret = vcf_parse(&kstr, header, record);
        if(ret != 0) {
            return {};
        }
        std::vector<int> ac(sz, 0);
        if(!calculate_ac(record, header, ac_which, &ac)) {
            return {};
        }
        return ac;
    };

    SECTION("when it can only read from info") {
        auto ac = parse_ac("1\t1\t.\tA\tC\t.\t.\tAC=4;AN=10\tGT\t0/0\t0/0\t1/1", BCF_UN_INFO, 2);
        std::vector<int> expected = {6, 4};
        CHECK(ac == expected);
    }
    SECTION("when it can only read from format") {
        auto ac = parse_ac("1\t1\t.\tA\tC\t.\t.\tAC=4;AN=10\tGT\t0/0\t0/0\t1/1", BCF_UN_FMT, 2);
        std::vector<int> expected = {4, 2};
        CHECK(ac == expected);
    }
    SECTION("when samples are missing") {
        auto ac = parse_ac("1\t1\t.\tA\tC\t.\t.\tAC=4;AN=10", BCF_UN_FMT | BCF_UN_INFO, 2);
        std::vector<int> expected = {6, 4};
        CHECK(ac == expected);
    }
    SECTION("when info is missing") {
        auto ac = parse_ac("1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/0\t0/0\t1/1", BCF_UN_FMT | BCF_UN_INFO, 2);
        std::vector<int> expected = {4, 2};
        CHECK(ac == expected);
    }
    SECTION("when ploidy is mixed") {
        auto ac = parse_ac("1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/0\t0/0\t1", BCF_UN_FMT | BCF_UN_INFO, 2);
        std::vector<int> expected = {4, 1};
        CHECK(ac == expected);
    }
    SECTION("when there is missing data") {
        auto ac = parse_ac("1\t1\t.\tA\tC\t.\t.\t.\tGT\t.\t0/.\t1", BCF_UN_FMT | BCF_UN_INFO, 2);
        std::vector<int> expected = {1, 1};
        CHECK(ac == expected);
    }
    SECTION("when there are three alleles") {
        auto ac = parse_ac("1\t1\t.\tA\tC,AA\t.\t.\tAC=4,6;AN=10", BCF_UN_FMT | BCF_UN_INFO, 3);
        std::vector<int> expected = {0, 4, 6};
        CHECK(ac == expected);
    }
    SECTION("when there is no data") {
        auto ac = parse_ac("1\t1\t.\tA\tC\t.\t.\t.", BCF_UN_FMT | BCF_UN_INFO, 2);
        std::vector<int> expected = {};
        CHECK(ac == expected);
    }
}
#endif

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

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("calculate_af() calculates allele counts from GP values") {
    using namespace Catch::literals;
    char header_str[] =
        "##fileformat=VCFv4.2\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
        "##FORMAT=<ID=GP,Number=G,Type=Float,Description=\"Genotype Posterior\">\n"
        "##contig=<ID=1,length=10000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n";

    std::unique_ptr<bcf1_t, detail::bcf_free_t> record_{bcf_init()};
    REQUIRE((record_));

    std::unique_ptr<bcf_hdr_t, detail::header_free_t> header_{bcf_hdr_init("w")};
    REQUIRE((header_));

    auto header = header_.get();
    auto record = record_.get();

    int ret = bcf_hdr_parse(header, header_str);
    REQUIRE(ret == 0);

    calculate_af_t calculate_af{3};

    auto parse_gp = [=, &calculate_af](const char *line, bool phread, int sz) -> std::vector<double> {
        kstring_obj_t kstr;
        if(kputs(line, &kstr) < 0) {
            return {};
        }
        int ret = vcf_parse(&kstr, header, record);
        if(ret != 0) {
            return {};
        }
        std::vector<double> af(sz, 0.0);
        if(!calculate_af(record, header, phread, &af)) {
            return {};
        }
        return af;
    };

    SECTION("when GP is standard scaled") {
        auto ac =
            parse_gp("1\t1\t.\tA\tC\t.\t.\t.\tGT:GP\t0/0:0.9,0.1,0.0\t0/0:0.8,0.1,0.1\t0/0:0.7,0.1,0.2", false, 2);
        std::vector<Approx> expected = {5.1_a, 0.9_a};
        CHECK(ac == expected);
    }

    SECTION("when GP is phread scaled") {
        auto ac = parse_gp(
            "1\t1\t.\tA\tC\t.\t.\t.\tGT:GP\t0/0:0.4575749,10.0,64.0\t0/0:0.9691001,10.0,10.0\t0/"
            "0:1.5490196,10.0,6.9897000",
            true, 2);
        std::vector<Approx> expected = {5.1_a, 0.9_a};
        CHECK(ac == expected);
    }

    SECTION("when GP has missing values") {
        auto ac = parse_gp("1\t1\t.\tA\tC\t.\t.\t.\tGT:GP\t0/0:1.0\t0/0:.\t0/0:.,.,1.0", false, 2);
        std::vector<Approx> expected = {2.0_a, 2.0_a};
        CHECK(ac == expected);
    }

    SECTION("when GP has mixed ploidy") {
        auto ac =
            parse_gp("1\t1\t.\tA\tC\t.\t.\t.\tGT:GP\t0/0/0:0.4,0.3,0.2,0.1\t0/0:0.8,0.1,0.1\t0:0.9,0.1", false, 2);
        std::vector<Approx> expected = {4.6_a, 1.4_a};
        CHECK(ac == expected);
    }

    SECTION("when GP is missing") {
        auto ac = parse_gp("1\t1\t.\tA\tC\t.\t.\t.\tGT\t0/0\t0/0\t0/0", false, 2);
        std::vector<Approx> expected = {};
        CHECK(ac == expected);
    }

    SECTION("when GT is missing") {
        auto ac = parse_gp("1\t1\t.\tA\tC\t.\t.\t.\tGP\t1,0,0\t1,0,0\t1,0,0", false, 2);
        std::vector<Approx> expected = {};
        CHECK(ac == expected);
    }

    SECTION("when GT has missing values") {
        auto ac = parse_gp("1\t1\t.\tA\tC\t.\t.\t.\tGT:GP\t./.:0.9,0.1,0.0\t./.:0.8,0.1,0.1\t.:0.9,0.1", false, 2);
        std::vector<Approx> expected = {4.5_a, 0.5_a};
        CHECK(ac == expected);
    }
}
#endif

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

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("timestamp() returns the current localtime and epoch") {
    auto now = std::chrono::system_clock::now();
    auto now_time = std::chrono::system_clock::to_time_t(now);
    auto epoch = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch());

    auto stamp = timestamp();
    auto stamp_epoch = std::stoull(stamp.second);

    std::tm stamp_tm = {};
    std::istringstream ss(stamp.first);
    ss >> std::get_time(&stamp_tm, "%Y-%m-%dT%H:%M:%S");
    REQUIRE(ss.good());

    auto stamp_time = std::mktime(&stamp_tm);

    REQUIRE(stamp_time != -1);
    double time_diff = fabs(std::difftime(stamp_time, now_time));
    CHECK(time_diff <= 1.0);

    CHECK(stamp_epoch <= epoch.count() + 1000);
}
#endif

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

/*
#Rcode for generating expectations
library(rmutil)
dbetabinom(0:n,n,a/(a+b),a+b)
*/
#ifdef SOFOS_UNIT_TESTS
TEST_CASE("update_counts adds resampled data to a vector") {
    using namespace Catch::literals;
    auto zero = Approx(0.0).margin(std::numeric_limits<float>::epsilon() * 100);
    std::vector<double> expected = {0.714285714, 0.219780220, 0.054945055, 0.009990010, 0.000999001};
    std::vector<Approx> expected_a;

    SECTION("when counts are zero") {
        std::vector<double> counts = {0, 0, 0, 0, 0};
        update_counts(1, 10, 1, &counts);
        for(auto &&x : expected) {
            expected_a.emplace_back(x);
        }
        CHECK(counts == expected_a);
    }
    SECTION("when counts are not zero") {
        std::vector<double> counts = {10, 20, 30, 40, 50};
        update_counts(1, 10, 1, &counts);
        double d = 10.0;
        for(auto &&x : expected) {
            expected_a.emplace_back(x + d);
            d += 10.0;
        }
        CHECK(counts == expected_a);
    }
    SECTION("when weight is not 1") {
        std::vector<double> counts = {0, 0, 0, 0, 0};
        update_counts(1, 10, 0.1, &counts);
        for(auto &&x : expected) {
            expected_a.emplace_back(x * 0.1);
        }
        CHECK(counts == expected_a);
    }
}
#endif

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

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("update_bins maps a/(a+b) to a bin in the range of [0,N]/N") {
    using namespace Catch::literals;

    std::vector<double> bins(4, 0.0);

    SECTION("a+b == 0.0") {
        update_bins(0, 0, 2, &bins);
        REQUIRE(bins.size() == 4);
        CHECK(bins[0] == 0.0);
        CHECK(bins[1] == 0.0);
        CHECK(bins[2] == 0.0);
        CHECK(bins[3] == 0.0);
    }
    SECTION("a == 1 && b == 0") {
        update_bins(1, 0, 2, &bins);
        REQUIRE(bins.size() == 4);
        CHECK(bins[0] == 0.0);
        CHECK(bins[1] == 0.0);
        CHECK(bins[2] == 0.0);
        CHECK(bins[3] == 2.0);
    }
    SECTION("a == 0 && b == 1") {
        update_bins(0, 1, 2, &bins);
        REQUIRE(bins.size() == 4);
        CHECK(bins[0] == 2.0);
        CHECK(bins[1] == 0.0);
        CHECK(bins[2] == 0.0);
        CHECK(bins[3] == 0.0);
    }
    SECTION("after many calls") {
        update_bins(0, 1, 2, &bins);
        update_bins(1, 0, 2, &bins);
        update_bins(1, 2, 2, &bins);
        update_bins(1.1, 2, 2, &bins);
        REQUIRE(bins.size() == 4);
        CHECK(bins[0] == 2.0);
        CHECK(bins[1] == 4.0);
        CHECK(bins[2] == 0.0);
        CHECK(bins[3] == 2.0);
    }
}
#endif

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("fold_histogram folds second half of a vector onto the first") {
    using namespace Catch::literals;
    SECTION("fold histogram with even number of elements") {
        std::vector<double> v = {1, 2, 3, 10, 20, 30};
        fold_histogram(&v);
        REQUIRE(v.size() == 3);
        CHECK(v[0] == 31.0_a);
        CHECK(v[1] == 22.0_a);
        CHECK(v[2] == 13.0_a);
    }
    SECTION("fold histogram with odd number of elements") {
        std::vector<double> v = {1, 2, 3, 10, 20};
        fold_histogram(&v);
        REQUIRE(v.size() == 3);
        CHECK(v[0] == 21.0_a);
        CHECK(v[1] == 12.0_a);
        CHECK(v[2] == 3.0_a);
    }
}
#endif

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

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("quality_to_p01 converts a quality value to a [0,1] value.") {
    using namespace Catch::literals;
    SECTION("quality_to_p01(int)") {
        CHECK(quality_to_p01(0) == 0.0_a);
        CHECK(quality_to_p01(10) == 0.9_a);
        CHECK(quality_to_p01(20) == 0.99_a);
        CHECK(quality_to_p01(30) == 0.999_a);
        CHECK(quality_to_p01(40) == 0.9999_a);
        CHECK(quality_to_p01(255) == 1.0_a);
    }
    SECTION("quality_to_p01(double)") {
        CHECK(quality_to_p01(0.0f) == 0.0_a);
        CHECK(quality_to_p01(10.0f) == 0.9_a);
        CHECK(quality_to_p01(20.0f) == 0.99_a);
        CHECK(quality_to_p01(30.0f) == 0.999_a);
        CHECK(quality_to_p01(40.0f) == 0.9999_a);
        CHECK(quality_to_p01(255.0f) == 1.0_a);
    }
}

TEST_CASE("phred_to_p01 converts a phred value to a [0,1] value.") {
    using namespace Catch::literals;
    auto zero = Approx(0.0).margin(std::numeric_limits<float>::epsilon() * 100);
    SECTION("phred_to_p01(int)") {
        CHECK(phred_to_p01(0) == 1.0_a);
        CHECK(phred_to_p01(10) == 0.1_a);
        CHECK(phred_to_p01(20) == 0.01_a);
        CHECK(phred_to_p01(30) == 0.001_a);
        CHECK(phred_to_p01(40) == 0.0001_a);
        CHECK(phred_to_p01(255) == zero);
    }
    SECTION("phred_to_p01(double)") {
        CHECK(phred_to_p01(0.0f) == 1.0_a);
        CHECK(phred_to_p01(10.0f) == 0.1_a);
        CHECK(phred_to_p01(20.0f) == 0.01_a);
        CHECK(phred_to_p01(30.0f) == 0.001_a);
        CHECK(phred_to_p01(40.0f) == 0.0001_a);
        CHECK(phred_to_p01(255.0f) == zero);
    }
}
#endif

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

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("output_header generates a header containing program information") {
    using Catch::Matchers::StartsWith;

    auto get_header_lines = [](const sofos_params_t &params,
                               const std::vector<const char *> &paths) -> std::vector<std::string> {
        Sofos sofos{params};
        std::stringstream str;
        output_header(str, sofos, paths);

        std::vector<std::string> lines;
        std::string token;
        while(std::getline(str, token)) {
            lines.push_back(token);
        }
        return lines;
    };

    SECTION("default output") {
        auto lines = get_header_lines({}, {});
        REQUIRE(lines.size() == 9);
        CHECK(lines[0] == "#SoFoS v2.0");
        for(auto &&line : lines) {
            CHECK_THAT(line, StartsWith("#"));
        }
    }
    SECTION("when there is one path") {
        auto lines = get_header_lines({}, {"-"});
        CHECK(lines[0] == "#SoFoS v2.0");
        REQUIRE(lines.size() == 10);
        for(auto &&line : lines) {
            CHECK_THAT(line, StartsWith("#"));
        }
    }
    SECTION("when everything is enabled") {
        sofos_params_t params;
        params.flag_use_gp = true;
        params.flag_phred_gp = true;
        params.error_rate = 1.0;
        params.ploidy = 2.1;
        params.zero_count = 100;

        auto lines = get_header_lines(params, {"-", "a", "b"});
        CHECK(lines[0] == "#SoFoS v2.0");
        REQUIRE(lines.size() == 16);
        for(auto &&line : lines) {
            CHECK_THAT(line, StartsWith("#"));
        }
    }
}
#endif

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

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("output_body generates a csv table containing column names and values") {
    sofos_params_t params;
    params.size = 4;
    Sofos sofos{params};
    std::stringstream str;
    output_body(str, sofos);

    std::vector<std::vector<std::string>> table;
    std::string line;
    while(std::getline(str, line, '\n')) {
        std::stringstream row;
        row.str(line);
        std::string token;
        std::vector<std::string> tokens;
        while(std::getline(row, token, ',')) {
            tokens.push_back(token);
        }
        table.push_back(tokens);
    }

    REQUIRE(table.size() == 6);
    CHECK(table[0].size() == 4);
    CHECK(table[1].size() == 4);
    CHECK(table[2].size() == 4);
    CHECK(table[3].size() == 4);
    CHECK(table[4].size() == 4);
    CHECK(table[5].size() == 4);
}
#endif

#ifdef SOFOS_UNIT_TESTS
TEST_CASE("Retrieving values from a VCF file") {
    using Catch::Matchers::Equals;
    char header_str[] =
        "##fileformat=VCFv4.2\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##INFO=<ID=TAG1,Number=.,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=TAG2,Number=.,Type=Integer,Description=\"\">\n"
        "##INFO=<ID=TAG3,Number=.,Type=String,Description=\"\">\n"
        "##contig=<ID=1,length=10000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n";

    std::unique_ptr<bcf1_t, detail::bcf_free_t> record_{bcf_init()};
    REQUIRE((record_));

    std::unique_ptr<bcf_hdr_t, detail::header_free_t> header_{bcf_hdr_init("w")};
    REQUIRE((header_));

    auto header = header_.get();
    auto record = record_.get();

    int n;
    int ret = bcf_hdr_parse(header, header_str);
    REQUIRE(ret == 0);
    kstring_obj_t kstr;

    SECTION("when reading INFO integer tag") {
        auto buffer = make_buffer<int32_t>(1);
        const char line1[] = "1\t1\t.\tA\tC\t.\t.\tTAG1=1;TAG2=1,1,1,1,1,1,1,1,1,1";
        ret = kputs(line1, &kstr);
        REQUIRE(ret > 0);
        ret = vcf_parse(&kstr, header, record);
        REQUIRE(ret == 0);

        n = get_info_int32(header, record, "TAG1", &buffer);
        REQUIRE(n == 1);
        CHECK(buffer.data[0] == 1);

        n = get_info_int32(header, record, "TAG2", &buffer);
        REQUIRE(n == 10);
        for(int i = 0; i < 10; ++i) {
            CHECK(buffer.data[i] == 1);
        }
    }
    SECTION("when reading INFO string tag") {
        auto buffer = make_buffer<char>(1);
        const char line1[] = "1\t1\t.\tA\tC\t.\t.\tTAG3=aaaaaaaaaaaaaaaaaaaaaaaa";
        ret = kputs(line1, &kstr);
        REQUIRE(ret > 0);
        ret = vcf_parse(&kstr, header, record);
        REQUIRE(ret == 0);

        n = get_info_string(header, record, "TAG3", &buffer);
        REQUIRE(n == 24);
        CHECK_THAT(buffer.data.get(), Equals("aaaaaaaaaaaaaaaaaaaaaaaa"));
    }
}
#endif
