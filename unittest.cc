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

namespace Catch {
namespace Detail {
class Approx;
}  // namespace Detail
}  // namespace Catch

template <typename T>
bool operator==(const std::vector<T> &a, const std::vector<Catch::Detail::Approx> &b) {
    if(a.size() != b.size()) {
        return false;
    }
    for(size_t i = 0; i < a.size(); ++i) {
        if(a[i] != b[i]) {
            return false;
        }
    }
    return true;
}

#define CATCH_CONFIG_MAIN
#include "contrib/catch.hpp"

#include <sstream>

struct kstring_obj_t : public kstring_t {  // NOLINT
    kstring_obj_t() : kstring_t{} {
        l = 0;
        m = 0;
        s = nullptr;
    }

    ~kstring_obj_t() {
        if(s != nullptr) {
            free(s);  // NOLINT
        }
    }
};

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

TEST_CASE("Sofos::RescaleBcf rescales SFS to a new sample size") {
    using Catch::literals::operator""_a;

    g_sofos_quiet = true;

    const char header_str[] =
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
        data += &header_str[0];
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
        const char *str = "1\t1\t.\tA\tC\t.\t.\tAA=C;AN=10;AC=9\n";
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
        const char *str = "1\t1\t.\tA\tC\t.\t.\tAN=10;AC=9\n";
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
        const char *str = "1\t1\t.\tA\tC\t.\t.\tAA=G;AN=10;AC=1\n";
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
        const char *str = "1\t1\t.\tA\tC\t.\t.\t.\n";
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
        const char *str =
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
        const char *str = "1\t1\t.\tA\tC\t.\t.\tAA=C;AN=10;AC=0\tGT:GP\t0/0:1,0,0\t0/1:0,1,0\t1/1:0,0,1\n";
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
        const char *str = "1\t1\t.\tA\tC\t.\t.\tAA=C;AN=10;AC=0\tGT\t0/0\t0/1\t1/1\n";
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
        const char *str = "1\t1\t.\tA\tC\t.\t.\tAA=C;AN=10;AC=0\tGP\t1,0,0\t0,1,0\t0,0,1\n";
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
        const char *str = "1\t1\t.\tA\tC\t.\t.\tAA=C;AN=0;AC=0\n";
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
        const char *str = "1\t1\t.\tA\tC\t.\t.\tAA=C;AN=10;AC=0\tGT:GP\t0/0:0,100,100\t0/1:100,0,100\t1/1:100,100,0\n";
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
        const char *str = "1\t1\t.\tA\tC\t.\t.\tAA=G;AN=10;AC=1\n";
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
        const char *str = "1\t1\t.\tA\tC\t.\t.\tAA=A;AN=10;AC=1\n";
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
        const char *str = "1\t1\t.\tA\tC\t.\t.\tAA=A;AN=10;AC=1;AAQ=3\n";
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
        const char *str = "1\t1\t.\tA\tC\t.\t.\tAA=A;AN=10;AC=1;AAQ=10\n";
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
        const char *str = "1\t1\t.\tA\tC\t.\t.\tAA=G;AN=10;AC=1\n";
        auto hist = rescale(params, str);
        std::vector<Approx> expected_prior = {0.6666667_a, 0.3333333_a};
        std::vector<Approx> expected_observed = {1.0_a, 0.0_a};
        std::vector<Approx> expected_posterior = {0.7435897_a, 0.25641026_a};
        CHECK(hist.col(0) == expected_prior);
        CHECK(hist.col(1) == expected_observed);
        CHECK(hist.col(2) == expected_posterior);
    }
}

TEST_CASE("calculate_aa_t calculates the ancestral allele id") {
    using Catch::literals::operator""_a;

    char header_str[] =
        "##fileformat=VCFv4.2\n"
        "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
        "##INFO=<ID=AA,Number=1,Type=String,Description=\"Ancestral Allele\">\n"
        "##INFO=<ID=AAQ,Number=1,Type=Integer,Description=\"AA Quality\">\n"
        "##contig=<ID=1,length=10000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tA\tB\tC\n";

    std::unique_ptr<bcf1_t, detail::bcf_free_t> record_{bcf_init()};
    REQUIRE((record_));  // NOLINT

    std::unique_ptr<bcf_hdr_t, detail::header_free_t> header_{bcf_hdr_init("w")};
    REQUIRE((header_));  // NOLINT

    auto header = header_.get();
    auto record = record_.get();

    int hret = bcf_hdr_parse(header, &header_str[0]);
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
    int anc = 0;
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

TEST_CASE("calculate_aaq_t calculates the ancestral allele error_rate") {
    using Catch::literals::operator""_a;

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

    int hret = bcf_hdr_parse(header, &header_str[0]);
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
    double error_rate = 0.0;
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

    int ret = bcf_hdr_parse(header, &header_str[0]);
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

TEST_CASE("calculate_af() calculates allele counts from GP values") {
    using Catch::literals::operator""_a;

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

    int ret = bcf_hdr_parse(header, &header_str[0]);
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

/*
#Rcode for generating expectations
library(rmutil)
dbetabinom(0:n,n,a/(a+b),a+b)
*/
TEST_CASE("update_counts adds resampled data to a vector") {
    using Catch::literals::operator""_a;

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

TEST_CASE("update_bins maps a/(a+b) to a bin in the range of [0,N]/N") {
    using Catch::literals::operator""_a;

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

TEST_CASE("fold_histogram folds second half of a vector onto the first") {
    using Catch::literals::operator""_a;

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

TEST_CASE("quality_to_p01 converts a quality value to a [0,1] value.") {
    using Catch::literals::operator""_a;

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
    using Catch::literals::operator""_a;

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
    REQUIRE((record_));  // NOLINT

    std::unique_ptr<bcf_hdr_t, detail::header_free_t> header_{bcf_hdr_init("w")};
    REQUIRE((header_));  // NOLINT

    auto header = header_.get();
    auto record = record_.get();

    int n;
    int ret = bcf_hdr_parse(header, &header_str[0]);
    REQUIRE(ret == 0);
    kstring_obj_t kstr;

    SECTION("when reading INFO integer tag") {
        auto buffer = make_buffer<int32_t>(1);
        const char *line1 = "1\t1\t.\tA\tC\t.\t.\tTAG1=1;TAG2=1,1,1,1,1,1,1,1,1,1";
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
        const char *line1 = "1\t1\t.\tA\tC\t.\t.\tTAG3=aaaaaaaaaaaaaaaaaaaaaaaa";
        ret = kputs(line1, &kstr);
        REQUIRE(ret > 0);
        ret = vcf_parse(&kstr, header, record);
        REQUIRE(ret == 0);

        n = get_info_string(header, record, "TAG3", &buffer);
        REQUIRE(n == 24);
        CHECK_THAT(buffer.data.get(), Equals("aaaaaaaaaaaaaaaaaaaaaaaa"));
    }
}
