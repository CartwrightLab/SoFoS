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
#ifndef SOFOS_VCF_HPP
#define SOFOS_VCF_HPP

#include <chrono>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <memory>

// The *_free_t classes are used enable RAII on pointers created by htslib.
namespace detail {
struct buffer_free_t {
    void operator()(void* ptr) const {
        free(ptr); // NOLINT
    }
};
struct file_free_t {
    void operator()(void* ptr) const {
        hts_close(static_cast<htsFile*>(ptr));
    }
};
struct header_free_t {
    void operator()(void* ptr) const {
        bcf_hdr_destroy(static_cast<bcf_hdr_t*>(ptr));
    }
};
struct bcf_free_t {
    void operator()(void* ptr) const {
        bcf_destroy(static_cast<bcf1_t*>(ptr));
    }
};
} // namespace detail

class BcfReader {
public:
    explicit BcfReader(const char *path) {
        input_.reset(hts_open(path,"r"));
        if(!input_) {
            throw std::runtime_error(std::string{"unable to open input file: '"} + path + "'.");
        }
        header_.reset(bcf_hdr_read(input_.get()));
        if(!header_) {
            throw std::invalid_argument("unable to read header from input.");
        }
    }

    const bcf_hdr_t* header() const { return header_.get(); }

    template<typename callback_t>
    void operator()(callback_t callback);

protected:
    std::unique_ptr<htsFile,detail::file_free_t> input_;
    std::unique_ptr<bcf_hdr_t,detail::header_free_t> header_;
};


template<typename callback_t>
void BcfReader::operator()(callback_t callback) {
    std::unique_ptr<bcf1_t,detail::bcf_free_t> record{bcf_init()};
    if(!record) {
        throw std::invalid_argument("unable to allocate vcf record.");
    }
    // process all sites
    while(bcf_read(input_.get(), header_.get(), record.get()) == 0) {
        callback(record.get(), header());
    }
}

// Templates and functions for handling buffers used by htslib
template<typename T>
struct buffer_t { // NOLINT(cppcoreguidelines-pro-type-member-init)
    std::unique_ptr<T[],detail::buffer_free_t> data;
    int capacity; 
};

template<typename T>
inline
buffer_t<T> make_buffer(int sz) {
    void *p = std::malloc(sizeof(T)*sz); // NOLINT(cppcoreguidelines-owning-memory, cppcoreguidelines-no-malloc)
    if(p == nullptr) {
        throw std::bad_alloc{};
    }
    return {std::unique_ptr<T[], detail::buffer_free_t>{static_cast<T*>(p)}, sz};
}

// htslib may call realloc on our pointer. When using a managed buffer,
// we need to check to see if it needs to be updated.
inline
int get_info_string(const bcf_hdr_t *header, bcf1_t *record,
    const char *tag, buffer_t<char>* buffer)
{
    char *p = buffer->data.get();
    int n = bcf_get_info_string(header, record, tag, &p, &buffer->capacity); // NOLINT
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
    int n = bcf_get_info_int32(header, record, tag, &p, &buffer->capacity); // NOLINT
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
    int n = bcf_get_format_float(header, record, tag, &p, &buffer->capacity); // NOLINT
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
    int n = bcf_get_genotypes(header, record, &p, &buffer->capacity); // NOLINT
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

#endif // SOFOS_VCF_HPP
