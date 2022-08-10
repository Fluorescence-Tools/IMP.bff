#ifndef IMPBFF_CNPORTCACHE_H
#define IMPBFF_CNPORTCACHE_H

#include <IMP/bff/bff_config.h>

IMPBFF_BEGIN_NAMESPACE


class CnPortCache{

private:

    /// Content of cache
    std::vector<unsigned char> cache_;

    /// Number of elements
    size_t n_elements_ = 0;

    /// size (in bytes) of elements
    size_t element_size = 1;

public:

    size_t size(){
        return cache_.size();
    }

    bool empty(){
        return cache_.empty();
    }

    std::vector<unsigned char>& get_bytes(){
        return cache_;
    }

    void set_bytes(std::vector<unsigned char>& v){
        cache_ = v;
    }

    void set_bytes(unsigned char *input, int n_input){
        cache_.resize(n_input);
        cache_.assign(input, input + n_input);
    }

    template<typename T>
    void get(T **output, int *n_output){
        T* o = reinterpret_cast<T *>(cache_.data());
        *output = o;
        *n_output = size() / sizeof(T);
    }

    template<typename T>
    void set(T *input, int n_input){
        size_t n_bytes = n_input * sizeof(T);
        cache_.resize(n_bytes, 0);
        memcpy(cache_.data(), input, n_bytes);
    }

};

IMPBFF_END_NAMESPACE


#endif //IMPBFF_CNPORTCACHE_H
