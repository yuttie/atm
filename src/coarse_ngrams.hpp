#ifndef COARSE_NGRAMS_H
#define COARSE_NGRAMS_H

#include "coarse_substrings.hpp"


template <class RandomAccessRange, class Index>
struct coarse_ngrams : public coarse_substrings<RandomAccessRange, Index> {
    typedef coarse_substrings<RandomAccessRange, Index> base_type;
    using typename base_type::index_type;
    using typename base_type::substr;
    using typename base_type::iterator;
    using typename base_type::const_iterator;

    coarse_ngrams(const RandomAccessRange& input, const size_t alphabet_size, const int r, const int n)
        : base_type(input, alphabet_size, r), r_(r), n_(n)
    {}

    iterator begin() { return iterator(this, r_, 0, n_); }
    iterator end()   { return iterator(this, r_, 0, n_ - 1); }
    const_iterator begin() const { return const_iterator(this, r_, 0, n_); }
    const_iterator end()   const { return const_iterator(this, r_, 0, n_ - 1); }

private:
    const int r_;
    const int n_;
};


#endif  /* COARSE_NGRAMS_H */
