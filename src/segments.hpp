#ifndef SEGMENTS_H
#define SEGMENTS_H

#include "coarse_substrings.hpp"


template <class RandomAccessRange, class Index>
struct Segments : public CoarseSubstrings<RandomAccessRange, Index> {
    typedef CoarseSubstrings<RandomAccessRange, Index> base_type;
    using typename base_type::index_type;
    using typename base_type::substr;
    using typename base_type::iterator;
    using typename base_type::const_iterator;

    Segments(const RandomAccessRange& input, const size_t alphabet_size, const int n)
        : base_type(input, alphabet_size, n)
    {}

    iterator begin() {
        return iterator(this, base_type::n_, 0, 1);
    }

    const_iterator begin() const {
        return const_iterator(this, base_type::n_, 0, 1);
    }
};


#endif  /* SEGMENTS_H */
