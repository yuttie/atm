#ifndef SEGMENTS_H
#define SEGMENTS_H

#include "atm/coarse_substrings.hpp"


namespace atm {

template <class RandomAccessRange, class Index>
struct segments : public coarse_substrings<RandomAccessRange, Index> {
private:
    using base_type = coarse_substrings<RandomAccessRange, Index>;

public:
    using typename base_type::range_type;
    using typename base_type::char_type;
    using typename base_type::index_type;
    using typename base_type::substr;
    using typename base_type::iterator;
    using typename base_type::const_iterator;

    segments(const RandomAccessRange& input, const size_t alphabet_size, const int n)
        : base_type(input, alphabet_size, n)
    {}

    iterator begin() { return iterator(this, base_type::n_, 0, 1); }
    const_iterator begin() const { return const_iterator(this, base_type::n_, 0, 1); }
};

}  // namespace atm


#endif  /* SEGMENTS_H */