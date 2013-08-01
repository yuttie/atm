#ifndef SINGLE_RANGE_H
#define SINGLE_RANGE_H

#include "atm/substrings_from_longest.hpp"


namespace atm {

template <class RandomAccessRange, class Index>
struct single_range : public substrings_from_longest<RandomAccessRange, Index> {
private:
    using base_type = substrings_from_longest<RandomAccessRange, Index>;

public:
    using typename base_type::range_type;
    using typename base_type::char_type;
    using typename base_type::index_type;
    using typename base_type::substr;
    using typename base_type::iterator;
    using typename base_type::const_iterator;

    single_range(const RandomAccessRange& input, const size_t alphabet_size, const int i, const int j)
        : base_type(input, alphabet_size), i_(i), j_(j)
    {}

    iterator begin() { return iterator(this, i_, j_); }
    iterator end()   { return begin() + 1; }
    const_iterator begin() const { return const_iterator(this, i_, j_); }
    const_iterator end()   const { return begin() + 1; }

private:
    const int i_;
    const int j_;
};

}  // namespace atm


#endif  /* SINGLE_RANGE_H */
