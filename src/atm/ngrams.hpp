#ifndef NGRAMS_H
#define NGRAMS_H

#include "atm/substrings_from_longest.hpp"


namespace atm {

template <class RandomAccessRange, class Index>
struct ngrams : public substrings_from_longest<RandomAccessRange, Index> {
private:
    using base_type = substrings_from_longest<RandomAccessRange, Index>;

public:
    using typename base_type::range_type;
    using typename base_type::char_type;
    using typename base_type::index_type;
    using typename base_type::substr;
    using typename base_type::iterator;
    using typename base_type::const_iterator;

    ngrams(const RandomAccessRange& input, const size_t alphabet_size, const int n)
        : base_type(input, alphabet_size), n_(n)
    {}

    iterator begin() { return iterator(this, 0, n_); }
    iterator end()   { return iterator(this, 0, n_ - 1); }
    const_iterator begin() const { return const_iterator(this, 0, n_); }
    const_iterator end()   const { return const_iterator(this, 0, n_ - 1); }

private:
    const int n_;
};

}  // namespace atm


#endif  /* NGRAMS_H */
