#ifndef NGRAMS_H
#define NGRAMS_H

#include "substrings_from_longest.hpp"


template <class Char, class Index>
struct NGrams : public SubstringsFromLongest<Char, Index> {
    typedef SubstringsFromLongest<Char, Index> base_type;
    using typename base_type::index_type;
    using typename base_type::substr;
    using typename base_type::iterator;
    using typename base_type::const_iterator;

    NGrams(const std::vector<Char>& input, const size_t alphabet_size, const int n)
        : base_type(input, alphabet_size), n_(n)
    {}

    iterator begin() {
        return iterator(this, 0, n_);
    }

    iterator end() {
        return iterator(this, 0, n_ - 1);
    }

    const_iterator begin() const {
        return const_iterator(this, 0, n_);
    }

    const_iterator end() const {
        return const_iterator(this, 0, n_ - 1);
    }

private:
    const int n_;
};


#endif  /* NGRAMS_H */
