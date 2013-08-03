#ifndef ATM_COARSE_NGRAMS_HPP
#define ATM_COARSE_NGRAMS_HPP

#include "atm/coarse_substrings.hpp"


namespace atm {

template <class RandomAccessRange, class Index>
struct coarse_ngrams : public coarse_substrings<RandomAccessRange, Index> {
private:
    using base_type = coarse_substrings<RandomAccessRange, Index>;

public:
    using typename base_type::range_type;
    using typename base_type::char_type;
    using typename base_type::index_type;
    using typename base_type::substr;
    using typename base_type::iterator;
    using typename base_type::const_iterator;

protected:
    using typename base_type::sast_type;

public:
    coarse_ngrams(const sast_type& sast, const int r, const int n)
        : base_type(sast, r), r_(r), n_(n)
    {}

    iterator begin() { return iterator(this, r_, 0, n_); }
    iterator end()   { return iterator(this, r_, 0, n_ - 1); }
    const_iterator begin() const { return const_iterator(this, r_, 0, n_); }
    const_iterator end()   const { return const_iterator(this, r_, 0, n_ - 1); }

private:
    const int r_;
    const int n_;
};

}  // namespace atm


#endif  /* ATM_COARSE_NGRAMS_HPP */
