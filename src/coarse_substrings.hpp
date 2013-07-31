#ifndef COARSE_SUBSTRINGS_H
#define COARSE_SUBSTRINGS_H

#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/size.hpp>

#include "substrings_from_longest.hpp"


template <class RandomAccessRange, class Index>
struct coarse_substrings : public substrings_from_longest<RandomAccessRange, Index> {
    typedef substrings_from_longest<RandomAccessRange, Index> base_type;
    using typename base_type::index_type;
    using typename base_type::substr;

private:
    template <class Value>
    struct substring_iterator
        : public boost::iterator_facade<
            substring_iterator<Value>,
            Value,
            boost::random_access_traversal_tag,
            Value,
            int>
    {
        substring_iterator()
            : parent_(0), n_(0), i_(-1), j_(-1)
        {}

        substring_iterator(const coarse_substrings* parent, int n, int i, int j)
            : parent_(parent), n_(n), i_(i), j_(j)
        {}

    private:
        friend class boost::iterator_core_access;

        void increment() {
            if (j_ == n_) {
                const auto width = (j_ - i_) - 1;
                i_ = 0;
                j_ = 0 + width;
            }
            else {
                ++i_;
                ++j_;
            }
        }

        void decrement() {
            if (i_ == 0) {
                const auto width = (j_ - i_) + 1;
                i_ = n_ - width;
                j_ = n_;
            }
            else {
                --i_;
                --j_;
            }
        }

        void advance(int n) {
            if (n >= 0) {
                auto width = j_ - i_;
                while (n > n_ - j_) {
                    n -= (n_ - j_) + 1;
                    --width;
                    i_ = 0;
                    j_ = 0 + width;
                }
                i_ += n;
                j_ += n;
            }
            else {
                n = -n;
                auto width = j_ - i_;
                while (n > i_ - 0) {
                    n -= (i_ - 0) + 1;
                    ++width;
                    i_ = n_ - width;
                    j_ = n_;
                }
                i_ -= n;
                j_ -= n;
            }
        }

        int distance_to(const substring_iterator<Value>& other) const {
            int d = 0;
            auto width = j_ - i_;
            auto owidth = other.j_ - other.i_;
            if (owidth < width || (owidth == width && other.i_ > this->i_)) {
                int j = this->j_;
                while (owidth < width) {
                    d += (n_ - j) + 1;
                    --width;
                    j = 0 + width;
                }
                d += other.j_ - j;
            }
            else {
                int j = other.j_;
                while (owidth > width) {
                    d += (n_ - j) + 1;
                    --owidth;
                    j = 0 + owidth;
                }
                d += this->j_ - j;
                d = -d;
            }
            return d;
        }

        bool equal(const substring_iterator<Value>& other) const {
            return this->parent_ == other.parent_
                && this->i_ == other.i_
                && this->j_ == other.j_;
        }

        Value dereference() const {
            const int i = boost::size(parent_->input_) * i_ / n_;
            const int j = boost::size(parent_->input_) * j_ / n_;
            return typename base_type::substr(parent_, i, j);
        }

        const coarse_substrings* parent_;
        int n_;
        int i_;
        int j_;
    };

public:
    typedef substring_iterator<typename base_type::substr> iterator;
    typedef substring_iterator<const typename base_type::substr> const_iterator;

    coarse_substrings(const RandomAccessRange& input, const size_t alphabet_size, const int n)
        : base_type(input, alphabet_size), n_(n)
    {}

    iterator begin() { return iterator(this, n_, 0, n_); }
    iterator end()   { return iterator(this, n_, 0, 0); }
    const_iterator begin() const { return const_iterator(this, n_, 0, n_); }
    const_iterator end()   const { return const_iterator(this, n_, 0, 0); }

protected:
    using base_type::input_;
    const int n_;
};


#endif  /* COARSE_SUBSTRINGS_H */
