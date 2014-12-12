#ifndef ATM_CHUNKED_SLIDING_WINDOW_HPP
#define ATM_CHUNKED_SLIDING_WINDOW_HPP

#include <algorithm>
#include <boost/range/begin.hpp>
#include <boost/range/size.hpp>
#include "atm/substrings_from_longest.hpp"


namespace atm {

template <class RandomAccessRange, class Index>
struct chunked_sliding_window : public substrings_from_longest<RandomAccessRange, Index> {
private:
    using base_type = substrings_from_longest<RandomAccessRange, Index>;

public:
    using typename base_type::range_type;
    using typename base_type::char_type;
    using typename base_type::index_type;
    using typename base_type::substr;

protected:
    using typename base_type::sast_type;

private:
    template <class> struct substring_iterator;

public:
    using iterator       = substring_iterator<typename base_type::substr>;
    using const_iterator = substring_iterator<const typename base_type::substr>;

    chunked_sliding_window(const sast_type& sast, const int l, const int n)
        : base_type(sast), poslens_()
    {
        const auto& input = sast.input();

        for (auto i = 0; i < input.size(); i += l) {
            poslens_.emplace_back(i, std::min(static_cast<std::size_t>(i + n * l), input.size()) - i);
        }
    }

    iterator begin() { return iterator(this, poslens_.begin()); }
    iterator end()   { return iterator(this, poslens_.end()); }
    const_iterator begin() const { return const_iterator(this, poslens_.begin()); }
    const_iterator end()   const { return const_iterator(this, poslens_.end()); }

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
            : parent_(0), it_()
        {}

        substring_iterator(const chunked_sliding_window* parent, const std::vector<std::pair<int, int>>::iterator& it)
            : parent_(parent), it_(it)
        {}

        template <class OtherValue>
        substring_iterator(const substring_iterator<OtherValue>& other)
            : parent_(other.parent_), it_(other.it_)
        {}

    private:
        friend class boost::iterator_core_access;
        template <class> friend struct substring_iterator;

        void increment() { ++it_; }
        void decrement() { --it_; }
        void advance(int n) { it_ += n; }

        int distance_to(const substring_iterator<Value>& other) const {
            return other.it_ - this->it_;
        }

        template <class OtherValue>
        bool equal(const substring_iterator<OtherValue>& other) const {
            return this->parent_ == other.parent_
                && this->it_ == other.it_;
        }

        Value dereference() const {
            const int i = it_->first;
            const int j = i + it_->second;
            return typename base_type::substr(parent_, i, j);
        }

        const chunked_sliding_window* parent_;
        std::vector<std::pair<int, int>>::iterator it_;
    };

private:
    std::vector<std::pair<int, int>> poslens_;
};

}  // namespace atm


#endif  /* ATM_CHUNKED_SLIDING_WINDOW_HPP */
