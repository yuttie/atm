#ifndef WORDS_H
#define WORDS_H

#include <boost/range/begin.hpp>
#include <boost/range/size.hpp>
#include "atm/substrings_from_longest.hpp"


namespace atm {

template <class RandomAccessRange, class Index>
struct words : public substrings_from_longest<RandomAccessRange, Index> {
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

    template <class Pred>
    words(const sast_type& sast, Pred is_word_char)
        : base_type(sast), poslens_()
    {
        const auto& input = sast.input();

        bool within_word = false;
        for (auto it = boost::const_begin(input);
             it != boost::const_end(input);
             ++it)
        {
            if (!within_word) {
                if (is_word_char(*it)) {
                    within_word = true;
                    const int pos = it - boost::begin(input);
                    poslens_.emplace_back(pos, 1);
                }
                else {
                    // do nothing
                }
            }
            else {
                if (is_word_char(*it)) {
                    ++poslens_.back().second;
                }
                else {
                    within_word = false;
                }
            }
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

        substring_iterator(const words* parent, const std::vector<std::pair<int, int>>::iterator& it)
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

        const words* parent_;
        std::vector<std::pair<int, int>>::iterator it_;
    };

private:
    std::vector<std::pair<int, int>> poslens_;
};

}  // namespace atm


#endif  /* WORDS_H */
