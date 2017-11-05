#ifndef ATM_BLUMER_SUBSTRINGS_HPP
#define ATM_BLUMER_SUBSTRINGS_HPP

#include <vector>
#include <boost/iterator/iterator_facade.hpp>
#include "atm/branching_substrings.hpp"


namespace atm {

template <class RandomAccessRange, class Index>
struct blumer_substrings : public branching_substrings<RandomAccessRange, Index> {
private:
    using base_type = branching_substrings<RandomAccessRange, Index>;

public:
    using typename base_type::range_type;
    using typename base_type::char_type;
    using typename base_type::index_type;
    using typename base_type::substr;

protected:
    using sast_type = typename base_type::sast_type;

private:
    template <class> struct substring_iterator;

public:
    using iterator       = substring_iterator<substr>;
    using const_iterator = substring_iterator<const substr>;

    blumer_substrings(const sast_type& sast)
        : base_type(sast),
          selected_node_indices_()
    {
        // split the nodes into Blumer's equivalence classes
        int num_classes = 0;
        std::vector<int> class_ids(sast_.size(), -1);
        std::vector<std::vector<int>> classes;
        for (int i = 0; i < sast_.size(); ++i) {
            const int cid = get_class_id(i, num_classes, class_ids);
            classes.resize(num_classes);
            classes[cid].push_back(i);
        }

        // find the node corresponding to the longest substring in each class
        for (int i = 0; i < num_classes; ++i) {
            auto j = std::max_element(classes[i].begin(), classes[i].end(), [&](int a, int b) {
                    return (sast_.begin() + a)->length() < (sast_.begin() + b)->length();
                });
            selected_node_indices_.push_back(*j);
        }
    }

    iterator begin() { return iterator(this, 0); }
    iterator end()   { return iterator(this, selected_node_indices_.size()); }
    const_iterator begin() const { return const_iterator(this, 0); }
    const_iterator end()   const { return const_iterator(this, selected_node_indices_.size()); }

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
            : parent_(0), i_(-1)
        {}

        substring_iterator(const blumer_substrings* parent, const int i)
            : parent_(parent), i_(i)
        {}

        template <class OtherValue>
        substring_iterator(substring_iterator<OtherValue> const& other)
            : parent_(other.parent_), i_(other.i_)
        {}

    private:
        friend class boost::iterator_core_access;
        template <class> friend struct substring_iterator;

        void increment() { ++i_; }

        void decrement() { --i_; }

        void advance(int n) { i_ += n; }

        int distance_to(const substring_iterator<Value>& other) const { return other.i_ - this->i_; }

        template <class OtherValue>
        bool equal(const substring_iterator<OtherValue>& other) const {
            return this->parent_ == other.parent_ && this->i_ == other.i_;
        }

        Value dereference() const {
            return substr(parent_, parent_->sast_.begin() + parent_->selected_node_indices_[i_]);
        }

        const blumer_substrings* parent_;
        int i_;
    };

protected:
    int get_class_id(const int i, int& num_classes, std::vector<int>& class_ids) const {
        if (class_ids[i] >= 0) {
            return class_ids[i];
        }
        else {
            typename sast_type::const_iterator n = sast_.begin() + i;
            const auto freq_substr = n->frequency();

            const auto m = n.suffix();
            const int j = m - sast_.begin();
            const auto freq_subsubstr = m->frequency();

            if (freq_subsubstr == freq_substr) {
                class_ids[i] = get_class_id(j, num_classes, class_ids);
                return class_ids[i];
            }
            else {
                class_ids[i] = num_classes++;
                return class_ids[i];
            }
        }
    }

    using base_type::sast_;
    std::vector<Index> selected_node_indices_;
};

}  // namespace atm


#endif  /* ATM_BLUMER_SUBSTRINGS_HPP */
