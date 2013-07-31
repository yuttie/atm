#ifndef PURITY_MAXIMAL_SUBSTRINGS_HPP
#define PURITY_MAXIMAL_SUBSTRINGS_HPP

#include <vector>
#include <boost/iterator/iterator_facade.hpp>
#include "branching_substrings.hpp"


template <class RandomAccessRange, class Index>
struct purity_maximal_substrings : public branching_substrings<RandomAccessRange, Index> {
private:
    using base_type = branching_substrings<RandomAccessRange, Index>;

public:
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
            : parent_(0), i_(-1)
        {}

        substring_iterator(const purity_maximal_substrings* parent, const int i)
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

        const purity_maximal_substrings* parent_;
        int i_;
    };

public:
    using iterator       = substring_iterator<substr>;
    using const_iterator = substring_iterator<const substr>;

    purity_maximal_substrings(const RandomAccessRange& input, const size_t alphabet_size)
        : base_type(input, alphabet_size),
          selected_node_indices_()
    {
        // traverse nodes by suffix links, find descending node sequences and
        // make a node group for each sequence
        int num_groups = 0;
        std::vector<int> group_ids(sast_.size(), -1);
        std::vector<std::vector<int>> groups;
        for (int i = 0; i < sast_.size(); ++i) {
            const int cid = get_group_id(i, num_groups, group_ids);
            groups.resize(num_groups);
            groups[cid].push_back(i);
        }

        // find the node with the best purity for each group
        for (int i = 0; i < num_groups; ++i) {
            double max_purity = -1;
            int index_max;
            for (int j = 0; j < groups[i].size(); ++j) {
                const double purity = strict_purity(sast_.begin() + groups[i][j]);
                if (purity > max_purity) {
                    max_purity = purity;
                    index_max = j;
                }
            }
            selected_node_indices_.push_back(groups[i][index_max]);
        }
    }

    iterator begin() { return iterator(this, 0); }
    iterator end()   { return iterator(this, selected_node_indices_.size()); }
    const_iterator begin() const { return const_iterator(this, 0); }
    const_iterator end()   const { return const_iterator(this, selected_node_indices_.size()); }

protected:
    using sast_type = typename base_type::sast_type;
    using base_type::strict_purity;

    int get_group_id(const int i, int& num_groups, std::vector<int>& group_ids) const {
        if (group_ids[i] >= 0) {
            return group_ids[i];
        }
        else {
            typename sast_type::const_iterator n = sast_.begin() + i;
            const auto freq_substr = n->frequency();

            const auto m = n.suffix();
            const int j = m - sast_.begin();
            const auto freq_subsubstr = m->frequency();

            if (freq_subsubstr == freq_substr && strict_purity(m) <= strict_purity(n)) {
                group_ids[i] = get_group_id(j, num_groups, group_ids);
                return group_ids[i];
            }
            else {
                group_ids[i] = num_groups++;
                return group_ids[i];
            }
        }
    }

    using base_type::sast_;
    std::vector<Index> selected_node_indices_;
};


#endif  /* PURITY_MAXIMAL_SUBSTRINGS_HPP */
