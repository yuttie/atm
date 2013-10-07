#ifndef ATM_BRANCHING_SUBSTRINGS_HPP
#define ATM_BRANCHING_SUBSTRINGS_HPP

#include <cmath>
#include <map>
#include <vector>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/range/value_type.hpp>
#include "sast/sast.hpp"


namespace atm {

template <class RandomAccessRange, class Index>
struct branching_substrings {
    using range_type = RandomAccessRange;
    using char_type = typename boost::range_value<RandomAccessRange>::type;
    using index_type = Index;

protected:
    using sast_type = sast::sast<RandomAccessRange, index_type>;

public:
    struct substr {
        using iterator       = typename boost::range_iterator<RandomAccessRange>::type;
        using const_iterator = typename boost::range_const_iterator<RandomAccessRange>::type;

        index_type pos()           const { return i_->pos(); }
        boost::sub_range<const std::vector<index_type>> allpos() const { return i_->allpos(); }
        index_type length()        const { return i_->length(); }
        index_type frequency()     const { return i_->frequency(); }
        double     spurity()       const { return parent_->strict_purity(i_); }
        double     lpurity()       const { return parent_->loose_purity(i_); }
        double     luniversality() const { return parent_->left_universality(i_); }
        double     runiversality() const { return parent_->right_universality(i_); }

        iterator begin() { return boost::begin(parent_->input_) + pos(); }
        iterator end()   { return boost::begin(parent_->input_) + pos() + length(); }
        const_iterator begin() const { return boost::const_begin(parent_->input_) + pos(); }
        const_iterator end()   const { return boost::const_begin(parent_->input_) + pos() + length(); }

        substr(const branching_substrings* parent, typename sast_type::const_iterator i)
            : parent_(parent), i_(i)
        {}

    private:
        const branching_substrings* parent_;
        typename sast_type::const_iterator i_;
    };

private:
    template <class> struct substring_iterator;

public:
    using iterator       = substring_iterator<substr>;
    using const_iterator = substring_iterator<const substr>;

    branching_substrings(const sast_type& sast)
        : sast_(sast),
          input_(sast_.input()),
          count_(sast_.size(), 0), // initialize the count table with an "undefined" value.
          recip_(sast_.size(), 0)  // 正数は計算結果、それ以外は未計算を表わす。
    {}

    iterator begin() { return iterator(this, sast_.begin()); }
    iterator end()   { return iterator(this, sast_.end()); }
    const_iterator begin() const { return const_iterator(this, sast_.begin()); }
    const_iterator end()   const { return const_iterator(this, sast_.end()); }

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
            : parent_(0), i_()
        {}

        substring_iterator(const branching_substrings* parent, typename sast_type::const_iterator i)
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
            return substr(parent_, i_);
        }

        const branching_substrings* parent_;
        typename sast_type::const_iterator i_;
    };

protected:
    uint64_t get_count(typename sast_type::const_iterator n) const {
        const int i = n - sast_.begin();

        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        // substrと同じ出現回数のsub-substrの数count[i]を数える。
        if (count_[i] > 0) {
            return count_[i];
        }
        else {
            const auto freq_substr = n->frequency();

            // substrと同じ出現回数のsub-substrを数える。
            uint64_t count = 0;
            {
                // substrの末尾を0文字以上削って得られるsub-substrについて考える。
                // ノードiに対応する部分文字列をsubstr[i]とすると、substr[i]の末尾
                // を削って得られる部分文字列の内で、substr[i]と同じ頻度をもつもの
                // はsuffix tree上ではノードiにまとめられている
                // （分岐が無い <=> 頻度が同じ）。

                // ノードiの親ノードjを見つける。
                const auto p = n.parent();

                // substrの末尾を0文字以上削って得られるsub-substrの内で、出現
                // 回数がsubstrと同じものの数はd[i] - d[j]である。
                count += n->length() - p->length();;
            }
            {
                // substrの先頭を1文字以上削ったsub-substrを考える。
                const auto m = n.suffix();
                const auto freq_subsubstr = m->frequency();
                if (freq_subsubstr == freq_substr) {
                    count += get_count(m);
                }
            }

            // memoize
            count_[i] = count;
            return count_[i];
        }
    }

    double strict_purity(typename sast_type::const_iterator n) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto len_substr  = n->length();

        // substrと同じ出現回数のsub-substrを数える。
        const uint64_t count = get_count(n);

        // strict purity of substr
        const uint64_t num_subsubstrs = static_cast<uint64_t>(1 + len_substr) * len_substr / 2;
        const double spurity = static_cast<double>(count) / num_subsubstrs;

        return spurity;
    }

    double get_reciprocal(typename sast_type::const_iterator n) const {
        const int i = n - sast_.begin();

        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        // substrの全部分文字列の頻度の逆数の総和recip[i]を求める。
        if (recip_[i] > 0) {
            return recip_[i];
        }
        else {
            double recip = 0;
            {
                // substrの末尾を0文字以上削って得られるsub-substrについて考える。
                for (auto m = n, p = n.parent(); m != sast_.end(); m = p, p = p.parent()) {
                    const auto num_subsubstrs_of_same_frequency = m->length() - p->length();
                    const auto freq_subsubstr = m->frequency();
                    const double r = 1.0 / freq_subsubstr;
                    recip += num_subsubstrs_of_same_frequency * r;
                }
            }
            {
                // substrの先頭を1文字以上削ったsub-substrを考える。
                const auto m = n.suffix();
                if (m != sast_.end()) {
                    recip += get_reciprocal(m);
                }
            }

            // memoize
            recip_[i] = recip;
            return recip_[i];
        }
    }

    double loose_purity(typename sast_type::const_iterator n) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto freq_substr = n->frequency();
        const auto len_substr  = n->length();

        const double rel = freq_substr * get_reciprocal(n);

        // loose purity of substr
        const uint64_t num_subsubstrs = static_cast<uint64_t>(1 + len_substr) * len_substr / 2;
        const double lpurity = rel / num_subsubstrs;

        return lpurity;
    }

    std::map<char_type, int> left_extensions(typename sast_type::const_iterator n) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。

        std::map<char_type, int> char_dist;
        for (auto pos : n->allpos()) {
            const auto& c = boost::const_begin(input_)[pos - 1];
            char_dist[c] += 1;
        }

        return char_dist;
    }

    std::map<char_type, int> right_extensions(typename sast_type::const_iterator n) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto len_substr = n->length();

        std::map<char_type, int> char_dist;
        for (auto pos : n->allpos()) {
            const auto& c = boost::const_begin(input_)[pos + len_substr];
            char_dist[c] += 1;
        }

        return char_dist;
    }

    double left_universality(typename sast_type::const_iterator n) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto freq_substr = n->frequency();

        std::map<char_type, int> char_dist = left_extensions(n);

        double e = 0;
        for (const auto& kv : char_dist) {
            const double p = static_cast<double>(kv.second) / freq_substr;
            e += -p * std::log(p);
        }
        const double u = 1 - std::exp(-e);

        return u;
    }

    double right_universality(typename sast_type::const_iterator n) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto freq_substr = n->frequency();

        std::map<char_type, int> char_dist = right_extensions(n);

        double e = 0;
        for (const auto& kv : char_dist) {
            const double p = static_cast<double>(kv.second) / freq_substr;
            e += -p * std::log(p);
        }
        const double u = 1 - std::exp(-e);

        return u;
    }

    const sast_type& sast_;
    const RandomAccessRange& input_;
    mutable std::vector<uint64_t> count_;
    mutable std::vector<double>   recip_;
};

}  // namespace atm


#endif  /* ATM_BRANCHING_SUBSTRINGS_HPP */
