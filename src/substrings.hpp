#ifndef SUBSTRINGS_HPP
#define SUBSTRINGS_HPP

#include <cmath>
#include <map>
#include <stack>
#include <stdexcept>
#include <vector>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/size.hpp>
#include <boost/range/value_type.hpp>
#include "sast.hpp"


template <class RandomAccessRange, class Index>
struct substrings {
protected:
    using char_type = typename boost::range_value<RandomAccessRange>::type;
    using sast_type = sast<RandomAccessRange, Index>;

public:
    typedef Index index_type;
    struct substr {
        using iterator       = typename boost::range_iterator<RandomAccessRange>::type;
        using const_iterator = typename boost::range_const_iterator<RandomAccessRange>::type;

        index_type              pos()           const { return i_->pos(); }
        std::vector<index_type> allpos()        const { return i_->allpos(); }
        index_type              length()        const { return i_->length() - ii_; }
        index_type              frequency()     const { return i_->frequency(); }
        double                  spurity()       const { return parent_->strict_purity(i_, ii_); }
        double                  lpurity()       const { return parent_->loose_purity(i_, ii_); }
        double                  luniversality() const { return parent_->left_universality(i_, ii_); }
        double                  runiversality() const { return parent_->right_universality(i_, ii_); }

        iterator begin() { return boost::begin(parent_->input_) + pos(); }
        iterator end()   { return boost::begin(parent_->input_) + pos() + length(); }
        const_iterator begin() const { return boost::begin(parent_->input_) + pos(); }
        const_iterator end()   const { return boost::begin(parent_->input_) + pos() + length(); }

        substr(const substrings* parent, typename sast_type::const_iterator i, int ii)
            : parent_(parent), i_(i), ii_(ii)
        {}

    private:
        const substrings* parent_;
        typename sast_type::const_iterator i_;
        int ii_;
    };

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
            : parent_(0), i_(), ii_(-1)
        {}

        substring_iterator(const substrings* parent, typename sast_type::const_iterator i, int ii)
            : parent_(parent), i_(i), ii_(ii)
        {}

        void next_branching() {
            ++i_;
            ii_ = 0;
        }

        void prev_branching() {
            --i_;
            auto j = i_.parent();
            ii_ = (i_->length() - j->length()) - 1;
        }

    private:
        friend class boost::iterator_core_access;

        void increment() {
            auto j = i_.parent();
            if (ii_ + 1 < i_->length() - j->length()) {
                ++ii_;
            }
            else {
                next_branching();
            }
        }

        void decrement() {
            if (ii_ - 1 >= 0) {
                --ii_;
            }
            else {
                prev_branching();
            }
        }

        void advance(int n) {
            if (n >= 0) {
                auto j = i_.parent();
                while (ii_ + n >= i_->length() - j->length()) {
                    n -= (i_->length() - j->length()) - ii_;
                    next_branching();
                    j = i_.parent();
                }
                ii_ += n;
            }
            else {
                while (ii_ + n < 0) {
                    n += ii_ + 1;
                    prev_branching();
                }
                ii_ += n;
            }
        }

        int distance_to(const substring_iterator<Value>& other) const {
            int d = 0;
            if (other.i_ - this->i_ > 0 || (other.i_ == this->i_ && other.ii_ >= this->ii_)) {
                auto i = this->i_;
                int ii = this->ii_;
                while (other.i_ - i > 0) {
                    auto j = i.parent();
                    d += (i->length() - j->length()) - ii;
                    ++i;
                    ii = 0;
                }
                d += other.ii_ - ii_;
            }
            else {
                auto i = other.i_;
                int ii = other.ii_;
                while (this->i_ - i > 0) {
                    auto j = i.parent();
                    d += (i->length() - j->length()) - ii;
                    ++i;
                    ii = 0;
                }
                d += other.ii_ - ii_;
                d = -d;
            }
            return d;
        }

        bool equal(const substring_iterator<Value>& other) const {
            return this->parent_ == other.parent_
                && this->i_ == other.i_
                && this->ii_ == other.ii_;
        }

        Value dereference() const {
            return substr(parent_, i_, ii_);
        }

        const substrings* parent_;
        typename sast_type::const_iterator i_;
        int ii_;
    };

public:
    typedef substring_iterator<substr> iterator;
    typedef substring_iterator<const substr> const_iterator;

    substrings(const RandomAccessRange& input, const size_t alphabet_size)
        : input_(input),
          sast_(input, alphabet_size)
    {}

    iterator begin() { return iterator(this, sast_.begin(), 0); }
    iterator end()   { return iterator(this, sast_.begin() + (sast_.size() - 1), 0); }
    const_iterator begin() const { return const_iterator(this, sast_.begin(), 0); }
    const_iterator end()   const { return const_iterator(this, sast_.begin() + (sast_.size() - 1), 0); }

private:
    double strict_purity(typename sast_type::const_iterator n, const int ii) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列の末尾をii文字削ったsubstrを扱う。
        // iiが満たさなければならない条件: 0 <= ii < d[i] - d[node_to_parent_node[i]]
        const auto freq_substr = n->frequency();
        const auto len_substr  = n->length() - ii;

        // substrと同じ出現回数のsub-substrを数える。
        uint64_t count = 0;
        {
            // substrの末尾を0文字以上削って得られるsub-substrについて考える。
            // ノードiに対応する部分文字列をsubstr[i]とすると、substr[i]の末尾
            // を削って得られる部分文字列の内で、substr[i]と同じ頻度をもつもの
            // はsuffix tree上ではノードiにまとめられている
            // （分岐が無い <=> 頻度が同じ）。

            // ノードiの親ノードjを見つける。
            const auto m = n.parent();

            // substrの末尾を0文字以上削って得られるsub-substrの内で、出現
            // 回数がsubstrと同じものの数はd[i] - d[j]である。
            count += n->length() - m->length() - ii;
        }
        for (int j = 1; j < len_substr; ++j) {
            // substrの先頭をj文字削ったsub-substrを考える。
            const auto len_subsubstr = len_substr - j;

            // sub-substrに対応するノードを見つける。
            // {内部,葉}ノードに対応する部分文字列から先頭の1文字を削って得られ
            // る文字列には、必ず対応する{内部,葉}ノードが存在する。
            // ii ＞ 0 の場合、d[k] == len_subsubstr を満たす「ノード」が必ずし
            // も存在するとは限らない。よって、
            // d[node_to_parent_node[k]] < len_subsubstr <= d[k] を満たす k を
            // 見つける。
            auto m = n;
            for (int l = 0; l < j; ++l) {
                m = m.suffix();
            }
            while (m.parent()->length() >= len_subsubstr) {
                m = m.parent();
            }
            const auto kk = m->length() - len_subsubstr;  // ノードiでii文字削ると、ノードkではkk文字削ったことに相当する。

            // sub-substrの出現回数をチェックする。
            const auto freq_subsubstr = m->frequency();
            if (freq_subsubstr == freq_substr) {
                // ノードkの親ノードmを見つける。
                const auto mp = m.parent();

                // sub-substrの末尾を0文字以上削って得られる
                // sub-sub-substrの内で、出現回数がsub-substrと同じもの
                // の数はd[k] - d[m]である。
                count += m->length() - mp->length() - kk;
            }
        }

        // strict purity of substr
        const uint64_t num_subsubstrs = static_cast<uint64_t>(1 + len_substr) * len_substr / 2;
        const double spurity = static_cast<double>(count) / num_subsubstrs;

        return spurity;
    }

    double loose_purity(typename sast_type::const_iterator n, const int ii) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列の末尾をii文字削ったsubstrを扱う。
        // iiが満たさなければならない条件: 0 <= ii < d[i] - d[node_to_parent_node[i]]
        const auto freq_substr = n->frequency();
        const auto len_substr  = n->length() - ii;

        double support = 0;
        {
            // substrの末尾を0文字以上削って得られるsub-substrについて考える。
            for (auto m = n, p = n.parent(); m->length() > 0; m = p, p = p.parent()) {
                const auto num_subsubstrs_of_same_frequency = m->length() - p->length() - (m == n ? ii : 0);
                const auto freq_subsubstr = m->frequency();
                const double sup = static_cast<double>(freq_substr) / freq_subsubstr;
                support += num_subsubstrs_of_same_frequency * sup;
            }
        }
        for (int j = 1; j < len_substr; ++j) {
            // substrの先頭をj文字削ったsub-substrを考える。
            const auto len_subsubstr = len_substr - j;

            // sub-substrに対応するノードを見つける。
            // ii ＞ 0 の場合、d[k] == len_subsubstr を満たす「ノード」が必ずし
            // も存在するとは限らない。よって、
            // d[node_to_parent_node[k]] < len_subsubstr <= d[k] を満たす k を
            // 見つける。
            auto m = n;
            for (int l = 0; l < j; ++l) {
                m = m.suffix();
            }
            while (m.parent()->length() >= len_subsubstr) {
                m = m.parent();
            }
            const auto kk = m->length() - len_subsubstr;  // ノードiでii文字削ると、ノードkではkk文字削ったことに相当する。

            // sub-substrの末尾を0文字以上削って得られるsub-substrについて考える。
            for (auto o = m, p = m.parent(); o->length() > 0; o = p, p = p.parent()) {
                const auto num_subsubstrs_of_same_frequency = o->length() - p->length() - (o == m ? kk : 0);
                const auto freq_subsubstr = o->frequency();
                const double sup = static_cast<double>(freq_substr) / freq_subsubstr;
                support += num_subsubstrs_of_same_frequency * sup;
            }
        }

        // loose purity of substr
        const uint64_t num_subsubstrs = static_cast<uint64_t>(1 + len_substr) * len_substr / 2;
        const double lpurity = support / num_subsubstrs;

        return lpurity;
    }

    std::map<char_type, int> left_extensions(typename sast_type::const_iterator n, const int ii) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列の末尾をii文字削ったsubstrを扱う。
        // iiが満たさなければならない条件: 0 <= ii < d[i] - d[node_to_parent_node[i]]

        std::map<char_type, int> char_dist;
        for (auto pos : n->allpos()) {
            const auto& c = input_[pos - 1];
            char_dist[c] += 1;
        }

        return char_dist;
    }

    std::map<char_type, int> right_extensions(typename sast_type::const_iterator n, const int ii) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列の末尾をii文字削ったsubstrを扱う。
        // iiが満たさなければならない条件: 0 <= ii < d[i] - d[node_to_parent_node[i]]
        const auto len_substr  = n->length() - ii;

        std::map<char_type, int> char_dist;
        for (auto pos : n->allpos()) {
            const auto& c = input_[pos + len_substr];
            char_dist[c] += 1;
        }

        return char_dist;
    }

    double left_universality(typename sast_type::const_iterator n, const int ii) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列の末尾をii文字削ったsubstrを扱う。
        // iiが満たさなければならない条件: 0 <= ii < d[i] - d[node_to_parent_node[i]]
        const auto freq_substr = n->frequency();

        std::map<char_type, int> char_dist = left_extensions(n, ii);

        double e = 0;
        for (const auto& kv : char_dist) {
            const double p = static_cast<double>(kv.second) / freq_substr;
            e += -p * std::log(p);
        }
        const double u = 1 - std::exp(-e);

        return u;
    }

    double right_universality(typename sast_type::const_iterator n, const int ii) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列の末尾をii文字削ったsubstrを扱う。
        // iiが満たさなければならない条件: 0 <= ii < d[i] - d[node_to_parent_node[i]]
        const auto freq_substr = n->frequency();

        std::map<char_type, int> char_dist = right_extensions(n, ii);

        double e = 0;
        for (const auto& kv : char_dist) {
            const double p = static_cast<double>(kv.second) / freq_substr;
            e += -p * std::log(p);
        }
        const double u = 1 - std::exp(-e);

        return u;
    }

    const RandomAccessRange& input_;
    sast_type sast_;
};


#endif  /* SUBSTRINGS_HPP */
