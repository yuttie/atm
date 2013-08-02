#ifndef SUBSTRINGS_FROM_LONGEST_H
#define SUBSTRINGS_FROM_LONGEST_H

#include <cmath>
#include <map>
#include <vector>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/range/size.hpp>
#include <boost/range/sub_range.hpp>
#include <boost/range/value_type.hpp>
#include <boost/utility.hpp>
#include "sast/sast.hpp"


namespace atm {

template <class RandomAccessRange, class Index>
struct substrings_from_longest {
    using range_type = RandomAccessRange;
    using char_type = typename boost::range_value<RandomAccessRange>::type;
    using index_type = Index;

protected:
    using sast_type = sast::sast<RandomAccessRange, index_type>;

public:
    struct substr {
        using iterator       = typename boost::range_iterator<RandomAccessRange>::type;
        using const_iterator = typename boost::range_const_iterator<RandomAccessRange>::type;

        index_type pos()           const { return parent_->pos(i_, j_); }
        boost::sub_range<const std::vector<index_type>> allpos() const { return parent_->allpos(i_, j_); }
        index_type length()        const { return j_ - i_; }
        index_type frequency()     const { return parent_->frequency(i_, j_); }
        double     spurity()       const { return parent_->strict_purity(i_, j_); }
        double     lpurity()       const { return parent_->loose_purity(i_, j_); }
        double     luniversality() const { return parent_->left_universality(i_, j_); }
        double     runiversality() const { return parent_->right_universality(i_, j_); }

        iterator begin() { return boost::begin(parent_->input_) + i_; }
        iterator end()   { return boost::begin(parent_->input_) + j_; }
        const_iterator begin() const { return boost::begin(parent_->input_) + i_; }
        const_iterator end()   const { return boost::begin(parent_->input_) + j_; }

        substr(const substrings_from_longest* parent, int i, int j)
            : parent_(parent), i_(i), j_(j)
        {}

    private:
        const substrings_from_longest* parent_;
        int i_;
        int j_;
    };

private:
    template <class> struct substring_iterator;

public:
    using iterator       = substring_iterator<substr>;
    using const_iterator = substring_iterator<const substr>;

    substrings_from_longest(const sast_type& sast)
        : sast_(sast),
          input_(sast_.input()),
          finder_(sast::make_positional_finder(sast_))
    {}

    iterator begin() { return iterator(this, 0, boost::size(input_)); }
    iterator end()   { return iterator(this, 0, 0); }
    const_iterator begin() const { return const_iterator(this, 0, boost::size(input_)); }
    const_iterator end()   const { return const_iterator(this, 0, 0); }

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
            : parent_(0), i_(-1), j_(-1)
        {}

        substring_iterator(const substrings_from_longest* parent, int i, int j)
            : parent_(parent), i_(i), j_(j)
        {}

    private:
        friend class boost::iterator_core_access;

        void increment() {
            if (static_cast<std::size_t>(j_) == boost::size(parent_->input_)) {
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
                i_ = boost::size(parent_->input_) - width;
                j_ = boost::size(parent_->input_);
            }
            else {
                --i_;
                --j_;
            }
        }

        void advance(int n) {
            if (n >= 0) {
                auto width = j_ - i_;
                while (static_cast<std::size_t>(n) > boost::size(parent_->input_) - j_) {
                    n -= (boost::size(parent_->input_) - j_) + 1;
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
                    i_ = boost::size(parent_->input_) - width;
                    j_ = boost::size(parent_->input_);
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
                    d += (boost::size(parent_->input_) - j) + 1;
                    --width;
                    j = 0 + width;
                }
                d += other.j_ - j;
            }
            else {
                int j = other.j_;
                while (owidth > width) {
                    d += (boost::size(parent_->input_) - j) + 1;
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
            return substr(parent_, i_, j_);
        }

        const substrings_from_longest* parent_;
        int i_;
        int j_;
    };

protected:
    index_type pos(const int i, const int j) const {
        const auto len_substr = j - i;
        // substrに対応する内部ノードを見つける。
        typename sast_type::const_iterator n = finder_.find(i, j);
        if (n->length() >= len_substr) {
            // substrは2回以上出現しており、対応する内部ノードが存在する。
            return n->pos();
        }
        else {
            // substrは1回しか出現しておらず、対応する内部ノードが存在しない。
            return i;
        }
    }

    boost::sub_range<const std::vector<index_type>> allpos(const int i, const int j) const {
        const auto len_substr = j - i;
        // substrに対応する内部ノードを見つける。
        typename sast_type::const_iterator n = finder_.find(i, j);
        if (n->length() >= len_substr) {
            // substrは2回以上出現しており、対応する内部ノードが存在する。
            return n->allpos();
        }
        else {
            // substrは1回しか出現しておらず、対応する内部ノードが存在しない。
            auto found = boost::find(n->allpos(), i);
            return boost::make_iterator_range(found, boost::next(found));
        }
    }

    index_type frequency(const int i, const int j) const {
        const auto len_substr = j - i;
        // substrに対応する内部ノードを見つける。
        typename sast_type::const_iterator n = finder_.find(i, j);
        if (n->length() >= len_substr) {
            // substrは2回以上出現しており、対応する内部ノードが存在する。
            return n->frequency();
        }
        else {
            // substrは1回しか出現しておらず、対応する内部ノードが存在しない。
            return 1;
        }
    }

    double strict_purity(const int i, const int j) const {
        const auto len_substr = j - i;
        // substrに対応する内部ノードを見つける。
        typename sast_type::const_iterator n = finder_.find(i, j);
        if (n->length() >= len_substr) {
            // substrは2回以上出現しており、対応する内部ノードが存在する。
            const auto kk = n->length() - len_substr;  // ノードkではkk文字削ったことに相当する。

            return strict_purity_(n, kk);
        }
        else {
            // substrは1回しか出現しておらず、対応する内部ノードが存在しない。
            // substrと同じ出現回数のsub-substrを数える。
            uint64_t count = 0;
            {
                // substrの末尾を0文字以上削って得られるsub-substrについて考える。
                count += len_substr - n->length();
            }
            for (int l = 1; l < len_substr; ++l) {
                // substrの先頭をl文字削ったsub-substrを考える。
                const auto len_subsubstr = len_substr - l;

                // sub-substrに対応するノードを見つける。
                // d[node_to_parent_node[m]] < len_subsubstr <= d[m] を満たす m を
                // 見つける。
                typename sast_type::const_iterator m = finder_.find(i + l, j);
                if (m->length() < len_subsubstr) {
                    count += len_subsubstr - m->length();
                }
            }

            // strict purity of substr
            const uint64_t num_subsubstrs = static_cast<uint64_t>(1 + len_substr) * len_substr / 2;
            const double spurity = static_cast<double>(count) / num_subsubstrs;

            return spurity;
        }
    }

    double strict_purity_(typename sast_type::const_iterator n, const int ii) const {
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

    double loose_purity(const int i, const int j) const {
        const auto len_substr = j - i;
        // substrに対応する内部ノードを見つける。
        typename sast_type::const_iterator n = finder_.find(i, j);
        if (n->length() >= len_substr) {
            // substrは2回以上出現しており、対応する内部ノードが存在する。
            const auto kk = n->length() - len_substr;  // ノードkではkk文字削ったことに相当する。

            return loose_purity_(n, kk);
        }
        else {
            // substrは1回しか出現しておらず、対応する内部ノードが存在しない。
            // substrと同じ出現回数のsub-substrを数える。
            double support = 0;
            {
                // substrの末尾を0文字以上削って得られるsub-substrについて考える。
                support += len_substr - n->length();
                for (auto m = n, p = n.parent(); m->length() > 0; m = p, p = p.parent()) {
                    const auto num_subsubstrs_of_same_frequency = m->length() - p->length();
                    const auto freq_subsubstr = m->frequency();
                    const double sup = 1.0 / freq_subsubstr;
                    support += num_subsubstrs_of_same_frequency * sup;
                }
            }
            for (int l = 1; l < len_substr; ++l) {
                // substrの先頭をl文字削ったsub-substrを考える。
                const auto len_subsubstr = len_substr - l;

                // sub-substrに対応するノードを見つける。
                // d[node_to_parent_node[m]] < len_subsubstr <= d[m] を満たす m を
                // 見つける。
                typename sast_type::const_iterator m = finder_.find(i + l, j);
                if (m->length() >= len_subsubstr) {
                    const auto mm = m->length() - len_subsubstr;  // ノードmではmm文字削ったことに相当する。

                    // sub-substrの末尾を0文字以上削って得られるsub-substrについて考える。
                    for (auto o = m, p = m.parent(); o->length() > 0; o = p, p = p.parent()) {
                        const auto num_subsubstrs_of_same_frequency = o->length() - p->length() - (o == m ? mm : 0);
                        const auto freq_subsubstr = o->frequency();
                        const double sup = 1.0 / freq_subsubstr;
                        support += num_subsubstrs_of_same_frequency * sup;
                    }
                }
                else {
                    support += len_subsubstr - m->length();
                    for (auto o = m, p = m.parent(); o->length() > 0; o = p, p = p.parent()) {
                        const auto num_subsubstrs_of_same_frequency = o->length() - p->length();
                        const auto freq_subsubstr = o->frequency();
                        const double sup = 1.0 / freq_subsubstr;
                        support += num_subsubstrs_of_same_frequency * sup;
                    }
                }
            }

            // loose purity of substr
            const uint64_t num_subsubstrs = static_cast<uint64_t>(1 + len_substr) * len_substr / 2;
            const double lpurity = support / num_subsubstrs;

            return lpurity;
        }
    }

    double loose_purity_(typename sast_type::const_iterator n, const int ii) const {
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

    std::map<char_type, int> left_extensions(const int i, const int j) const {
        std::map<char_type, int> char_dist;
        for (const auto pos : allpos(i, j)) {
            const auto& c = boost::begin(input_)[pos - 1];
            char_dist[c] += 1;
        }

        return char_dist;
    }

    std::map<char_type, int> right_extensions(const int i, const int j) const {
        const auto len_substr = j - i;

        std::map<char_type, int> char_dist;
        for (const auto pos : allpos(i, j)) {
            const auto& c = boost::begin(input_)[pos + len_substr];
            char_dist[c] += 1;
        }

        return char_dist;
    }

    double left_universality(const int i, const int j) const {
        const auto freq_substr = frequency(i, j);

        std::map<char_type, int> char_dist = left_extensions(i, j);

        double e = 0;
        for (const auto& kv : char_dist) {
            const double p = static_cast<double>(kv.second) / freq_substr;
            e += -p * std::log(p);
        }
        const double u = 1 - std::exp(-e);

        return u;
    }

    double right_universality(const int i, const int j) const {
        const auto freq_substr = frequency(i, j);

        std::map<char_type, int> char_dist = right_extensions(i, j);

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
    sast::positional_finder<RandomAccessRange, index_type> finder_;
};

}  // namespace atm


#endif  /* SUBSTRINGS_FROM_LONGEST_H */
