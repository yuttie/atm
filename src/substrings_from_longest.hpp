#ifndef SUBSTRINGS_FROM_LONGEST_H
#define SUBSTRINGS_FROM_LONGEST_H

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
#include "esa.hxx"


template <class RandomAccessRange, class Index>
struct SubstringsFromLongest {
protected:
    using char_type = typename boost::range_value<RandomAccessRange>::type;

public:
    typedef Index index_type;
    struct substr {
        using iterator       = typename boost::range_iterator<RandomAccessRange>::type;
        using const_iterator = typename boost::range_const_iterator<RandomAccessRange>::type;

        index_type              pos()           const { return i_; }
        std::vector<index_type> allpos()        const { return parent_->allpos(i_, j_); }
        index_type              length()        const { return j_ - i_; }
        index_type              frequency()     const { return parent_->frequency(i_, j_); }
        double                  spurity()       const { return parent_->strict_purity(i_, j_); }
        double                  lpurity()       const { return parent_->loose_purity(i_, j_); }
        double                  luniversality() const { return parent_->left_universality(i_, j_); }
        double                  runiversality() const { return parent_->right_universality(i_, j_); }

        iterator begin() {
            return boost::begin(parent_->input_) + i_;
        }

        iterator end() {
            return boost::begin(parent_->input_) + j_;
        }

        const_iterator begin() const {
            return boost::begin(parent_->input_) + i_;
        }

        const_iterator end() const {
            return boost::begin(parent_->input_) + j_;
        }

        substr(const SubstringsFromLongest* parent, int i, int j)
            : parent_(parent), i_(i), j_(j)
        {}

    private:
        const SubstringsFromLongest* parent_;
        int i_;
        int j_;
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
            : parent_(0), i_(-1), j_(-1)
        {}

        substring_iterator(const SubstringsFromLongest* parent, int i, int j)
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

        const SubstringsFromLongest* parent_;
        int i_;
        int j_;
    };

public:
    typedef substring_iterator<substr> iterator;
    typedef substring_iterator<const substr> const_iterator;

    SubstringsFromLongest(const RandomAccessRange& input, const size_t alphabet_size)
        : input_(input),
          sa_(input.size()),
          l_(input.size()),
          r_(input.size()),
          d_(input.size()),
          suffix_to_parent_node_(input.size(), -1),
          node_to_parent_node_()
    {
        // suffix array
        int err = esaxx(boost::begin(input_),
                        sa_.begin(),
                        l_.begin(), r_.begin(), d_.begin(),
                        static_cast<index_type>(boost::size(input_)),
                        static_cast<index_type>(alphabet_size),
                        num_nodes_);
        if (err) throw std::runtime_error("saisxx failed to construct a suffix array.");

        // suffix_to_parent_node[k]: 接尾辞input[k..$]に対応する葉ノードの、親ノードのpost-order順の番号。
        // post-order巡回により、直接の親が最初に値を設定する（最初かどうかは-1かどうかで判定する）。
        for (int i = 0; i < num_nodes_; ++i) {
            // ノードi直下の全ての葉ノードjについて、接尾辞input[k..$]からノードiへのリンクを張る
            for (int j = l_[i]; j < r_[i]; ++j) {
                const auto k = sa_[j];
                if (suffix_to_parent_node_[k] < 0) {
                    suffix_to_parent_node_[k] = i;
                }
            }
        }

        // node_to_parent_node[i]: ノードiの親ノードの番号（post-order）。
        node_to_parent_node_.resize(num_nodes_ - 1);
        std::stack<index_type> stk;
        stk.push(num_nodes_ - 1);
        for (int i = num_nodes_ - 2; i >= 0; --i) {
            while (!(l_[stk.top()] <= l_[i] && r_[i] <= r_[stk.top()])) {
                stk.pop();
            }
            node_to_parent_node_[i] = stk.top();
            stk.push(i);
        }
    }

    iterator begin() {
        return iterator(this, 0, boost::size(input_));
    }

    iterator end() {
        return iterator(this, 0, 0);
    }

    const_iterator begin() const {
        return const_iterator(this, 0, boost::size(input_));
    }

    const_iterator end() const {
        return const_iterator(this, 0, 0);
    }

protected:
    std::vector<index_type> allpos(const int i, const int j) const {
        const auto len_substr = j - i;
        const auto pos_substr = i;
        // substrに対応する内部ノードを見つける。
        auto k = suffix_to_parent_node_[pos_substr];  // 接尾辞input[pos_substr..$]に対応する葉ノードの親ノード
        if (d_[k] >= len_substr) {
            // substrは2回以上出現しており、対応する内部ノードが存在する。
            // d[node_to_parent_node[k]] < len_substr <= d[k] を満たす k を
            // 見つける。
            while (d_[node_to_parent_node_[k]] >= len_substr) k = node_to_parent_node_[k];

            return std::vector<index_type>(sa_.begin() + l_[k], sa_.begin() + r_[k]);
        }
        else {
            // substrは1回しか出現しておらず、対応する内部ノードが存在しない。
            return {i};
        }
    }

    index_type frequency(const int i, const int j) const {
        const auto len_substr = j - i;
        const auto pos_substr = i;
        // substrに対応する内部ノードを見つける。
        auto k = suffix_to_parent_node_[pos_substr];  // 接尾辞input[pos_substr..$]に対応する葉ノードの親ノード
        if (d_[k] >= len_substr) {
            // substrは2回以上出現しており、対応する内部ノードが存在する。
            // d[node_to_parent_node[k]] < len_substr <= d[k] を満たす k を
            // 見つける。
            while (d_[node_to_parent_node_[k]] >= len_substr) k = node_to_parent_node_[k];

            return r_[k] - l_[k];
        }
        else {
            // substrは1回しか出現しておらず、対応する内部ノードが存在しない。
            return 1;
        }
    }

    double strict_purity(const int i, const int j) const {
        const auto len_substr = j - i;
        const auto pos_substr = i;
        // substrに対応する内部ノードを見つける。
        auto k = suffix_to_parent_node_[pos_substr];  // 接尾辞input[pos_substr..$]に対応する葉ノードの親ノード
        if (d_[k] >= len_substr) {
            // substrは2回以上出現しており、対応する内部ノードが存在する。
            // d[node_to_parent_node[k]] < len_substr <= d[k] を満たす k を
            // 見つける。
            while (d_[node_to_parent_node_[k]] >= len_substr) k = node_to_parent_node_[k];
            const auto kk = d_[k] - len_substr;  // ノードkではkk文字削ったことに相当する。

            return strict_purity_(k, kk);
        }
        else {
            // substrは1回しか出現しておらず、対応する内部ノードが存在しない。
            // substrと同じ出現回数のsub-substrを数える。
            uint64_t count = 0;
            {
                // substrの末尾を0文字以上削って得られるsub-substrについて考える。
                count += len_substr - d_[k];
            }
            for (int l = 1; l < len_substr; ++l) {
                // substrの先頭をl文字削ったsub-substrを考える。
                const auto len_subsubstr = len_substr - l;

                // sub-substrに対応するノードを見つける。
                // d[node_to_parent_node[m]] < len_subsubstr <= d[m] を満たす m を
                // 見つける。
                auto m = suffix_to_parent_node_[pos_substr + l];  // 接尾辞input[(pos_substr + l)..$]に対応する葉ノードの親ノード
                if (d_[m] < len_subsubstr) {
                    count += len_subsubstr - d_[m];
                }
            }

            // strict purity of substr
            const uint64_t num_subsubstrs = static_cast<uint64_t>(1 + len_substr) * len_substr / 2;
            const double spurity = static_cast<double>(count) / num_subsubstrs;

            return spurity;
        }
    }

    double strict_purity_(const int i, const int ii) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列の末尾をii文字削ったsubstrを扱う。
        // iiが満たさなければならない条件: 0 <= ii < d[i] - d[node_to_parent_node[i]]
        const auto freq_substr = r_[i] - l_[i];
        const auto len_substr  = d_[i] - ii;
        const auto pos_substr  = sa_[l_[i]];

        // substrと同じ出現回数のsub-substrを数える。
        uint64_t count = 0;
        {
            // substrの末尾を0文字以上削って得られるsub-substrについて考える。
            // ノードiに対応する部分文字列をsubstr[i]とすると、substr[i]の末尾
            // を削って得られる部分文字列の内で、substr[i]と同じ頻度をもつもの
            // はsuffix tree上ではノードiにまとめられている
            // （分岐が無い <=> 頻度が同じ）。

            // ノードiの親ノードjを見つける。
            const auto j = node_to_parent_node_[i];

            // substrの末尾を0文字以上削って得られるsub-substrの内で、出現
            // 回数がsubstrと同じものの数はd[i] - d[j]である。
            count += d_[i] - d_[j] - ii;
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
            auto k = suffix_to_parent_node_[pos_substr + j];  // 接尾辞input[(pos_substr + j)..$]に対応する葉ノードの親ノード
            while (d_[node_to_parent_node_[k]] >= len_subsubstr) k = node_to_parent_node_[k];
            const auto kk = d_[k] - len_subsubstr;  // ノードiでii文字削ると、ノードkではkk文字削ったことに相当する。

            // sub-substrの出現回数をチェックする。
            const auto freq_subsubstr = r_[k] - l_[k];
            if (freq_subsubstr == freq_substr) {
                // ノードkの親ノードmを見つける。
                const auto m = node_to_parent_node_[k];

                // sub-substrの末尾を0文字以上削って得られる
                // sub-sub-substrの内で、出現回数がsub-substrと同じもの
                // の数はd[k] - d[m]である。
                count += d_[k] - d_[m] - kk;
            }
        }

        // strict purity of substr
        const uint64_t num_subsubstrs = static_cast<uint64_t>(1 + len_substr) * len_substr / 2;
        const double spurity = static_cast<double>(count) / num_subsubstrs;

        return spurity;
    }

    double loose_purity(const int i, const int j) const {
        const auto len_substr = j - i;
        const auto pos_substr = i;
        // substrに対応する内部ノードを見つける。
        auto k = suffix_to_parent_node_[pos_substr];  // 接尾辞input[pos_substr..$]に対応する葉ノードの親ノード
        if (d_[k] >= len_substr) {
            // substrは2回以上出現しており、対応する内部ノードが存在する。
            // d[node_to_parent_node[k]] < len_substr <= d[k] を満たす k を
            // 見つける。
            while (d_[node_to_parent_node_[k]] >= len_substr) k = node_to_parent_node_[k];
            const auto kk = d_[k] - len_substr;  // ノードkではkk文字削ったことに相当する。

            return loose_purity_(k, kk);
        }
        else {
            // substrは1回しか出現しておらず、対応する内部ノードが存在しない。
            // substrと同じ出現回数のsub-substrを数える。
            double support = 0;
            {
                // substrの末尾を0文字以上削って得られるsub-substrについて考える。
                support += len_substr - d_[k];
                for (index_type l = k, m = node_to_parent_node_[k]; d_[l] > 0; l = m, m = node_to_parent_node_[m]) {
                    const auto num_subsubstrs_of_same_frequency = d_[l] - d_[m];
                    const auto freq_subsubstr = r_[l] - l_[l];
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
                auto m = suffix_to_parent_node_[pos_substr + l];  // 接尾辞input[(pos_substr + l)..$]に対応する葉ノードの親ノード
                if (d_[m] >= len_subsubstr) {
                    while (d_[node_to_parent_node_[m]] >= len_subsubstr) m = node_to_parent_node_[m];
                    const auto mm = d_[m] - len_subsubstr;  // ノードmではmm文字削ったことに相当する。

                    // sub-substrの末尾を0文字以上削って得られるsub-substrについて考える。
                    for (index_type n = m, o = node_to_parent_node_[m]; d_[n] > 0; n = o, o = node_to_parent_node_[o]) {
                        const auto num_subsubstrs_of_same_frequency = d_[n] - d_[o] - (n == m ? mm : 0);
                        const auto freq_subsubstr = r_[n] - l_[n];
                        const double sup = 1.0 / freq_subsubstr;
                        support += num_subsubstrs_of_same_frequency * sup;
                    }
                }
                else {
                    support += len_subsubstr - d_[m];
                    for (index_type n = m, o = node_to_parent_node_[m]; d_[n] > 0; n = o, o = node_to_parent_node_[o]) {
                        const auto num_subsubstrs_of_same_frequency = d_[n] - d_[o];
                        const auto freq_subsubstr = r_[n] - l_[n];
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

    double loose_purity_(const int i, const int ii) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列の末尾をii文字削ったsubstrを扱う。
        // iiが満たさなければならない条件: 0 <= ii < d[i] - d[node_to_parent_node[i]]
        const auto freq_substr = r_[i] - l_[i];
        const auto len_substr  = d_[i] - ii;
        const auto pos_substr  = sa_[l_[i]];

        double support = 0;
        {
            // substrの末尾を0文字以上削って得られるsub-substrについて考える。
            for (index_type j = i, k = node_to_parent_node_[i]; d_[j] > 0; j = k, k = node_to_parent_node_[k]) {
                const auto num_subsubstrs_of_same_frequency = d_[j] - d_[k] - (j == i ? ii : 0);
                const auto freq_subsubstr = r_[j] - l_[j];
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
            auto k = suffix_to_parent_node_[pos_substr + j];  // 接尾辞input[(pos_substr + j)..$]に対応する葉ノードの親ノード
            while (d_[node_to_parent_node_[k]] >= len_subsubstr) k = node_to_parent_node_[k];
            const auto kk = d_[k] - len_subsubstr;  // ノードiでii文字削ると、ノードkではkk文字削ったことに相当する。

            // sub-substrの末尾を0文字以上削って得られるsub-substrについて考える。
            for (index_type m = k, n = node_to_parent_node_[k]; d_[m] > 0; m = n, n = node_to_parent_node_[n]) {
                const auto num_subsubstrs_of_same_frequency = d_[m] - d_[n] - (m == k ? kk : 0);
                const auto freq_subsubstr = r_[m] - l_[m];
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
            const auto& c = input_[pos - 1];
            char_dist[c] += 1;
        }

        return char_dist;
    }

    std::map<char_type, int> right_extensions(const int i, const int j) const {
        const auto len_substr = j - i;

        std::map<char_type, int> char_dist;
        for (const auto pos : allpos(i, j)) {
            const auto& c = input_[pos + len_substr];
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

    const RandomAccessRange& input_;
    std::vector<index_type> sa_;
    std::vector<index_type>  l_;
    std::vector<index_type>  r_;
    std::vector<index_type>  d_;
    index_type  num_nodes_;
    std::vector<index_type> suffix_to_parent_node_;
    std::vector<index_type> node_to_parent_node_;
};


#endif  /* SUBSTRINGS_FROM_LONGEST_H */
