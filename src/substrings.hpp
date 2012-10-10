#ifndef SUBSTRINGS_HPP
#define SUBSTRINGS_HPP

#include <stack>
#include <stdexcept>
#include <vector>
#include <boost/iterator/iterator_facade.hpp>
#include "esa.hxx"


template <class Char, class Index>
struct Substrings {
    typedef Index index_type;
    struct substr {
        typedef typename std::vector<Char>::const_iterator iterator;
        typedef typename std::vector<Char>::const_iterator const_iterator;

        index_type pos()       const { return parent_->sa_[parent_->l_[i_]]; }
        index_type length()    const { return parent_->d_[i_] - ii_; }
        index_type frequency() const { return parent_->r_[i_] - parent_->l_[i_]; }
        double     spurity()   const { return parent_->strict_purity(i_, ii_); }
        double     lpurity()   const { return parent_->loose_purity(i_, ii_); }

        iterator begin() const {
            return parent_->input_.begin() + pos();
        }

        iterator end() const {
            return parent_->input_.begin() + pos() + length();
        }

        substr(const Substrings* parent, int i, int ii)
            : parent_(parent), i_(i), ii_(ii)
        {}

    private:
        const Substrings* parent_;
        int i_;
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
            : parent_(0), i_(-1), ii_(-1)
        {}

        substring_iterator(const Substrings* parent, int i, int ii)
            : parent_(parent), i_(i), ii_(ii)
        {}

    private:
        friend class boost::iterator_core_access;

        void increment() {
            const auto j = parent_->node_to_parent_node_[i_];
            if (ii_ + 1 < parent_->d_[i_] - parent_->d_[j]) {
                ++ii_;
            }
            else {
                ++i_;
                ii_ = 0;
            }
        }

        void decrement() {
            if (ii_ - 1 >= 0) {
                --ii_;
            }
            else {
                --i_;
                const auto j = parent_->node_to_parent_node_[i_];
                ii_ = (parent_->d_[i_] - parent_->d_[j]) - 1;
            }
        }

        void advance(int n) {
            if (n >= 0) {
                auto j = parent_->node_to_parent_node_[i_];
                while (ii_ + n >= parent_->d_[i_] - parent_->d_[j]) {
                    n -= (parent_->d_[i_] - parent_->d_[j]) - ii_;
                    ++i_;
                    j = parent_->node_to_parent_node_[i_];
                    ii_ = 0;
                }
                ii_ += n;
            }
            else {
                while (ii_ + n < 0) {
                    n += ii_ + 1;
                    --i_;
                    const auto j = parent_->node_to_parent_node_[i_];
                    ii_ = (parent_->d_[i_] - parent_->d_[j]) - 1;
                }
                ii_ += n;
            }
        }

        int distance_to(const substring_iterator<Value>& other) const {
            int d = 0;
            if (other.i_ > this->i_ || other.i_ == this->i_ && other.ii_ >= this->ii_) {
                int i = this->i_;
                int ii = this->ii_;
                while (other.i_ > i) {
                    const auto j = parent_->node_to_parent_node_[i];
                    d += (parent_->d_[i] - parent_->d_[j]) - ii;
                    ++i;
                    ii = 0;
                }
                d += other.ii_ - ii_;
            }
            else {
                int i = other.i_;
                int ii = other.ii_;
                while (this->i_ > i) {
                    const auto j = parent_->node_to_parent_node_[i];
                    d += (parent_->d_[i] - parent_->d_[j]) - ii;
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

        const Substrings* parent_;
        int i_;
        int ii_;
    };

public:
    typedef substring_iterator<substr> iterator;
    typedef substring_iterator<const substr> const_iterator;

    Substrings(const std::vector<Char>& input, const size_t alphabet_size)
        : input_(input),
          sa_(input.size()),
          l_(input.size()),
          r_(input.size()),
          d_(input.size()),
          suffix_to_parent_node_(input.size(), -1),
          node_to_parent_node_()
    {
        // suffix array
        int err = esaxx(input_.begin(),
                        sa_.begin(),
                        l_.begin(), r_.begin(), d_.begin(),
                        static_cast<index_type>(input_.size()),
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
        node_to_parent_node_.resize(num_nodes_);
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

    iterator begin() const {
        return iterator(this, 0, 0);
    }

    iterator end() const {
        return iterator(this, num_nodes_ - 1, 0);
    }

private:
    double strict_purity(const int i, const int ii) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列の末尾をii文字削ったsubstrを扱う。
        // iiが満たさなければならない条件: 0 <= ii < d[i] - d[node_to_parent_node[i]]
        const auto freq_substr = r_[i] - l_[i];
        const auto len_substr  = d_[i] - ii;
        const auto pos_substr  = sa_[l_[i]];

        // substrと同じ出現回数のsub-substrを数える。
        int count = 0;
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
        const int num_subsubstrs = (1 + len_substr) * len_substr / 2;
        const double spurity = static_cast<double>(count) / num_subsubstrs;

        return spurity;
    }

    double loose_purity(const int i, const int ii) const {
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
        const int num_subsubstrs = (1 + len_substr) * len_substr / 2;
        const double lpurity = support / num_subsubstrs;

        return lpurity;
    }

    const std::vector<Char>& input_;
    std::vector<index_type> sa_;
    std::vector<index_type>  l_;
    std::vector<index_type>  r_;
    std::vector<index_type>  d_;
    index_type  num_nodes_;
    std::vector<index_type> suffix_to_parent_node_;
    std::vector<index_type> node_to_parent_node_;
};


#endif  /* SUBSTRINGS_HPP */
