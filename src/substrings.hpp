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
        index_type length()    const { return parent_->d_[i_]; }
        index_type frequency() const { return parent_->r_[i_] - parent_->l_[i_]; }
        double     purity()    const { return parent_->purity(i_); }

        iterator begin() const {
            return parent_->input_.begin() + pos();
        }

        iterator end() const {
            return parent_->input_.begin() + pos() + length();
        }

        substr(const Substrings* parent, int i)
            : parent_(parent), i_(i)
        {}

    private:
        const Substrings* parent_;
        int i_;
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
            : parent_(0), i_(-1)
        {}

        substring_iterator(const Substrings* parent, int i)
            : parent_(parent), i_(i)
        {}

    private:
        friend class boost::iterator_core_access;

        void increment() { ++i_; }

        void decrement() { --i_; }

        void advance(int n) { i_ += n; }

        int distance_to(const substring_iterator<Value>& other) const { return other.i_ - this->i_; }

        bool equal(const substring_iterator<Value>& other) const {
            return this->parent_ == other.parent_ && this->i_ == other.i_;
        }

        Value dereference() const {
            return substr(parent_, i_);
        }

        const Substrings* parent_;
        int i_;
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
        return iterator(this, 0);
    }

    iterator end() const {
        return iterator(this, num_nodes_ - 1);
    }

private:
    double purity(const int i) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto freq_substr = r_[i] - l_[i];
        const auto len_substr  = d_[i];
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
            auto j = node_to_parent_node_[i];

            // substrの末尾を0文字以上削って得られるsub-substrの内で、出現
            // 回数がsubstrと同じものの数はd[i] - d[j]である。
            count += d_[i] - d_[j];
        }
        for (int j = 1; j < len_substr; ++j) {
            // substrの先頭をj文字削ったsub-substrを考える。
            const auto len_subsubstr = len_substr - j;

            // sub-substrに対応するノードを見つける。
            auto k = suffix_to_parent_node_[pos_substr + j];  // 接尾辞input[(pos_substr + j)..$]に対応する葉ノードの親ノード
            while (d_[k] > len_subsubstr) k = node_to_parent_node_[k];  // d[k] == len_subsubstr ならば、ノードkはsub-substrに対応するノード。

            if (d_[k] < len_subsubstr) {
                // このsub-substrは1回しか出現していない
                // （対応するノードがノードkの子の葉ノードである）。
                // 今考えているsubstrは出現頻度が2以上の、内部ノードに対応する
                // substrなので、そのsub-substrの出現頻度が2未満であるはずがな
                // い。
                // よって、このような場合はあり得ない。
            }
            else {
                // このsub-substrは2回以上出現している。
                const auto freq_subsubstr = r_[k] - l_[k];
                if (freq_subsubstr == freq_substr) {
                    // ノードkの親ノードmを見つける。
                    auto m = node_to_parent_node_[k];

                    // sub-substrの末尾を0文字以上削って得られる
                    // sub-sub-substrの内で、出現回数がsub-substrと同じもの
                    // の数はd[k] - d[m]である。
                    count += d_[k] - d_[m];
                }
            }
        }

        // purity of substr
        const int num_subsubstrs = (1 + len_substr) * len_substr / 2;
        const double purity = static_cast<double>(count) / num_subsubstrs;

        return purity;
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
