#ifndef BLUMER_SUBSTRINGS_H
#define BLUMER_SUBSTRINGS_H

#include <cmath>
#include <map>
#include <stack>
#include <stdexcept>
#include <vector>
#include <boost/iterator/iterator_facade.hpp>
#include "esa.hxx"


template <class Char, class Index>
struct BlumerSubstrings {
    typedef Index index_type;
    struct substr {
        typedef typename std::vector<Char>::const_iterator iterator;
        typedef typename std::vector<Char>::const_iterator const_iterator;

        index_type              pos()           const { return parent_->sa_[parent_->l_[i_]]; }
        std::vector<index_type> allpos()        const { return parent_->allpos(i_); }
        index_type              length()        const { return parent_->d_[i_]; }
        index_type              frequency()     const { return parent_->r_[i_] - parent_->l_[i_]; }
        double                  spurity()       const { return parent_->strict_purity(i_); }
        double                  lpurity()       const { return parent_->loose_purity(i_); }
        double                  luniversality() const { return parent_->left_universality(i_); }
        double                  runiversality() const { return parent_->right_universality(i_); }

        iterator begin() {
            return parent_->input_.begin() + pos();
        }

        iterator end() {
            return parent_->input_.begin() + pos() + length();
        }

        const_iterator begin() const {
            return parent_->input_.begin() + pos();
        }

        const_iterator end() const {
            return parent_->input_.begin() + pos() + length();
        }

        substr(const BlumerSubstrings* parent, int i)
            : parent_(parent), i_(i)
        {}

    private:
        const BlumerSubstrings* parent_;
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

        substring_iterator(const BlumerSubstrings* parent, int i)
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
            return substr(parent_, parent_->selected_node_indices_[i_]);
        }

        const BlumerSubstrings* parent_;
        int i_;
    };

public:
    typedef substring_iterator<substr> iterator;
    typedef substring_iterator<const substr> const_iterator;

    BlumerSubstrings(const std::vector<Char>& input, const size_t alphabet_size)
        : input_(input),
          sa_(input.size()),
          l_(input.size()),
          r_(input.size()),
          d_(input.size()),
          count_(),
          recip_(),
          node_to_parent_node_(),
          selected_node_indices_()
    {
        // suffix array
        int err = esaxx(input_.begin(),
                        sa_.begin(),
                        l_.begin(), r_.begin(), d_.begin(),
                        static_cast<index_type>(input_.size()),
                        static_cast<index_type>(alphabet_size),
                        num_nodes_);
        if (err) throw std::runtime_error("saisxx failed to construct a suffix array.");

        // Dummy node
        // These values are designed so that they can be used as well as those
        // of the other nodes.
        l_[num_nodes_] = 0;
        r_[num_nodes_] = 0;
        d_[num_nodes_] = 0;

        // initialize the count table with an "undefined" value.
        // 正数は計算結果、それ以外は未計算を表わす。
        count_.assign(num_nodes_, 0);
        recip_.assign(num_nodes_, 0);

        // suffix_to_parent_node[k]: 接尾辞input[k..$]に対応する葉ノードの、親ノードのpost-order順の番号。
        // 逆向きpost-order巡回により、直接の親が最後に値を設定（上書き）する。
        std::vector<index_type> suffix_to_parent_node(input.size() + 1);
        suffix_to_parent_node[input.size()] = num_nodes_;  // 接尾辞input[$..$]
        for (int i = num_nodes_ - 1; i >= 0; --i) {
            // ノードi直下の全ての葉ノードjについて、接尾辞input[k..$]からノードiへのリンクを張る
            for (int j = l_[i]; j < r_[i]; ++j) {
                const auto k = sa_[j];
                suffix_to_parent_node[k] = i;
            }
        }

        // node_to_parent_node[i]: ノードiの親ノードの番号（post-order）。
        node_to_parent_node_.resize(num_nodes_);
        node_to_parent_node_[num_nodes_ - 1] = num_nodes_;  // transfers to the dummy node.
        std::stack<index_type> stk;
        stk.push(num_nodes_ - 1);
        for (int i = num_nodes_ - 2; i >= 0; --i) {
            while (!(l_[stk.top()] <= l_[i] && r_[i] <= r_[stk.top()])) {
                stk.pop();
            }
            node_to_parent_node_[i] = stk.top();
            stk.push(i);
        }

        // suffix_link_
        suffix_link_.resize(num_nodes_);
        suffix_link_[num_nodes_ - 1] = num_nodes_;  // transfers to the dummy node.
        for (int i = 0; i < num_nodes_ - 1; ++i) {
            // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
            const auto len_substr  = d_[i];
            const auto pos_substr  = sa_[l_[i]];
            // substrの先頭を1文字削ったsub-substrを考える。
            const auto len_subsubstr = len_substr - 1;
            // sub-substrに対応するノードを見つける。
            // {内部,葉}ノードに対応する部分文字列から先頭の1文字を削って得られ
            // る文字列には、必ず対応する{内部,葉}ノードが存在する。
            auto j = suffix_to_parent_node[pos_substr + 1];  // 接尾辞input[(pos_substr + j)..$]に対応する葉ノードの親ノード
            while (d_[j] > len_subsubstr) j = node_to_parent_node_[j];  // d[j] == len_subsubstr ならば、ノードjはsub-substrに対応するノード。
            suffix_link_[i] = j;
        }

        // split the nodes into Blumer's equivalence classes
        int num_classes = 0;
        std::vector<int> class_ids(num_nodes_, -1);
        std::vector<std::vector<int>> classes;
        for (int i = 0; i < num_nodes_; ++i) {
            const int cid = get_class_id(i, num_classes, class_ids);
            classes.resize(num_classes);
            classes[cid].push_back(i);
        }

        // find the node corresponding to the longest substring in each class
        for (int i = 0; i < num_classes; ++i) {
            int max_length = -1;
            int index_max;
            for (int j = 0; j < classes[i].size(); ++j) {
                const int length = d_[classes[i][j]];
                if (length > max_length) {
                    max_length = length;
                    index_max = j;
                }
            }
            selected_node_indices_.push_back(classes[i][index_max]);
        }
    }

    iterator begin() {
        return iterator(this, 0);
    }

    iterator end() {
        return iterator(this, selected_node_indices_.size());
    }

    const_iterator begin() const {
        return const_iterator(this, 0);
    }

    const_iterator end() const {
        return const_iterator(this, selected_node_indices_.size());
    }

private:
    std::vector<index_type> allpos(const int i) const {
        return std::vector<index_type>(sa_.begin() + l_[i], sa_.begin() + r_[i]);
    }

    int get_class_id(const int i, int& num_classes, std::vector<int>& class_ids) const {
        if (class_ids[i] >= 0) {
            return class_ids[i];
        }
        else {
            const auto freq_substr = r_[i] - l_[i];

            const int j = suffix_link_[i];
            const auto freq_subsubstr = r_[j] - l_[j];

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

    uint64_t get_count(const int i) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        // substrと同じ出現回数のsub-substrの数count[i]を数える。
        if (count_[i] > 0) {
            return count_[i];
        }
        else {
            const auto freq_substr = r_[i] - l_[i];

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
                count += d_[i] - d_[j];
            }
            {
                // substrの先頭を1文字以上削ったsub-substrを考える。
                const auto j = suffix_link_[i];
                const auto freq_subsubstr = r_[j] - l_[j];
                if (freq_subsubstr == freq_substr) {
                    count += get_count(j);
                }
            }

            // memoize
            count_[i] = count;
            return count_[i];
        }
    }

    double strict_purity(const int i) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto len_substr  = d_[i];

        // substrと同じ出現回数のsub-substrを数える。
        const uint64_t count = get_count(i);

        // strict purity of substr
        const uint64_t num_subsubstrs = static_cast<uint64_t>(1 + len_substr) * len_substr / 2;
        const double spurity = static_cast<double>(count) / num_subsubstrs;

        return spurity;
    }

    double get_reciprocal(const int i) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        // substrの全部分文字列の頻度の逆数の総和recip[i]を求める。
        if (recip_[i] > 0) {
            return recip_[i];
        }
        else {
            double recip = 0;
            {
                // substrの末尾を0文字以上削って得られるsub-substrについて考える。
                for (index_type j = i, k = node_to_parent_node_[i]; j < num_nodes_; j = k, k = node_to_parent_node_[k]) {
                    const auto num_subsubstrs_of_same_frequency = d_[j] - d_[k];
                    const auto freq_subsubstr = r_[j] - l_[j];
                    const double r = 1.0 / freq_subsubstr;
                    recip += num_subsubstrs_of_same_frequency * r;
                }
            }
            {
                // substrの先頭を1文字以上削ったsub-substrを考える。
                const auto j = suffix_link_[i];
                if (j < num_nodes_) {
                    recip += get_reciprocal(j);
                }
            }

            // memoize
            recip_[i] = recip;
            return recip_[i];
        }
    }

    double loose_purity(const int i) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto freq_substr = r_[i] - l_[i];
        const auto len_substr  = d_[i];

        const double rel = freq_substr * get_reciprocal(i);

        // loose purity of substr
        const uint64_t num_subsubstrs = static_cast<uint64_t>(1 + len_substr) * len_substr / 2;
        const double lpurity = rel / num_subsubstrs;

        return lpurity;
    }

    std::map<Char, int> left_extensions(const int i) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。

        std::map<Char, int> char_dist;
        for (auto j = l_[i]; j < r_[i]; ++j) {
            const auto pos = sa_[j];
            const auto& c = input_[pos - 1];
            char_dist[c] += 1;
        }

        return char_dist;
    }

    std::map<Char, int> right_extensions(const int i) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto len_substr = d_[i];

        std::map<Char, int> char_dist;
        for (auto j = l_[i]; j < r_[i]; ++j) {
            const auto pos = sa_[j];
            const auto& c = input_[pos + len_substr];
            char_dist[c] += 1;
        }

        return char_dist;
    }

    double left_universality(const int i) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto freq_substr = r_[i] - l_[i];

        std::map<Char, int> char_dist = left_extensions(i);

        double e = 0;
        for (const auto& kv : char_dist) {
            const double p = static_cast<double>(kv.second) / freq_substr;
            e += -p * std::log(p);
        }
        const double u = 1 - std::exp(-e);

        return u;
    }

    double right_universality(const int i) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto freq_substr = r_[i] - l_[i];

        std::map<Char, int> char_dist = right_extensions(i);

        double e = 0;
        for (const auto& kv : char_dist) {
            const double p = static_cast<double>(kv.second) / freq_substr;
            e += -p * std::log(p);
        }
        const double u = 1 - std::exp(-e);

        return u;
    }

    const std::vector<Char>& input_;
    std::vector<index_type> sa_;
    std::vector<index_type>  l_;
    std::vector<index_type>  r_;
    std::vector<index_type>  d_;
    mutable std::vector<uint64_t> count_;
    mutable std::vector<double>   recip_;
    index_type  num_nodes_;
    std::vector<index_type> node_to_parent_node_;
    std::vector<index_type> suffix_link_;
    std::vector<index_type> selected_node_indices_;
};


#endif  /* BLUMER_SUBSTRINGS_H */