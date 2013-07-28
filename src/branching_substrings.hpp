#ifndef BRANCHING_SUBSTRINGS_H
#define BRANCHING_SUBSTRINGS_H

#include <cmath>
#include <map>
#include <stack>
#include <stdexcept>
#include <vector>
#include <boost/iterator/iterator_facade.hpp>
#include "esa.hxx"


template <class Char, class Index>
struct core {
    typedef Index index_type;
    struct substr {
        typedef typename std::vector<Char>::const_iterator iterator;
        typedef typename std::vector<Char>::const_iterator const_iterator;

        index_type pos()       const { return parent_->sa_[parent_->l_[i_]]; }
        std::vector<index_type> allpos() const {
            return std::vector<index_type>(parent_->sa_.begin() + parent_->l_[i_], parent_->sa_.begin() + parent_->r_[i_]);
        }
        index_type length()    const { return parent_->d_[i_]; }
        index_type frequency() const { return parent_->r_[i_] - parent_->l_[i_]; }

        iterator begin() { return parent_->input_.begin() + pos(); }
        iterator end()   { return parent_->input_.begin() + pos() + length(); }
        const_iterator begin() const { return parent_->input_.begin() + pos(); }
        const_iterator end()   const { return parent_->input_.begin() + pos() + length(); }

        substr(const core* parent, int i)
            : parent_(parent), i_(i)
        {}

    private:
        const core* parent_;
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

        substring_iterator(const core* parent, int i)
            : parent_(parent), i_(i)
        {}

        substring_iterator<Value> parent() {
            return substring_iterator(parent_, parent_->node_to_parent_node_[i_]);
        }

        substring_iterator<Value> suffix() {
            return substring_iterator(parent_, parent_->suffix_link_[i_]);
        }

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

        const core* parent_;
        int i_;
    };

public:
    typedef substring_iterator<substr> iterator;
    typedef substring_iterator<const substr> const_iterator;

    core(const std::vector<Char>& input, const size_t alphabet_size)
        : input_(input),
          sa_(input.size()),
          l_(input.size()),
          r_(input.size()),
          d_(input.size()),
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

        // Dummy node
        // These values are designed so that they can be used as well as those
        // of the other nodes.
        l_[num_nodes_] = 0;
        r_[num_nodes_] = 0;
        d_[num_nodes_] = 0;

        // node_to_parent_node[i]: ノードiの親ノードの番号（post-order）。
        // suffix_to_parent_node[k]: 接尾辞input[k..$]に対応する葉ノードの、親ノードのpost-order順の番号。
        node_to_parent_node_.resize(num_nodes_);
        std::vector<index_type> suffix_to_parent_node(input.size() + 1);
        suffix_to_parent_node[input.size()] = num_nodes_;  // 接尾辞input[$..$]
        {
            std::stack<index_type> stk;  // the top of the stack is a current parent node
            stk.push(num_nodes_);  // put the dummy node, which will be the parent of the root node
            index_type next_node = num_nodes_ - 1;  // a node to consider next
            index_type i = input.size() - 1;  // a current suffix, the i-th suffix in the suffix array
            // narrow the range [l, r) to find the immediate parent of the i-th node
            while (next_node >= 0 && l_[next_node] <= i && i < r_[next_node]) {
                node_to_parent_node_[next_node] = stk.top();
                stk.push(next_node);
                --next_node;
            }
            while (i >= 0) {
                // widen the range [l, r) to find the lowest ancestor of the i-th node
                while (!(l_[stk.top()] <= i && i < r_[stk.top()])) {
                    stk.pop();
                }
                // narrow the range [l, r) to find the immediate parent of the i-th node
                while (next_node >= 0 && l_[next_node] <= i && i < r_[next_node]) {
                    node_to_parent_node_[next_node] = stk.top();
                    stk.push(next_node);
                    --next_node;
                }
                suffix_to_parent_node[i] = stk.top();
                --i;
            }
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
    }

    iterator begin() { return iterator(this, 0); }
    iterator end()   { return iterator(this, num_nodes_); }
    const_iterator begin() const { return const_iterator(this, 0); }
    const_iterator end()   const { return const_iterator(this, num_nodes_); }

    index_type size() const { return num_nodes_; }

private:
    const std::vector<Char>& input_;
    std::vector<index_type> sa_;
    std::vector<index_type>  l_;
    std::vector<index_type>  r_;
    std::vector<index_type>  d_;
    index_type  num_nodes_;
    std::vector<index_type> node_to_parent_node_;
    std::vector<index_type> suffix_link_;
};


template <class Char, class Index>
struct BranchingSubstrings {
    using core = core<Char, Index>;

    typedef Index index_type;
    struct substr {
        typedef typename std::vector<Char>::const_iterator iterator;
        typedef typename std::vector<Char>::const_iterator const_iterator;

        index_type              pos()           const { return i_->pos(); }
        std::vector<index_type> allpos()        const { return i_->allpos(); }
        index_type              length()        const { return i_->length(); }
        index_type              frequency()     const { return i_->frequency(); }
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

        substr(const BranchingSubstrings* parent, typename core::const_iterator i)
            : parent_(parent), i_(i)
        {}

    private:
        const BranchingSubstrings* parent_;
        typename core::const_iterator i_;
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
            : parent_(0), i_()
        {}

        substring_iterator(const BranchingSubstrings* parent, typename core::const_iterator i)
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

        const BranchingSubstrings* parent_;
        typename core::const_iterator i_;
    };

public:
    typedef substring_iterator<substr> iterator;
    typedef substring_iterator<const substr> const_iterator;

    BranchingSubstrings(const std::vector<Char>& input, const size_t alphabet_size)
        : input_(input),
          core_(input, alphabet_size),
          count_(core_.size(), 0), // initialize the count table with an "undefined" value.
          recip_(core_.size(), 0)  // 正数は計算結果、それ以外は未計算を表わす。
    {}

    iterator begin() { return iterator(this, core_.begin()); }
    iterator end()   { return iterator(this, core_.end()); }
    const_iterator begin() const { return const_iterator(this, core_.begin()); }
    const_iterator end()   const { return const_iterator(this, core_.end()); }

private:
    uint64_t get_count(typename core::const_iterator n) const {
        const int i = n - core_.begin();

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

    double strict_purity(typename core::const_iterator n) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto len_substr  = n->length();

        // substrと同じ出現回数のsub-substrを数える。
        const uint64_t count = get_count(n);

        // strict purity of substr
        const uint64_t num_subsubstrs = static_cast<uint64_t>(1 + len_substr) * len_substr / 2;
        const double spurity = static_cast<double>(count) / num_subsubstrs;

        return spurity;
    }

    double get_reciprocal(typename core::const_iterator n) const {
        const int i = n - core_.begin();

        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        // substrの全部分文字列の頻度の逆数の総和recip[i]を求める。
        if (recip_[i] > 0) {
            return recip_[i];
        }
        else {
            double recip = 0;
            {
                // substrの末尾を0文字以上削って得られるsub-substrについて考える。
                for (typename core::const_iterator m = n, p = n.parent(); m != core_.end(); m = p, p = p.parent()) {
                    const auto num_subsubstrs_of_same_frequency = m->length() - p->length();
                    const auto freq_subsubstr = m->frequency();
                    const double r = 1.0 / freq_subsubstr;
                    recip += num_subsubstrs_of_same_frequency * r;
                }
            }
            {
                // substrの先頭を1文字以上削ったsub-substrを考える。
                const auto m = n.suffix();
                if (m != core_.end()) {
                    recip += get_reciprocal(m);
                }
            }

            // memoize
            recip_[i] = recip;
            return recip_[i];
        }
    }

    double loose_purity(typename core::const_iterator n) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto freq_substr = n->frequency();
        const auto len_substr  = n->length();

        const double rel = freq_substr * get_reciprocal(n);

        // loose purity of substr
        const uint64_t num_subsubstrs = static_cast<uint64_t>(1 + len_substr) * len_substr / 2;
        const double lpurity = rel / num_subsubstrs;

        return lpurity;
    }

    std::map<Char, int> left_extensions(typename core::const_iterator n) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。

        std::map<Char, int> char_dist;
        for (auto pos : n->allpos()) {
            const auto& c = input_[pos - 1];
            char_dist[c] += 1;
        }

        return char_dist;
    }

    std::map<Char, int> right_extensions(typename core::const_iterator n) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto len_substr = n->length();

        std::map<Char, int> char_dist;
        for (auto pos : n->allpos()) {
            const auto& c = input_[pos + len_substr];
            char_dist[c] += 1;
        }

        return char_dist;
    }

    double left_universality(typename core::const_iterator n) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto freq_substr = n->frequency();

        std::map<Char, int> char_dist = left_extensions(n);

        double e = 0;
        for (const auto& kv : char_dist) {
            const double p = static_cast<double>(kv.second) / freq_substr;
            e += -p * std::log(p);
        }
        const double u = 1 - std::exp(-e);

        return u;
    }

    double right_universality(typename core::const_iterator n) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto freq_substr = n->frequency();

        std::map<Char, int> char_dist = right_extensions(n);

        double e = 0;
        for (const auto& kv : char_dist) {
            const double p = static_cast<double>(kv.second) / freq_substr;
            e += -p * std::log(p);
        }
        const double u = 1 - std::exp(-e);

        return u;
    }

    const std::vector<Char>& input_;
    core core_;
    mutable std::vector<uint64_t> count_;
    mutable std::vector<double>   recip_;
};


#endif  /* BRANCHING_SUBSTRINGS_H */
