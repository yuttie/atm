#ifndef SAST_HPP
#define SAST_HPP

#include <stack>
#include <stdexcept>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/range/begin.hpp>
#include <boost/range/iterator.hpp>
#include <boost/range/size.hpp>
#include "esa.hxx"


template <class RandomAccessRange, class Index>
struct sast {
    typedef Index index_type;
    struct substr {
        using iterator       = typename boost::range_iterator<RandomAccessRange>::type;
        using const_iterator = typename boost::range_const_iterator<RandomAccessRange>::type;

        index_type pos()       const { return parent_->sa_[parent_->l_[i_]]; }
        std::vector<index_type> allpos() const {
            return std::vector<index_type>(parent_->sa_.begin() + parent_->l_[i_], parent_->sa_.begin() + parent_->r_[i_]);
        }
        index_type length()    const { return parent_->d_[i_]; }
        index_type frequency() const { return parent_->r_[i_] - parent_->l_[i_]; }

        iterator begin() { return boost::begin(parent_->input_) + pos(); }
        iterator end()   { return boost::begin(parent_->input_) + pos() + length(); }
        const_iterator begin() const { return boost::begin(parent_->input_) + pos(); }
        const_iterator end()   const { return boost::begin(parent_->input_) + pos() + length(); }

        substr(const sast* parent, int i)
            : parent_(parent), i_(i)
        {}

    private:
        const sast* parent_;
        int i_;
    };

private:
    template <class Value>
    struct node_iterator
        : public boost::iterator_facade<
            node_iterator<Value>,
            Value,
            boost::random_access_traversal_tag,
            Value,
            int>
    {
        node_iterator()
            : parent_(0), i_(-1)
        {}

        node_iterator(const sast* parent, int i)
            : parent_(parent), i_(i)
        {}

        template <class OtherValue>
        node_iterator(node_iterator<OtherValue> const& other)
            : parent_(other.parent_), i_(other.i_)
        {}

        node_iterator<Value> parent() {
            return node_iterator(parent_, parent_->node_to_parent_node_[i_]);
        }

        node_iterator<Value> suffix() {
            return node_iterator(parent_, parent_->suffix_link_[i_]);
        }

    private:
        friend class boost::iterator_core_access;
        template <class> friend struct node_iterator;

        void increment() { ++i_; }
        void decrement() { --i_; }
        void advance(int n) { i_ += n; }
        int distance_to(const node_iterator<Value>& other) const { return other.i_ - this->i_; }

        template <class OtherValue>
        bool equal(const node_iterator<OtherValue>& other) const {
            return this->parent_ == other.parent_ && this->i_ == other.i_;
        }

        Value dereference() const {
            return substr(parent_, i_);
        }

        const sast* parent_;
        int i_;
    };

public:
    typedef node_iterator<substr> iterator;
    typedef node_iterator<const substr> const_iterator;

    sast(const RandomAccessRange& input, const size_t alphabet_size)
        : input_(input),
          sa_(input.size()),
          l_(input.size()),
          r_(input.size()),
          d_(input.size()),
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
        suffix_to_parent_node[input.size()] = num_nodes_ - 1;  // 接尾辞input[$..$]
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
                suffix_to_parent_node[sa_[i]] = stk.top();
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
    const RandomAccessRange& input_;
    std::vector<index_type> sa_;
    std::vector<index_type>  l_;
    std::vector<index_type>  r_;
    std::vector<index_type>  d_;
    index_type  num_nodes_;
    std::vector<index_type> node_to_parent_node_;
    std::vector<index_type> suffix_link_;
};


#endif  /* SAST_HPP */
