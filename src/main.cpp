#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <stack>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/bind.hpp>
#include <boost/cstdint.hpp>
#include <boost/function.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include "pstade/oven/algorithm.hpp"
#include "pstade/oven/any_output_iterator.hpp"
#include "pstade/oven/copied.hpp"
#include "pstade/oven/file_range.hpp"
#include "pstade/oven/make_range.hpp"
#include "pstade/oven/transformed.hpp"
#include "pstade/oven/transformer.hpp"
#include "pstade/oven/utf8_decoded.hpp"
#include "pstade/oven/utf8_encoded.hpp"
#include "pstade/oven/utf8_encoder.hpp"
#include "cmdline.h"
#include "esa.hxx"


using namespace std;
namespace oven = pstade::oven;


typedef boost::uint8_t byte_type;

template <class K, class V>
typename boost::function<V (K)> lookup_by(const map<K, V>& m) {
    struct lookup {
        typedef typename map<K, V>::mapped_type result_type;

        result_type operator()(const map<K, V>& m, const typename map<K, V>::key_type& k) const {
            return m.find(k)->second;
        }
    };

    return boost::bind(lookup(), ref(m), _1);
}

template <class V>
typename boost::function<V (typename vector<V>::size_type)> lookup_by(const vector<V>& v) {
    struct lookup {
        typedef typename vector<V>::value_type result_type;

        result_type operator()(const vector<V>& v, const typename vector<V>::size_type& k) const {
            return v[k];
        }
    };

    return boost::bind(lookup(), ref(v), _1);
}

template <class Char, class Index = int32_t>
struct SubStrings {
private:
    typedef Index index_type;

public:
    struct substr {
        typedef typename vector<Char>::const_iterator iterator;
        typedef typename vector<Char>::const_iterator const_iterator;

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

        substr(const SubStrings* parent, int i)
            : parent_(parent), i_(i)
        {}

    private:
        const SubStrings* parent_;
        int i_;
    };

    struct iterator
        : public boost::iterator_facade<
            iterator,
            substr,
            boost::random_access_traversal_tag,
            substr,
            int>
    {
        iterator()
            : parent_(0), i_(-1)
        {}

        iterator(const SubStrings* parent, int i)
            : parent_(parent), i_(i)
        {}

    private:
        friend class boost::iterator_core_access;

        void increment() { ++i_; }

        void decrement() { --i_; }

        void advance(int n) { i_ += n; }

        int distance_to(const iterator& other) const { return other.i_ - this->i_; }

        bool equal(const iterator& other) const {
            return this->parent_ == other.parent_ && this->i_ == other.i_;
        }

        substr dereference() const {
            return substr(parent_, i_);
        }

        const SubStrings* parent_;
        int i_;
    };

    SubStrings(const vector<Char>& input, const size_t alphabet_size)
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
        if (err) throw runtime_error("saisxx failed to construct a suffix array.");

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
        stack<index_type> stk;
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
            // ノードiの親ノードjを見つける。
            auto j = node_to_parent_node_[i];

            // substrの末尾を0文字以上削って得られるsub-substrの内で、出現
            // 回数がsubstrと同じものの数はd[i] - d[j]である。
            count += d_[i] - d_[j];
        }
        for (int j = 1; j < len_substr; ++j) {
            // substrの先頭をj文字削ったsub-substrを考える。

            // sub-substrに対応するノードを見つける。
            auto k = suffix_to_parent_node_[pos_substr + j];  // 接尾辞input[(pos_substr + j)..$]に対応する葉ノードの親ノード
            const auto len_subsubstr = len_substr - j;
            while (d_[k] > len_subsubstr) k = node_to_parent_node_[k];  // d[k] == len_subsubstr ならば、ノードkはsub-substrに対応するノード。

            if (d_[k] < len_subsubstr) {
                // このsub-substrは1回しか出現していない。
                // 今考えているsubstrは内部ノードに対応しており、出現頻度が
                // 2以上なので、purityを考える場合は、このsub-substrを無視
                // してよい。
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

    const vector<Char>& input_;
    vector<index_type> sa_;
    vector<index_type>  l_;
    vector<index_type>  r_;
    vector<index_type>  d_;
    index_type  num_nodes_;
    vector<index_type> suffix_to_parent_node_;
    vector<index_type> node_to_parent_node_;
};

struct ResultPrinter {
    ResultPrinter(std::ostream& os, bool show_substr)
        : os_(os),
          outit_(std::ostream_iterator<byte_type>(os)),
          show_substr_(show_substr)
    {}

    template <class F>
    ResultPrinter(std::ostream& os, F to_unicode_char, bool show_substr)
        : os_(os),
          outit_(oven::transformer(to_unicode_char) |= oven::utf8_encoder |= std::ostream_iterator<byte_type>(os)),
          show_substr_(show_substr)
    {}

    void print_header() {
        os_ << "position" << "\t" << "length" << "\t" << "frequency" << "\t" << "purity" << "\n";
    }

    template <class S>
    void print(const S& substr) {
        os_ << substr.pos() << "\t" << substr.length() << "\t" << substr.frequency() << "\t" << substr.purity();
        if (show_substr_) {
            os_ << "\t";
            oven::copy(substr, outit_);
        }
        os_ << "\n";
    }

private:
    std::ostream& os_;
    oven::any_output_iterator<byte_type> outit_;
    bool show_substr_;
};

int main(int argc, char* argv[]) {
    // command line
    cmdline::parser p;
    p.add("help", 'h', "");
    p.add<string>("mode", 'm', "", false, "binary", cmdline::oneof<string>("binary", "text"));
    p.add("show-substring", 's', "");
    if (!p.parse(argc, argv) || p.exist("help")) {
        cout << p.error_full() << p.usage();
        return 0;
    }

    vector<string> rest_args = p.rest();
    if (rest_args.size() == 0) {
        throw runtime_error("no input filename is given.");
    }

    // an input file
    string fp = rest_args[0];
    oven::file_range<byte_type> is(fp);

    if (p.get<string>("mode") == "binary") {
        typedef boost::uint8_t char_type;
        typedef boost::uint8_t id_type;

        // alphabets
        const size_t alphabet_size = 0x100;

        // input
        const vector<id_type> input = is | oven::copied;

        // enumerate substrings
        SubStrings<id_type> substrs(input, alphabet_size);

        ResultPrinter printer(std::cout, p.exist("show-substring"));
        printer.print_header();
        for (auto substr : substrs) {
            printer.print(substr);
        }
    }
    else {
        typedef boost::uint32_t char_type;

        // alphabets
        const set<char_type> alphabets = is | oven::utf8_decoded | oven::copied;
        const size_t alphabet_size = alphabets.size();

        // map: id -> char
        const vector<char_type> id2char = alphabets | oven::copied;

        if (alphabet_size <= 0x100) {
            typedef boost::uint8_t id_type;

            // map: char -> id
            map<char_type, id_type> char2id;
            for (size_t id = 0; id < alphabet_size; ++id) {
                char2id[id2char[id]] = id;
            }

            // input
            const vector<id_type> input = is | oven::utf8_decoded | oven::transformed(lookup_by(char2id)) | oven::copied;

            // enumerate substrings
            SubStrings<id_type> substrs(input, alphabet_size);

            ResultPrinter printer(std::cout, lookup_by(id2char), p.exist("show-substring"));
            printer.print_header();
            for (auto substr : substrs) {
                printer.print(substr);
            }
        }
        else if (alphabet_size <= 0x10000) {
            typedef boost::uint16_t id_type;

            // map: char -> id
            map<char_type, id_type> char2id;
            for (size_t id = 0; id < alphabet_size; ++id) {
                char2id[id2char[id]] = id;
            }

            // input
            const vector<id_type> input = is | oven::utf8_decoded | oven::transformed(lookup_by(char2id)) | oven::copied;

            // enumerate substrings
            SubStrings<id_type> substrs(input, alphabet_size);

            ResultPrinter printer(std::cout, lookup_by(id2char), p.exist("show-substring"));
            printer.print_header();
            for (auto substr : substrs) {
                printer.print(substr);
            }
        }
        else {
            typedef boost::uint32_t id_type;

            // map: char -> id
            map<char_type, id_type> char2id;
            for (size_t id = 0; id < alphabet_size; ++id) {
                char2id[id2char[id]] = id;
            }

            // input
            const vector<id_type> input = is | oven::utf8_decoded | oven::transformed(lookup_by(char2id)) | oven::copied;

            // enumerate substrings
            SubStrings<id_type> substrs(input, alphabet_size);

            ResultPrinter printer(std::cout, lookup_by(id2char), p.exist("show-substring"));
            printer.print_header();
            for (auto substr : substrs) {
                printer.print(substr);
            }
        }
    }
}
