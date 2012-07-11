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
#include <boost/lambda/lambda.hpp>
#include <boost/optional/optional.hpp>
#include "pstade/oven/algorithm.hpp"
#include "pstade/oven/copied.hpp"
#include "pstade/oven/file_range.hpp"
#include "pstade/oven/filtered.hpp"
#include "pstade/oven/filterer.hpp"
#include "pstade/oven/make_range.hpp"
#include "pstade/oven/transformed.hpp"
#include "pstade/oven/utf8_decoded.hpp"
#include "pstade/oven/utf8_encoded.hpp"
#include "cmdline.h"
#include "esa.hxx"


using namespace std;
namespace lambda = boost::lambda;
namespace oven = pstade::oven;


typedef boost::uint8_t byte_type;
typedef boost::int32_t index_type;

template <class K, class V>
typename boost::function<V (K)> tr_by(const map<K, V>& m) {
    struct translate {
        typedef typename map<K, V>::mapped_type result_type;

        result_type operator()(const map<K, V>& m, const typename map<K, V>::key_type& k) const {
            return m.at(k);
        }
    };

    return boost::bind(translate(), ref(m), _1);
}

template <class V>
typename boost::function<V (typename vector<V>::size_type)> tr_by(const vector<V>& v) {
    struct translate {
        typedef typename vector<V>::value_type result_type;

        result_type operator()(const vector<V>& v, const typename vector<V>::size_type& k) const {
            return v[k];
        }
    };

    return boost::bind(translate(), ref(v), _1);
}

template <class Char>
struct SubStrings {
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

        substring_iterator(const SubStrings* parent, int i)
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

        const SubStrings* parent_;
        int i_;
    };

public:
    typedef substring_iterator<substr> iterator;
    typedef substring_iterator<const substr> const_iterator;

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
    ResultPrinter(std::ostream& os, bool show_substr, bool exclude_newline)
        : os_(os),
          to_unicode_char_(),
          show_substr_(show_substr),
          exclude_newline_(exclude_newline)
    {}

    template <class F>
    ResultPrinter(std::ostream& os, F to_unicode_char, bool show_substr, bool exclude_newline)
        : os_(os),
          to_unicode_char_(to_unicode_char),
          show_substr_(show_substr),
          exclude_newline_(exclude_newline)
    {}

    void print_header() {
        os_ << "position" << "\t" << "length" << "\t" << "frequency" << "\t" << "purity" << "\n";
    }

    template <class S>
    void print(const S& substr) {
        os_ << substr.pos() << "\t" << substr.length() << "\t" << substr.frequency() << "\t" << substr.purity();
        if (show_substr_) {
            os_ << "\t";
            if (to_unicode_char_) {
                auto encoded = substr | oven::transformed(*to_unicode_char_) | oven::utf8_encoded;
                if (exclude_newline_) {
                    oven::copy(encoded, oven::filterer(lambda::_1 != '\n') |= std::ostream_iterator<byte_type>(os_));
                }
                else {
                    oven::copy(encoded, std::ostream_iterator<byte_type>(os_));
                }
            }
            else {
                if (exclude_newline_) {
                    oven::copy(substr, oven::filterer(lambda::_1 != '\n') |= std::ostream_iterator<byte_type>(os_));
                }
                else {
                    oven::copy(substr, std::ostream_iterator<byte_type>(os_));
                }
            }
        }
        os_ << "\n";
    }

private:
    typedef boost::uint32_t unicode_char_type;
    typedef boost::uint32_t largest_id_type;
    std::ostream& os_;
    boost::optional<boost::function<unicode_char_type (largest_id_type)>> to_unicode_char_;
    bool show_substr_;
    bool exclude_newline_;
};

struct substring_constraint {
    boost::optional<index_type> min_length;
    boost::optional<index_type> max_length;
    boost::optional<index_type> min_frequency;
    boost::optional<index_type> max_frequency;
    boost::optional<double>     min_purity;
    boost::optional<double>     max_purity;
};

template <class Substring>
struct satisfy {
    satisfy(const substring_constraint& c)
        : c_(c)
    {}

    bool operator()(const Substring& substr) const {
        return (!c_.min_length    || substr.length()    >= *c_.min_length)
            && (!c_.max_length    || substr.length()    <= *c_.max_length)
            && (!c_.min_frequency || substr.frequency() >= *c_.min_frequency)
            && (!c_.max_frequency || substr.frequency() <= *c_.max_frequency)
            && (!c_.min_purity    || substr.purity()    >= *c_.min_purity)
            && (!c_.max_purity    || substr.purity()    <= *c_.max_purity);
    }

private:
    const substring_constraint c_;
};

int main(int argc, char* argv[]) {
    // command line
    cmdline::parser p;
    p.add("help", 'h', "");
    p.add<string>("number-format", 'F', "", false, "fixed", cmdline::oneof<string>("fixed", "scientific"));
    p.add<string>("mode", 'm', "", false, "binary", cmdline::oneof<string>("binary", "text"));
    p.add("show-substring", 's', "");
    p.add("exclude-newline", 'N', "");
    p.add<index_type>("longer",  0, "", false);
    p.add<index_type>("shorter", 0, "", false);
    p.add<index_type>("more-frequent",   0, "", false);
    p.add<index_type>("more-infrequent", 0, "", false);
    p.add<double>("purer",   0, "", false);
    p.add<double>("impurer", 0, "", false);
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

    // number format
    cout.setf(p.get<string>("number-format") == "fixed" ? ios::fixed : ios::scientific,
              ios::floatfield);

    // substring filter
    substring_constraint constraint = {
        p.exist("longer")          ? boost::make_optional(p.get<index_type>("longer"))          : boost::none,
        p.exist("shorter")         ? boost::make_optional(p.get<index_type>("shorter"))         : boost::none,
        p.exist("more-frequent")   ? boost::make_optional(p.get<index_type>("more-frequent"))   : boost::none,
        p.exist("more-infrequent") ? boost::make_optional(p.get<index_type>("more-infrequent")) : boost::none,
        p.exist("purer")           ? boost::make_optional(p.get<double>("purer"))               : boost::none,
        p.exist("impurer")         ? boost::make_optional(p.get<double>("impurer"))             : boost::none
    };

    if (p.get<string>("mode") == "binary") {
        typedef boost::uint8_t char_type;
        typedef boost::uint8_t id_type;
        typedef typename SubStrings<id_type>::substr substr_type;

        // alphabets
        const size_t alphabet_size = 0x100;

        // input
        const vector<id_type> input = is | oven::copied;

        // printer
        ResultPrinter printer(std::cout, p.exist("show-substring"), p.exist("exclude-newline"));

        // enumerate substrings
        SubStrings<id_type> substrs(input, alphabet_size);

        printer.print_header();
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
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

        // printer
        ResultPrinter printer(std::cout, tr_by(id2char), p.exist("show-substring"), p.exist("exclude-newline"));

        if (alphabet_size <= 0x100) {
            typedef boost::uint8_t id_type;
            typedef typename SubStrings<id_type>::substr substr_type;

            // map: char -> id
            map<char_type, id_type> char2id;
            for (size_t id = 0; id < alphabet_size; ++id) {
                char2id[id2char[id]] = id;
            }

            // input
            const vector<id_type> input = is | oven::utf8_decoded | oven::transformed(tr_by(char2id)) | oven::copied;

            // enumerate substrings
            SubStrings<id_type> substrs(input, alphabet_size);

            printer.print_header();
            for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
                printer.print(substr);
            }
        }
        else if (alphabet_size <= 0x10000) {
            typedef boost::uint16_t id_type;
            typedef typename SubStrings<id_type>::substr substr_type;

            // map: char -> id
            map<char_type, id_type> char2id;
            for (size_t id = 0; id < alphabet_size; ++id) {
                char2id[id2char[id]] = id;
            }

            // input
            const vector<id_type> input = is | oven::utf8_decoded | oven::transformed(tr_by(char2id)) | oven::copied;

            // enumerate substrings
            SubStrings<id_type> substrs(input, alphabet_size);

            printer.print_header();
            for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
                printer.print(substr);
            }
        }
        else {
            typedef boost::uint32_t id_type;
            typedef typename SubStrings<id_type>::substr substr_type;

            // map: char -> id
            map<char_type, id_type> char2id;
            for (size_t id = 0; id < alphabet_size; ++id) {
                char2id[id2char[id]] = id;
            }

            // input
            const vector<id_type> input = is | oven::utf8_decoded | oven::transformed(tr_by(char2id)) | oven::copied;

            // enumerate substrings
            SubStrings<id_type> substrs(input, alphabet_size);

            printer.print_header();
            for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
                printer.print(substr);
            }
        }
    }
}
