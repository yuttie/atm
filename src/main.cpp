#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/bind.hpp>
#include <boost/cstdint.hpp>
#include <boost/function.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include "pstade/oven/algorithm.hpp"
#include "pstade/oven/copied.hpp"
#include "pstade/oven/file_range.hpp"
#include "pstade/oven/make_range.hpp"
#include "pstade/oven/utf8_decoded.hpp"
#include "pstade/oven/utf8_encoded.hpp"
#include "cmdline.h"
#include "esa.hxx"


using namespace std;
namespace oven = pstade::oven;


typedef boost::int32_t index_type;
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

template <class index_type>
inline
int find_parent_of_node(const int i,
                        const vector<index_type>& l,
                        const vector<index_type>& r,
                        const vector<index_type>& d)
{
    int j = i + 1;
    while (!(d[j] < d[i] && l[j] <= l[i] && r[i] <= r[j])) {
        ++j;
    }

    return j;
}

template <class Char, class Index = int32_t>
struct SubStrings {
private:
    typedef Index index_type;

public:
    struct substr {
        typedef typename vector<Char>::const_iterator iterator;

        index_type pos()       const { return pos_; }
        index_type length()    const { return length_; }
        index_type frequency() const { return frequency_; }
        double     purity()    const { return purity_; }

        iterator begin() const {
            return first_;
        }

        iterator end() const {
            return first_ + length_;
        }

        substr(iterator first,
               index_type pos,
               index_type length,
               index_type frequency,
               double purity)
            : first_(first),
              pos_(pos),
              length_(length),
              frequency_(frequency),
              purity_(purity)
        {}

    private:
        iterator first_;
        index_type pos_;
        index_type length_;
        index_type frequency_;
        double     purity_;
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
            return parent_->at(i_);
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
          suffix_to_parent_node_(input.size(), -1)
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
    }

    iterator begin() const {
        return iterator(this, 0);
    }

    iterator end() const {
        return iterator(this, num_nodes_ - 1);
    }

private:
    substr at(const int i) const {
        // ここではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
        const auto freq_substr = r_[i] - l_[i];
        const auto len_substr  = d_[i];
        const auto pos_substr  = sa_[l_[i]];

        // substrと同じ出現回数のsub-substrを数える。
        int count = 0;
        {
            // ノードiの親ノードjを見つける。
            auto j = find_parent_of_node(i, l_, r_, d_);

            // substrの末尾を0文字以上削って得られるsub-substrの内で、出現
            // 回数がsubstrと同じものの数はd[i] - d[j]である。
            count += d_[i] - d_[j];
        }
        for (int j = 1; j < len_substr; ++j) {
            // substrの先頭をj文字削ったsub-substrを考える。

            // sub-substrに対応するノードを見つける。
            auto k = suffix_to_parent_node_[pos_substr + j];  // 接尾辞input[(pos_substr + j)..$]に対応する葉ノードの親ノード
            const auto len_subsubstr = len_substr - j;
            while (d_[k] > len_subsubstr) ++k;  // d[k] == len_subsubstr ならば、ノードkはsub-substrに対応するノード。

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
                    auto m = find_parent_of_node(k, l_, r_, d_);

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

        // return
        return substr(input_.begin() + pos_substr,
                      pos_substr,
                      len_substr,
                      freq_substr,
                      purity);
    }

    const vector<Char>& input_;
    vector<index_type> sa_;
    vector<index_type>  l_;
    vector<index_type>  r_;
    vector<index_type>  d_;
    index_type  num_nodes_;
    vector<index_type> suffix_to_parent_node_;
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

        // suffix array
        vector<index_type> sa(input.size());
        vector<index_type>  l(input.size());
        vector<index_type>  r(input.size());
        vector<index_type>  d(input.size());
        index_type  num_nodes;
        int err = esaxx(input.begin(),
                        sa.begin(),
                        l.begin(), r.begin(), d.begin(),
                        static_cast<index_type>(input.size()),
                        static_cast<index_type>(alphabet_size),
                        num_nodes);
        if (err) throw runtime_error("saisxx failed to construct a suffix array.");

        // suffix_to_parent_node[k]: 接尾辞input[k..$]に対応する葉ノードの、親ノードのpost-order順の番号。
        // post-order巡回により、直接の親が最初に値を設定する（最初かどうかは-1かどうかで判定する）。
        vector<index_type> suffix_to_parent_node(input.size(), -1);
        for (int i = 0; i < num_nodes; ++i) {
            // ノードi直下の全ての葉ノードjについて、接尾辞input[k..$]からノードiへのリンクを張る
            for (int j = l[i]; j < r[i]; ++j) {
                const auto k = sa[j];
                if (suffix_to_parent_node[k] < 0) {
                    suffix_to_parent_node[k] = i;
                }
            }
        }

        // process all internal nodes
        for (int i = 0; i < num_nodes - 1; ++i) {  // ルート以外のノードをpost-order巡回する。
            // このループではノードi（iはpost-orderでの番号）に対応する部分文字列substrを扱う。
            const auto freq_substr = r[i] - l[i];
            const auto len_substr  = d[i];
            const auto pos_substr  = sa[l[i]];

            // substrと同じ出現回数のsub-substrを数える。
            int count = 0;
            {
                // ノードiの親ノードjを見つける。
                auto j = find_parent_of_node(i, l, r, d);

                // substrの末尾を0文字以上削って得られるsub-substrの内で、出現
                // 回数がsubstrと同じものの数はd[i] - d[j]である。
                count += d[i] - d[j];
            }
            for (int j = 1; j < len_substr; ++j) {
                // substrの先頭をj文字削ったsub-substrを考える。

                // sub-substrに対応するノードを見つける。
                auto k = suffix_to_parent_node[pos_substr + j];  // 接尾辞input[(pos_substr + j)..$]に対応する葉ノードの親ノード
                const auto len_subsubstr = len_substr - j;
                while (d[k] > len_subsubstr) ++k;  // d[k] == len_subsubstr ならば、ノードkはsub-substrに対応するノード。

                if (d[k] < len_subsubstr) {
                    // このsub-substrは1回しか出現していない。
                    // 今考えているsubstrは内部ノードに対応しており、出現頻度が
                    // 2以上なので、purityを考える場合は、このsub-substrを無視
                    // してよい。
                }
                else {
                    // このsub-substrは2回以上出現している。
                    const auto freq_subsubstr = r[k] - l[k];
                    if (freq_subsubstr == freq_substr) {
                        // ノードkの親ノードmを見つける。
                        auto m = find_parent_of_node(k, l, r, d);

                        // sub-substrの末尾を0文字以上削って得られる
                        // sub-sub-substrの内で、出現回数がsub-substrと同じもの
                        // の数はd[k] - d[m]である。
                        count += d[k] - d[m];
                    }
                }
            }

            // purity of substr
            const int num_subsubstrs = (1 + len_substr) * len_substr / 2;
            const double purity = static_cast<double>(count) / num_subsubstrs;

            // 出力
            std::cout << pos_substr << "\t" << len_substr << "\t" << freq_substr << "\t" << purity;
            if (p.exist("show-substring")) {
                auto first = input.begin() + pos_substr;
                auto substr = oven::make_range(first, first + len_substr);
                std::cout << "\t";
                oven::copy(substr, ostream_iterator<byte_type>(std::cout, ""));
            }
            std::cout << "\n";
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

            // suffix array
            vector<index_type> sa(input.size());
            int err = saisxx(input.begin(),
                             sa.begin(),
                             static_cast<index_type>(input.size()),
                             static_cast<index_type>(alphabet_size));
            if (err) throw runtime_error("saisxx failed to construct a suffix array.");

            // output suffixes
            for (auto i : sa) {
                auto suffix = oven::make_range(input.begin() + i, input.end());
                oven::copy(suffix | oven::transformed(lookup_by(id2char)) | oven::utf8_encoded,
                           ostream_iterator<byte_type>(std::cout, ""));
            }
        }
        else if (alphabet_size <= 0x10000) {
            typedef boost::uint16_t id_type;
            throw runtime_error("Not implemented yet.");
        }
        else {
            typedef boost::uint32_t id_type;
            throw runtime_error("Not implemented yet.");
        }
    }
}
