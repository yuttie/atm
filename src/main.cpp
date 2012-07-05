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
#include "pstade/oven/algorithm.hpp"
#include "pstade/oven/copied.hpp"
#include "pstade/oven/file_range.hpp"
#include "pstade/oven/make_range.hpp"
#include "pstade/oven/utf8_decoded.hpp"
#include "pstade/oven/utf8_encoded.hpp"
#include "cmdline.h"
#include "sais.hxx"


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

int main(int argc, char* argv[]) {
    // command line
    cmdline::parser p;
    p.add("help", 'h', "");
    p.add<string>("mode", 'm', "", false, "binary", cmdline::oneof<string>("binary", "text"));
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
        int err = saisxx(input.begin(),
                         sa.begin(),
                         static_cast<index_type>(input.size()),
                         static_cast<index_type>(alphabet_size));
        if (err) throw runtime_error("saisxx failed to construct a suffix array.");

        // output suffixes
        for (auto i : sa) {
            auto suffix = oven::make_range(input.begin() + i, input.end());
            oven::copy(suffix,
                       ostream_iterator<byte_type>(std::cout, ""));
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
            typedef map<char_type, id_type> char2id_t;
            char2id_t char2id;
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
