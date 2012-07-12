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
#include "substrings.hpp"


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

    // turn off the synchronization of iostream and cstdio.
    ios::sync_with_stdio(false);

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
        typedef typename Substrings<id_type, index_type>::substr substr_type;

        // alphabets
        const size_t alphabet_size = 0x100;

        // input
        const vector<id_type> input = is | oven::copied;

        // printer
        ResultPrinter printer(std::cout, p.exist("show-substring"), p.exist("exclude-newline"));

        // enumerate substrings
        Substrings<id_type, index_type> substrs(input, alphabet_size);

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
            typedef typename Substrings<id_type, index_type>::substr substr_type;

            // map: char -> id
            map<char_type, id_type> char2id;
            for (size_t id = 0; id < alphabet_size; ++id) {
                char2id[id2char[id]] = id;
            }

            // input
            const vector<id_type> input = is | oven::utf8_decoded | oven::transformed(tr_by(char2id)) | oven::copied;

            // enumerate substrings
            Substrings<id_type, index_type> substrs(input, alphabet_size);

            printer.print_header();
            for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
                printer.print(substr);
            }
        }
        else if (alphabet_size <= 0x10000) {
            typedef boost::uint16_t id_type;
            typedef typename Substrings<id_type, index_type>::substr substr_type;

            // map: char -> id
            map<char_type, id_type> char2id;
            for (size_t id = 0; id < alphabet_size; ++id) {
                char2id[id2char[id]] = id;
            }

            // input
            const vector<id_type> input = is | oven::utf8_decoded | oven::transformed(tr_by(char2id)) | oven::copied;

            // enumerate substrings
            Substrings<id_type, index_type> substrs(input, alphabet_size);

            printer.print_header();
            for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
                printer.print(substr);
            }
        }
        else {
            typedef boost::uint32_t id_type;
            typedef typename Substrings<id_type, index_type>::substr substr_type;

            // map: char -> id
            map<char_type, id_type> char2id;
            for (size_t id = 0; id < alphabet_size; ++id) {
                char2id[id2char[id]] = id;
            }

            // input
            const vector<id_type> input = is | oven::utf8_decoded | oven::transformed(tr_by(char2id)) | oven::copied;

            // enumerate substrings
            Substrings<id_type, index_type> substrs(input, alphabet_size);

            printer.print_header();
            for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
                printer.print(substr);
            }
        }
    }
}
