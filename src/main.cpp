#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/cstdint.hpp>
#include <boost/function.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/optional/optional.hpp>
#include "pstade/oven/algorithm.hpp"
#include "pstade/oven/copied.hpp"
#include "pstade/oven/filtered.hpp"
#include "pstade/oven/filterer.hpp"
#include "pstade/oven/make_range.hpp"
#include "pstade/oven/memoized.hpp"
#include "pstade/oven/stream_read.hpp"
#include "pstade/oven/transformed.hpp"
#include "pstade/oven/utf8_decoded.hpp"
#include "pstade/oven/utf8_encoded.hpp"
#include "cmdline.h"
#include "substrings.hpp"
#include "../../config.h"


namespace oven = pstade::oven;


typedef boost::uint8_t byte_type;
typedef boost::int32_t index_type;

template <class K, class V>
typename boost::function<V (K)> tr_by(const std::map<K, V>& m) {
    using namespace boost::lambda;

    struct translate {
        typedef typename std::map<K, V>::mapped_type result_type;

        result_type operator()(const std::map<K, V>& m, const typename std::map<K, V>::key_type& k) const {
            return m.at(k);
        }
    };

    return boost::lambda::bind(translate(), ref(m), _1);
}

template <class V>
typename boost::function<V (typename std::vector<V>::size_type)> tr_by(const std::vector<V>& v) {
    using namespace boost::lambda;

    struct translate {
        typedef typename std::vector<V>::value_type result_type;

        result_type operator()(const std::vector<V>& v, const typename std::vector<V>::size_type& k) const {
            return v[k];
        }
    };

    return boost::lambda::bind(translate(), ref(v), _1);
}

struct TsvResultPrinter {
    TsvResultPrinter(std::ostream& os, bool show_all_pos, bool show_substr, bool escape)
        : os_(os),
          to_unicode_char_(),
          show_all_pos_(show_all_pos),
          show_substr_(show_substr),
          escape_(escape)
    {}

    template <class F>
    TsvResultPrinter(std::ostream& os, F to_unicode_char, bool show_all_pos, bool show_substr, bool escape)
        : os_(os),
          to_unicode_char_(to_unicode_char),
          show_all_pos_(show_all_pos),
          show_substr_(show_substr),
          escape_(escape)
    {}

    void print_header() {
        os_ << "position"       << "\t"
            << "length"         << "\t"
            << "frequency"      << "\t"
            << "s-purity"       << "\t"
            << "l-purity"       << "\t"
            << "l-universality" << "\t"
            << "r-universality";
        if (show_substr_) {
            os_ << "\t" << "substring";
        }
        os_ << "\n";
    }

    void print_footer() {
    }

    template <class S>
    void print(const S& substr) {
        if (show_all_pos_) {
            const auto ps = substr.allpos();
            for (const auto p : ps) {
                os_ << p << "\t";
                print_rest(substr);
            }
        }
        else {
            os_ << substr.pos() << "\t";
            print_rest(substr);
        }
    }

private:
    template <class S>
    void print_rest(const S& substr) {
        using boost::lambda::_1;

        os_ << substr.length()         << "\t"
            << substr.frequency()      << "\t"
            << substr.spurity()        << "\t"
            << substr.lpurity()        << "\t"
            << substr.luniversality()  << "\t"
            << substr.runiversality();
        if (show_substr_) {
            os_ << "\t";
            if (to_unicode_char_) {
                auto encoded = substr | oven::transformed(*to_unicode_char_) | oven::utf8_encoded;
                if (escape_) {
                    print_escaped_substr(encoded);
                }
                else {
                    print_substr(encoded);
                }
            }
            else {
                if (escape_) {
                    print_escaped_substr(substr);
                }
                else {
                    print_substr(substr);
                }
            }
        }
        os_ << "\n";
    }

    template <class S>
    void print_substr(const S& substr) {
        os_ << '"';
        for (const auto& c : substr) {
            if (c == '"')  os_ << "\"\"";
            else           os_ << c;
        }
        os_ << '"';
    }

    template <class S>
    void print_escaped_substr(const S& substr) {
        os_ << '"';
        for (const auto& c : substr) {
            if      (c == '"')   os_ << "\"\"";
            else if (c == '\\')  os_ << "\\\\";
            else if (c == '\r')  os_ << "\\r";
            else if (c == '\n')  os_ << "\\n";
            else if (c == '\t')  os_ << "\\t";
            else                 os_ << c;
        }
        os_ << '"';
    }

private:
    typedef boost::uint32_t unicode_char_type;
    typedef boost::uint32_t largest_id_type;
    std::ostream& os_;
    boost::optional<boost::function<unicode_char_type (largest_id_type)>> to_unicode_char_;
    bool show_all_pos_;
    bool show_substr_;
    bool escape_;
};

struct JsonResultPrinter {
    JsonResultPrinter(std::ostream& os, bool show_all_pos, bool show_substr)
        : os_(os),
          to_unicode_char_(),
          show_all_pos_(show_all_pos),
          show_substr_(show_substr),
          first_element_(true)
    {}

    template <class F>
    JsonResultPrinter(std::ostream& os, F to_unicode_char, bool show_all_pos, bool show_substr)
        : os_(os),
          to_unicode_char_(to_unicode_char),
          show_all_pos_(show_all_pos),
          show_substr_(show_substr),
          first_element_(true)
    {}

    void print_header() {
        os_ << "{\n"
            << "  \"substrings\": [\n";
    }

    void print_footer() {
        if (!first_element_)  os_ << "\n";
        os_ << "  ]\n"
            << "}\n";
    }

    template <class S>
    void print(const S& substr) {
        if (show_all_pos_) {
            const auto ps = substr.allpos();
            for (const auto p : ps) {
                if (first_element_) {
                    first_element_ = false;
                }
                else {
                    os_ << ",\n";
                }

                os_ << "    { "
                    << "\"position\": " << p << ", ";
                print_rest(substr);
            }
        }
        else {
            if (first_element_) {
                first_element_ = false;
            }
            else {
                os_ << ",\n";
            }

            os_ << "    { "
                << "\"position\": " << substr.pos() << ", ";
            print_rest(substr);
        }
    }

private:
    template <class S>
    void print_rest(const S& substr) {
        using boost::lambda::_1;

        os_ << "\"length\": "             << substr.length()         << ", "
            << "\"frequency\": "          << substr.frequency()      << ", "
            << "\"strict_purity\": "      << substr.spurity()        << ", "
            << "\"loose_purity\": "       << substr.lpurity()        << ", "
            << "\"left_universality\": "  << substr.luniversality()  << ", "
            << "\"right_universality\": " << substr.runiversality();
        if (show_substr_) {
            os_ << ", "
                << "\"substring\": ";
            if (to_unicode_char_) {
                auto encoded = substr | oven::transformed(*to_unicode_char_) | oven::utf8_encoded;
                print_substr(encoded);
            }
            else {
                print_substr(substr);
            }
        }
        os_ << " }";
    }

    template <class S>
    void print_substr(const S& substr) {
        os_ << '"';
        for (const auto& c : substr) {
            if      (c == '"')       os_ << "\\\"";
            else if (c == '\\')      os_ << "\\\\";
            else if (c == '\u0008')  os_ << "\\b";
            else if (c == '\u000C')  os_ << "\\f";
            else if (c == '\u000A')  os_ << "\\n";
            else if (c == '\u000D')  os_ << "\\r";
            else if (c == '\u0009')  os_ << "\\t";
            else if (c >= '\u0000' && c <= '\u001F') {
                os_ << "\\u" << std::setfill('0') << std::setw(4) << std::hex << static_cast<int>(c) << std::dec;
            }
            else  os_ << c;
        }
        os_ << '"';
    }

private:
    typedef boost::uint32_t unicode_char_type;
    typedef boost::uint32_t largest_id_type;
    std::ostream& os_;
    boost::optional<boost::function<unicode_char_type (largest_id_type)>> to_unicode_char_;
    bool show_all_pos_;
    bool show_substr_;
    bool first_element_;
};

enum class PurityType {
    StrictPurity,
    LoosePurity
};

enum class EnumerationType {
    BranchingEnumeration,
    FrequentEnumeration,
    LongestEnumeration,
    CoarseEnumeration,
    SegmentEnumeration,
    NGramEnumeration
};

struct substring_constraint {
    PurityType purity_type;
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
            && (!c_.min_purity    || get_purity(substr) >= *c_.min_purity)
            && (!c_.max_purity    || get_purity(substr) <= *c_.max_purity);
    }

private:
    double get_purity(const Substring& substr) const {
        if (c_.purity_type == PurityType::StrictPurity) {
            return substr.spurity();
        }
        else {
            return substr.lpurity();
        }
    }

    const substring_constraint c_;
};

template <class ResultPrinter>
void do_rest_of_binary_mode(const std::size_t& alphabet_size, std::ifstream& is, ResultPrinter& printer, EnumerationType enum_type, const int resolution, const int ngram, const substring_constraint& constraint);

template<class Char, class ID, class ResultPrinter>
void do_rest_of_text_mode(const std::size_t& alphabet_size, const std::vector<Char>& id2char, std::ifstream& is, ResultPrinter& printer, EnumerationType enum_type, const int resolution, const int ngram, const substring_constraint& constraint);

int main(int argc, char* argv[]) {
    using namespace std;

    // command line
    cmdline::parser p;
    p.add("help", 'h', "");
    p.add("version", 'V', "");
    p.add<string>("number-format", 'F', "", false, "fixed", cmdline::oneof<string>("fixed", "scientific"));
    p.add<string>("format", 0, "", false, "tsv", cmdline::oneof<string>("tsv", "json"));
    p.add<string>("mode", 'm', "", false, "binary", cmdline::oneof<string>("binary", "text"));
    p.add("show-all-position", 'a', "");
    p.add("show-substring", 's', "");
    p.add("escape", 'e', "");
    p.add<string>("purity", 'p', "", false, "strict", cmdline::oneof<string>("strict", "loose"));
    p.add<string>("enum", 0, "", false, "frequent", cmdline::oneof<string>("branching", "frequent", "longest", "coarse", "segment", "ngram"));
    p.add<int>("resolution", 'r', "", false, 1);
    p.add<int>("ngram", 'n', "", false, 1);
    p.add<index_type>("longer",  0, "", false, -1);
    p.add<index_type>("shorter", 0, "", false, -1);
    p.add<index_type>("more-frequent",   0, "", false, -1);
    p.add<index_type>("more-infrequent", 0, "", false, -1);
    p.add<double>("purer",   0, "", false, -1);
    p.add<double>("impurer", 0, "", false, -1);
    if (!p.parse(argc, argv) || p.exist("help")) {
        cout << p.error_full() << p.usage();
        return 0;
    }
    else if (p.exist("version")) {
        cout << APP_NAME " " APP_VERSION << endl;
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
    ifstream is(fp);
    if (!is)  throw runtime_error("Failed to open the input file.");

    // number format
    cout.setf(p.get<string>("number-format") == "fixed" ? ios::fixed : ios::scientific,
              ios::floatfield);

    // substring enumeration type
    const EnumerationType enum_type = p.get<string>("enum") == "branching" ? EnumerationType::BranchingEnumeration
                                    : p.get<string>("enum") == "frequent"  ? EnumerationType::FrequentEnumeration
                                    : p.get<string>("enum") == "longest"   ? EnumerationType::LongestEnumeration
                                    : p.get<string>("enum") == "coarse"    ? EnumerationType::CoarseEnumeration
                                    : p.get<string>("enum") == "segment"   ? EnumerationType::SegmentEnumeration
                                    : p.get<string>("enum") == "ngram"     ? EnumerationType::NGramEnumeration
                                    : throw runtime_error("Invalid enumeration type was specified.");
    const int resolution = p.get<int>("resolution");
    const int ngram = p.get<int>("ngram");

    // substring filter
    substring_constraint constraint = {
        p.get<string>("purity") == "strict" ? PurityType::StrictPurity : PurityType::LoosePurity,
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

        // alphabets
        const size_t alphabet_size = 0x100;

        if (p.get<string>("format") == "tsv") {
            TsvResultPrinter printer(std::cout, p.exist("show-all-position"), p.exist("show-substring"), p.exist("escape"));
            do_rest_of_binary_mode(alphabet_size, is, printer, enum_type, resolution, ngram, constraint);
        }
        else if (p.get<string>("format") == "json") {
            JsonResultPrinter printer(std::cout, p.exist("show-all-position"), p.exist("show-substring"));
            do_rest_of_binary_mode(alphabet_size, is, printer, enum_type, resolution, ngram, constraint);
        }
        else {
            throw runtime_error("Unsupported output format is specified.");
        }
    }
    else {
        typedef boost::uint32_t char_type;

        // alphabets
        is.seekg(0);
        const set<char_type> alphabets = oven::streambuf_read(is) | oven::memoized | oven::utf8_decoded | oven::copied;
        const size_t alphabet_size = alphabets.size();

        // map: id -> char
        const vector<char_type> id2char = alphabets | oven::copied;

        if (p.get<string>("format") == "tsv") {
            TsvResultPrinter printer(std::cout, tr_by(id2char), p.exist("show-all-position"), p.exist("show-substring"), p.exist("escape"));
            if (alphabet_size <= 0x100) {
                do_rest_of_text_mode<char_type, boost::uint8_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, constraint);
            }
            else if (alphabet_size <= 0x10000) {
                do_rest_of_text_mode<char_type, boost::uint16_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, constraint);
            }
            else {
                do_rest_of_text_mode<char_type, boost::uint32_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, constraint);
            }
        }
        else if (p.get<string>("format") == "json") {
            JsonResultPrinter printer(std::cout, tr_by(id2char), p.exist("show-all-position"), p.exist("show-substring"));
            if (alphabet_size <= 0x100) {
                do_rest_of_text_mode<char_type, boost::uint8_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, constraint);
            }
            else if (alphabet_size <= 0x10000) {
                do_rest_of_text_mode<char_type, boost::uint16_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, constraint);
            }
            else {
                do_rest_of_text_mode<char_type, boost::uint32_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, constraint);
            }
        }
        else {
            throw runtime_error("Unsupported output format is specified.");
        }
    }
}

template <class ResultPrinter>
void do_rest_of_binary_mode(const std::size_t& alphabet_size, std::ifstream& is, ResultPrinter& printer, EnumerationType enum_type, const int resolution, const int ngram, const substring_constraint& constraint)
{
    using namespace std;

    typedef boost::uint8_t char_type;
    typedef boost::uint8_t id_type;

    // input
    is.seekg(0);
    const vector<id_type> input = oven::streambuf_read(is) | oven::copied;

    // enumerate substrings
    printer.print_header();
    switch (enum_type) {
    case EnumerationType::BranchingEnumeration: {
        typedef typename BranchingSubstrings<id_type, index_type>::substr substr_type;

        BranchingSubstrings<id_type, index_type> substrs(input, alphabet_size);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::FrequentEnumeration: {
        typedef typename Substrings<id_type, index_type>::substr substr_type;

        Substrings<id_type, index_type> substrs(input, alphabet_size);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::LongestEnumeration: {
        typedef typename SubstringsFromLongest<id_type, index_type>::substr substr_type;

        SubstringsFromLongest<id_type, index_type> substrs(input, alphabet_size);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::CoarseEnumeration: {
        typedef typename CoarseSubstrings<id_type, index_type>::substr substr_type;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }

        CoarseSubstrings<id_type, index_type> substrs(input, alphabet_size, resolution);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::SegmentEnumeration: {
        typedef typename Segments<id_type, index_type>::substr substr_type;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }

        Segments<id_type, index_type> substrs(input, alphabet_size, resolution);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::NGramEnumeration: {
        typedef typename NGrams<id_type, index_type>::substr substr_type;

        if (static_cast<std::size_t>(ngram) > input.size()) {
            throw runtime_error("Specified N for N-grams is out of the size of the input.");
        }

        NGrams<id_type, index_type> substrs(input, alphabet_size, ngram);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    }
    printer.print_footer();
}

template<class Char, class ID, class ResultPrinter>
void do_rest_of_text_mode(const std::size_t& alphabet_size, const std::vector<Char>& id2char, std::ifstream& is, ResultPrinter& printer, EnumerationType enum_type, const int resolution, const int ngram, const substring_constraint& constraint)
{
    using namespace std;

    typedef Char char_type;
    typedef ID id_type;

    // map: char -> id
    map<char_type, id_type> char2id;
    for (size_t id = 0; id < alphabet_size; ++id) {
        char2id[id2char[id]] = id;
    }

    // input
    is.seekg(0);
    const vector<id_type> input = oven::streambuf_read(is) | oven::memoized | oven::utf8_decoded | oven::transformed(tr_by(char2id)) | oven::copied;

    // enumerate substrings
    printer.print_header();
    switch (enum_type) {
    case EnumerationType::BranchingEnumeration: {
        typedef typename BranchingSubstrings<id_type, index_type>::substr substr_type;

        BranchingSubstrings<id_type, index_type> substrs(input, alphabet_size);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::FrequentEnumeration: {
        typedef typename Substrings<id_type, index_type>::substr substr_type;

        Substrings<id_type, index_type> substrs(input, alphabet_size);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::LongestEnumeration: {
        typedef typename SubstringsFromLongest<id_type, index_type>::substr substr_type;

        SubstringsFromLongest<id_type, index_type> substrs(input, alphabet_size);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::CoarseEnumeration: {
        typedef typename CoarseSubstrings<id_type, index_type>::substr substr_type;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }

        CoarseSubstrings<id_type, index_type> substrs(input, alphabet_size, resolution);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::SegmentEnumeration: {
        typedef typename Segments<id_type, index_type>::substr substr_type;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }

        Segments<id_type, index_type> substrs(input, alphabet_size, resolution);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::NGramEnumeration: {
        typedef typename NGrams<id_type, index_type>::substr substr_type;

        if (static_cast<std::size_t>(ngram) > input.size()) {
            throw runtime_error("Specified N for N-grams is out of the size of the input.");
        }

        NGrams<id_type, index_type> substrs(input, alphabet_size, ngram);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    }
    printer.print_footer();
}
