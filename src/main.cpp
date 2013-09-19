#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/function.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/optional/optional.hpp>
#include <boost/timer/timer.hpp>
#include "pstade/oven/algorithm.hpp"
#include "pstade/oven/converted.hpp"
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

#include "atm/branching_substrings.hpp"
#include "atm/blumer_substrings.hpp"
#include "atm/purity_maximal_substrings.hpp"
#include "atm/substrings.hpp"
#include "atm/substrings_from_longest.hpp"
#include "atm/coarse_substrings.hpp"
#include "atm/coarse_ngrams.hpp"
#include "atm/ngrams.hpp"
#include "atm/segments.hpp"
#include "atm/single_range.hpp"
#include "sast/sast.hpp"

#include "../../config.h"


namespace oven = pstade::oven;


using byte_type = std::uint8_t;
using index_type = std::int32_t;

template <class K, class V>
typename boost::function<V (K)> tr_by(const std::map<K, V>& m) {
    using namespace boost::lambda;

    struct translate {
        using result_type = typename std::map<K, V>::mapped_type;

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
        using result_type = typename std::vector<V>::value_type;

        result_type operator()(const std::vector<V>& v, const typename std::vector<V>::size_type& k) const {
            return v[k];
        }
    };

    return boost::lambda::bind(translate(), ref(v), _1);
}

using column_set_t = unsigned int;
constexpr unsigned int COLUMN_STRICT_PURITY      = 1 << 0;
constexpr unsigned int COLUMN_LOOSE_PURITY       = 1 << 1;
constexpr unsigned int COLUMN_LEFT_UNIVERSALITY  = 1 << 2;
constexpr unsigned int COLUMN_RIGHT_UNIVERSALITY = 1 << 3;

struct TsvResultPrinter {
    TsvResultPrinter(std::ostream& os, const column_set_t cs, bool show_all_pos, bool show_substr, bool escape)
        : os_(os),
          to_unicode_char_(),
          column_set_(cs),
          show_all_pos_(show_all_pos),
          show_substr_(show_substr),
          escape_(escape)
    {}

    template <class F>
    TsvResultPrinter(std::ostream& os, F to_unicode_char, const column_set_t cs, bool show_all_pos, bool show_substr, bool escape)
        : os_(os),
          to_unicode_char_(to_unicode_char),
          column_set_(cs),
          show_all_pos_(show_all_pos),
          show_substr_(show_substr),
          escape_(escape)
    {}

    void print_header() {
        os_ << (show_all_pos_ ? "positions" : "position")
            << "\t" << "length"
            << "\t" << "frequency";
        if (column_set_ & COLUMN_STRICT_PURITY)      os_ << "\t" << "s-purity";
        if (column_set_ & COLUMN_LOOSE_PURITY)       os_ << "\t" << "l-purity";
        if (column_set_ & COLUMN_LEFT_UNIVERSALITY)  os_ << "\t" << "l-universality";
        if (column_set_ & COLUMN_RIGHT_UNIVERSALITY) os_ << "\t" << "r-universality";
        if (show_substr_)                            os_ << "\t" << "substring";
        os_ << "\n";
    }

    void print_footer() {
    }

    template <class S>
    void print(const S& substr) {
        if (show_all_pos_) {
            const auto ps = substr.allpos();
            os_ << ps[0];
            for (size_t i = 1; i < ps.size(); ++i) {
                os_ << ',' << ps[i];
            }
        }
        else {
            os_ << substr.pos();
        }
        print_rest(substr);
    }

private:
    template <class S>
    void print_rest(const S& substr) {
        using boost::lambda::_1;

        os_ << "\t" << substr.length()
            << "\t" << substr.frequency();
        if (column_set_ & COLUMN_STRICT_PURITY)      os_ << "\t" << substr.spurity();
        if (column_set_ & COLUMN_LOOSE_PURITY)       os_ << "\t" << substr.lpurity();
        if (column_set_ & COLUMN_LEFT_UNIVERSALITY)  os_ << "\t" << substr.luniversality();
        if (column_set_ & COLUMN_RIGHT_UNIVERSALITY) os_ << "\t" << substr.runiversality();
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
            if (c == '"')  os_ << "\\\"";
            else           os_ << c;
        }
        os_ << '"';
    }

    template <class S>
    void print_escaped_substr(const S& substr) {
        os_ << '"';
        for (const auto& c : substr) {
            if      (c == '"')   os_ << "\\\"";
            else if (c == '\\')  os_ << "\\\\";
            else if (c == '\r')  os_ << "\\r";
            else if (c == '\n')  os_ << "\\n";
            else if (c == '\t')  os_ << "\\t";
            else                 os_ << c;
        }
        os_ << '"';
    }

private:
    using unicode_char_type = std::uint32_t;
    using largest_id_type = std::uint32_t;
    std::ostream& os_;
    boost::optional<boost::function<unicode_char_type (largest_id_type)>> to_unicode_char_;
    column_set_t column_set_;
    bool show_all_pos_;
    bool show_substr_;
    bool escape_;
};

struct JsonResultPrinter {
    JsonResultPrinter(std::ostream& os, const column_set_t cs, bool show_all_pos, bool show_substr)
        : os_(os),
          to_unicode_char_(),
          column_set_(cs),
          show_all_pos_(show_all_pos),
          show_substr_(show_substr),
          first_element_(true)
    {}

    template <class F>
    JsonResultPrinter(std::ostream& os, F to_unicode_char, const column_set_t cs, bool show_all_pos, bool show_substr)
        : os_(os),
          to_unicode_char_(to_unicode_char),
          column_set_(cs),
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
        if (first_element_) {
            first_element_ = false;
        }
        else {
            os_ << ",\n";
        }

        if (show_all_pos_) {
            os_ << "    { "
                << "\"positions\": [";

            const auto ps = substr.allpos();
            os_ << ps[0];
            for (size_t i = 1; i < ps.size(); ++i) {
                os_ << ", " << ps[i];
            }
            os_ << "]";
        }
        else {
            os_ << "    { "
                << "\"position\": " << substr.pos();
        }
        print_rest(substr);
    }

private:
    template <class S>
    void print_rest(const S& substr) {
        using boost::lambda::_1;

        os_ << ", " << "\"length\": "    << substr.length()
            << ", " << "\"frequency\": " << substr.frequency();
        if (column_set_ & COLUMN_STRICT_PURITY)      os_ << ", " << "\"strict_purity\": "      << substr.spurity();
        if (column_set_ & COLUMN_LOOSE_PURITY)       os_ << ", " << "\"loose_purity\": "       << substr.lpurity();
        if (column_set_ & COLUMN_LEFT_UNIVERSALITY)  os_ << ", " << "\"left_universality\": "  << substr.luniversality();
        if (column_set_ & COLUMN_RIGHT_UNIVERSALITY) os_ << ", " << "\"right_universality\": " << substr.runiversality();
        if (show_substr_) {
            os_ << ", " << "\"substring\": ";
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
            else if (c <= '\u001F') {
                os_ << "\\u" << std::setfill('0') << std::setw(4) << std::hex << static_cast<int>(c) << std::dec;
            }
            else  os_ << c;
        }
        os_ << '"';
    }

private:
    using unicode_char_type = std::uint32_t;
    using largest_id_type = std::uint32_t;
    std::ostream& os_;
    boost::optional<boost::function<unicode_char_type (largest_id_type)>> to_unicode_char_;
    column_set_t column_set_;
    bool show_all_pos_;
    bool show_substr_;
    bool first_element_;
};

struct BenchmarkPrinter {
    BenchmarkPrinter(std::ostream& os, const column_set_t cs, bool, bool)
        : os_(os),
          column_set_(cs),
          timer_()
    {}

    template <class F>
    BenchmarkPrinter(std::ostream& os, F, const column_set_t cs, bool, bool)
        : os_(os),
          column_set_(cs),
          timer_()
    {}

    void print_header() {
        timer_.start();
    }

    void print_footer() {
        timer_.stop();
        os_ << timer_.format(boost::timer::default_places, "%w\t%u\t%s\t%t\t%p") << std::endl;
    }

    template <class S>
    void print(const S& substr) {
        if (column_set_ & COLUMN_STRICT_PURITY)      substr.spurity();
        if (column_set_ & COLUMN_LOOSE_PURITY)       substr.lpurity();
        if (column_set_ & COLUMN_LEFT_UNIVERSALITY)  substr.luniversality();
        if (column_set_ & COLUMN_RIGHT_UNIVERSALITY) substr.runiversality();
    }

private:
    std::ostream& os_;
    column_set_t column_set_;
    boost::timer::cpu_timer timer_;
};

enum class PurityType {
    StrictPurity,
    LoosePurity
};

enum class EnumerationType {
    BlumerEnumeration,
    PurityMaximalEnumeration,
    BranchingEnumeration,
    FrequentEnumeration,
    LongestEnumeration,
    CoarseEnumeration,
    SegmentEnumeration,
    NGramEnumeration,
    CoarseNGramEnumeration,
    SingleRangeEnumeration
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
void do_rest_of_binary_mode(const std::size_t& alphabet_size, std::ifstream& is, ResultPrinter& printer, EnumerationType enum_type, const int resolution, const int ngram, std::pair<int, int> range, const substring_constraint& constraint);

template<class ID, class ResultPrinter>
void do_rest_of_text_mode(const std::size_t& alphabet_size, const std::vector<std::uint32_t>& id2char, std::ifstream& is, ResultPrinter& printer, EnumerationType enum_type, const int resolution, const int ngram, std::pair<int, int> range, const substring_constraint& constraint);

int main(int argc, char* argv[]) {
    using namespace std;

    // command line
    cmdline::parser p;
    p.add("help", 'h', "");
    p.add("version", 'V', "");
    p.add<string>("number-format", 'F', "", false, "fixed", cmdline::oneof<string>("fixed", "scientific"));
    p.add<string>("format", 0, "", false, "tsv", cmdline::oneof<string>("tsv", "json", "benchmark"));
    p.add<string>("mode", 'm', "", false, "binary", cmdline::oneof<string>("binary", "text"));
    p.add("strict-purity", 0, "");
    p.add("loose-purity", 0, "");
    p.add("left-universality", 0, "");
    p.add("right-universality", 0, "");
    p.add("show-all-positions", 'a', "");
    p.add("show-substring", 's', "");
    p.add("escape", 'E', "");
    p.add<string>("purity", 'p', "", false, "strict", cmdline::oneof<string>("strict", "loose"));
    p.add<string>("enum", 0, "", false, "frequent", cmdline::oneof<string>("blumer", "purity-maximal", "branching", "frequent", "longest", "coarse", "segment", "ngram", "coarse-ngram", "single-range"));
    p.add<int>("resolution", 'r', "", false, 1);
    p.add<int>("ngram", 'n', "", false, 1);
    p.add<int>("range-begin", 'b', "inclusive", false, 0);
    p.add<int>("range-end", 'e', "inclusive", false, -1);
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
    const EnumerationType enum_type = p.get<string>("enum") == "blumer"         ? EnumerationType::BlumerEnumeration
                                    : p.get<string>("enum") == "purity-maximal" ? EnumerationType::PurityMaximalEnumeration
                                    : p.get<string>("enum") == "branching"      ? EnumerationType::BranchingEnumeration
                                    : p.get<string>("enum") == "frequent"       ? EnumerationType::FrequentEnumeration
                                    : p.get<string>("enum") == "longest"        ? EnumerationType::LongestEnumeration
                                    : p.get<string>("enum") == "coarse"         ? EnumerationType::CoarseEnumeration
                                    : p.get<string>("enum") == "segment"        ? EnumerationType::SegmentEnumeration
                                    : p.get<string>("enum") == "ngram"          ? EnumerationType::NGramEnumeration
                                    : p.get<string>("enum") == "coarse-ngram"   ? EnumerationType::CoarseNGramEnumeration
                                    : p.get<string>("enum") == "single-range"   ? EnumerationType::SingleRangeEnumeration
                                    : throw runtime_error("Invalid enumeration type was specified.");
    const int resolution = p.get<int>("resolution");
    const int ngram = p.get<int>("ngram");
    const pair<int, int> range = make_pair(p.get<int>("range-begin"), p.get<int>("range-end"));

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

    // column set
    const column_set_t cs = (p.exist("strict-purity")      ? COLUMN_STRICT_PURITY      : 0)
                          | (p.exist("loose-purity")       ? COLUMN_LOOSE_PURITY       : 0)
                          | (p.exist("left-universality")  ? COLUMN_LEFT_UNIVERSALITY  : 0)
                          | (p.exist("right-universality") ? COLUMN_RIGHT_UNIVERSALITY : 0);

    if (p.get<string>("mode") == "binary") {
        // alphabets
        const size_t alphabet_size = 0x100;

        if (p.get<string>("format") == "tsv") {
            TsvResultPrinter printer(std::cout, cs, p.exist("show-all-positions"), p.exist("show-substring"), p.exist("escape"));
            do_rest_of_binary_mode(alphabet_size, is, printer, enum_type, resolution, ngram, range, constraint);
        }
        else if (p.get<string>("format") == "json") {
            JsonResultPrinter printer(std::cout, cs, p.exist("show-all-positions"), p.exist("show-substring"));
            do_rest_of_binary_mode(alphabet_size, is, printer, enum_type, resolution, ngram, range, constraint);
        }
        else if (p.get<string>("format") == "benchmark") {
            BenchmarkPrinter printer(std::cout, cs, p.exist("show-all-positions"), p.exist("show-substring"));
            do_rest_of_binary_mode(alphabet_size, is, printer, enum_type, resolution, ngram, range, constraint);
        }
        else {
            throw runtime_error("Unsupported output format is specified.");
        }
    }
    else {
        using char_type = std::uint32_t;

        // alphabets
        is.seekg(0);
        const set<char_type> alphabets = oven::streambuf_read(is) | oven::memoized | oven::utf8_decoded | oven::copied;
        const size_t alphabet_size = alphabets.size();

        // map: id -> char
        const vector<char_type> id2char = alphabets | oven::copied;

        if (p.get<string>("format") == "tsv") {
            TsvResultPrinter printer(std::cout, tr_by(id2char), cs, p.exist("show-all-positions"), p.exist("show-substring"), p.exist("escape"));
            if (alphabet_size <= 0x100) {
                do_rest_of_text_mode<std::uint8_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, range, constraint);
            }
            else if (alphabet_size <= 0x10000) {
                do_rest_of_text_mode<std::uint16_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, range, constraint);
            }
            else {
                do_rest_of_text_mode<std::uint32_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, range, constraint);
            }
        }
        else if (p.get<string>("format") == "json") {
            JsonResultPrinter printer(std::cout, tr_by(id2char), cs, p.exist("show-all-positions"), p.exist("show-substring"));
            if (alphabet_size <= 0x100) {
                do_rest_of_text_mode<std::uint8_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, range, constraint);
            }
            else if (alphabet_size <= 0x10000) {
                do_rest_of_text_mode<std::uint16_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, range, constraint);
            }
            else {
                do_rest_of_text_mode<std::uint32_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, range, constraint);
            }
        }
        else if (p.get<string>("format") == "benchmark") {
            BenchmarkPrinter printer(std::cout, tr_by(id2char), cs, p.exist("show-all-positions"), p.exist("show-substring"));
            if (alphabet_size <= 0x100) {
                do_rest_of_text_mode<std::uint8_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, range, constraint);
            }
            else if (alphabet_size <= 0x10000) {
                do_rest_of_text_mode<std::uint16_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, range, constraint);
            }
            else {
                do_rest_of_text_mode<std::uint32_t>(alphabet_size, id2char, is, printer, enum_type, resolution, ngram, range, constraint);
            }
        }
        else {
            throw runtime_error("Unsupported output format is specified.");
        }
    }
}

template <class ResultPrinter>
void do_rest_of_binary_mode(const std::size_t& alphabet_size, std::ifstream& is, ResultPrinter& printer, EnumerationType enum_type, const int resolution, const int ngram, std::pair<int, int> range, const substring_constraint& constraint)
{
    using namespace std;

    using id_type = std::uint8_t;

    // input
    is.seekg(0);
    const vector<id_type> input = oven::streambuf_read(is) | oven::converted<id_type>() | oven::copied;

    // sast
    sast::sast<decltype(input), index_type> sast(input, alphabet_size);

    // enumerate substrings
    printer.print_header();
    switch (enum_type) {
    case EnumerationType::BlumerEnumeration: {
        using substr_type = typename atm::blumer_substrings<decltype(input), index_type>::substr;

        atm::blumer_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::PurityMaximalEnumeration: {
        using substr_type = typename atm::purity_maximal_substrings<decltype(input), index_type>::substr;

        atm::purity_maximal_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::BranchingEnumeration: {
        using substr_type = typename atm::branching_substrings<decltype(input), index_type>::substr;

        atm::branching_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::FrequentEnumeration: {
        using substr_type = typename atm::substrings<decltype(input), index_type>::substr;

        atm::substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::LongestEnumeration: {
        using substr_type = typename atm::substrings_from_longest<decltype(input), index_type>::substr;

        atm::substrings_from_longest<decltype(input), index_type> substrs(sast);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::CoarseEnumeration: {
        using substr_type = typename atm::coarse_substrings<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }

        atm::coarse_substrings<decltype(input), index_type> substrs(sast, resolution);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::SegmentEnumeration: {
        using substr_type = typename atm::segments<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }

        atm::segments<decltype(input), index_type> substrs(sast, resolution);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::NGramEnumeration: {
        using substr_type = typename atm::ngrams<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(ngram) > input.size()) {
            throw runtime_error("Specified N for N-grams is out of the size of the input.");
        }

        atm::ngrams<decltype(input), index_type> substrs(sast, ngram);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::CoarseNGramEnumeration: {
        using substr_type = typename atm::coarse_ngrams<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }
        if (static_cast<std::size_t>(ngram) > input.size()) {
            throw runtime_error("Specified N for N-grams is out of the size of the input.");
        }

        atm::coarse_ngrams<decltype(input), index_type> substrs(sast, resolution, ngram);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::SingleRangeEnumeration: {
        using substr_type = typename atm::single_range<decltype(input), index_type>::substr;

        while (range.first  < 0) range.first  += input.size();
        while (range.second < 0) range.second += input.size();

        if (static_cast<std::size_t>(range.first) > input.size()) {
            throw runtime_error("Specified range beginning position is out of the size of the input.");
        }
        if (static_cast<std::size_t>(range.second) > input.size()) {
            throw runtime_error("Specified range ending position is out of the size of the input.");
        }
        if (range.first > range.second) {
            throw runtime_error("Specified range is invalid.");
        }


        atm::single_range<decltype(input), index_type> substrs(sast, range.first, range.second + 1);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    }
    printer.print_footer();
}

template<class ID, class ResultPrinter>
void do_rest_of_text_mode(const std::size_t& alphabet_size, const std::vector<std::uint32_t>& id2char, std::ifstream& is, ResultPrinter& printer, EnumerationType enum_type, const int resolution, const int ngram, std::pair<int, int> range, const substring_constraint& constraint)
{
    using namespace std;

    using char_type = std::uint32_t;
    using id_type = ID;

    // map: char -> id
    map<char_type, id_type> char2id;
    for (size_t id = 0; id < alphabet_size; ++id) {
        char2id[id2char[id]] = id;
    }

    // input
    is.seekg(0);
    const vector<id_type> input = oven::streambuf_read(is) | oven::memoized | oven::utf8_decoded | oven::transformed(tr_by(char2id)) | oven::copied;

    // sast
    sast::sast<decltype(input), index_type> sast(input, alphabet_size);

    // enumerate substrings
    printer.print_header();
    switch (enum_type) {
    case EnumerationType::BlumerEnumeration: {
        using substr_type = typename atm::blumer_substrings<decltype(input), index_type>::substr;

        atm::blumer_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::PurityMaximalEnumeration: {
        using substr_type = typename atm::purity_maximal_substrings<decltype(input), index_type>::substr;

        atm::purity_maximal_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::BranchingEnumeration: {
        using substr_type = typename atm::branching_substrings<decltype(input), index_type>::substr;

        atm::branching_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::FrequentEnumeration: {
        using substr_type = typename atm::substrings<decltype(input), index_type>::substr;

        atm::substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::LongestEnumeration: {
        using substr_type = typename atm::substrings_from_longest<decltype(input), index_type>::substr;

        atm::substrings_from_longest<decltype(input), index_type> substrs(sast);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::CoarseEnumeration: {
        using substr_type = typename atm::coarse_substrings<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }

        atm::coarse_substrings<decltype(input), index_type> substrs(sast, resolution);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::SegmentEnumeration: {
        using substr_type = typename atm::segments<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }

        atm::segments<decltype(input), index_type> substrs(sast, resolution);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::NGramEnumeration: {
        using substr_type = typename atm::ngrams<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(ngram) > input.size()) {
            throw runtime_error("Specified N for N-grams is out of the size of the input.");
        }

        atm::ngrams<decltype(input), index_type> substrs(sast, ngram);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::CoarseNGramEnumeration: {
        using substr_type = typename atm::coarse_ngrams<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }
        if (static_cast<std::size_t>(ngram) > input.size()) {
            throw runtime_error("Specified N for N-grams is out of the size of the input.");
        }

        atm::coarse_ngrams<decltype(input), index_type> substrs(sast, resolution, ngram);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case EnumerationType::SingleRangeEnumeration: {
        using substr_type = typename atm::single_range<decltype(input), index_type>::substr;

        while (range.first  < 0) range.first  += input.size();
        while (range.second < 0) range.second += input.size();

        if (static_cast<std::size_t>(range.first) > input.size()) {
            throw runtime_error("Specified range beginning position is out of the size of the input.");
        }
        if (static_cast<std::size_t>(range.second) > input.size()) {
            throw runtime_error("Specified range ending position is out of the size of the input.");
        }
        if (range.first > range.second) {
            throw runtime_error("Specified range is invalid.");
        }

        atm::single_range<decltype(input), index_type> substrs(sast, range.first, range.second + 1);
        for (auto substr : oven::make_filtered(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    }
    printer.print_footer();
}
