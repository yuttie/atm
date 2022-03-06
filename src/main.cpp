#include <cctype>
#include <cstdint>
#include <codecvt>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <locale>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <boost/function.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/optional/optional.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/iterator_range_core.hpp>
#include <boost/timer/timer.hpp>
#include "nlohmann/json.hpp"
#include "cmdline.h"

#include "atm/branching_substrings.hpp"
#include "atm/blumer_substrings.hpp"
#include "atm/purity_maximal_substrings.hpp"
#include "atm/substrings.hpp"
#include "atm/substrings_from_longest.hpp"
#include "atm/coarse_substrings.hpp"
#include "atm/coarse_ngrams.hpp"
#include "atm/ngrams.hpp"
#include "atm/chunked_sliding_window.hpp"
#include "atm/words.hpp"
#include "atm/single_range.hpp"
#include "sast/sast.hpp"

#include "config.h"


using json = nlohmann::json;


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
    TsvResultPrinter(std::ostream& os, const column_set_t cs, bool show_header, bool show_all_pos, bool show_substr, bool escape)
        : os_(os),
          to_unicode_char_(),
          column_set_(cs),
          show_header_(show_header),
          show_all_pos_(show_all_pos),
          show_substr_(show_substr),
          escape_(escape)
    {}

    template <class F>
    TsvResultPrinter(std::ostream& os, F to_unicode_char, const column_set_t cs, bool show_header, bool show_all_pos, bool show_substr, bool escape)
        : os_(os),
          to_unicode_char_(to_unicode_char),
          column_set_(cs),
          show_header_(show_header),
          show_all_pos_(show_all_pos),
          show_substr_(show_substr),
          escape_(escape)
    {}

    void print_header() {
        if (show_header_) {
            os_ << "position"
                << "\t" << "length"
                << "\t" << "frequency";
            if (column_set_ & COLUMN_STRICT_PURITY)      os_ << "\t" << "s-purity";
            if (column_set_ & COLUMN_LOOSE_PURITY)       os_ << "\t" << "l-purity";
            if (column_set_ & COLUMN_LEFT_UNIVERSALITY)  os_ << "\t" << "l-universality";
            if (column_set_ & COLUMN_RIGHT_UNIVERSALITY) os_ << "\t" << "r-universality";
            if (show_substr_)                            os_ << "\t" << "substring";
            os_ << "\n";
        }
    }

    void print_footer() {
    }

    template <class S>
    void print(const S& substr) {
        if (show_all_pos_) {
            const auto ps = substr.allpos();
            for (size_t i = 0; i < ps.size(); ++i) {
                os_ << ps[i];
                print_rest(substr);
            }
        }
        else {
            os_ << substr.pos();
            print_rest(substr);
        }
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
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                std::string encoded = cvt.to_bytes(
                        boost::copy_range<std::u32string>(
                            substr | boost::adaptors::transformed(*to_unicode_char_)));
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
    bool show_header_;
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
        if (show_all_pos_) {
            const auto ps = substr.allpos();
            for (size_t i = 0; i < ps.size(); ++i) {
                if (first_element_) {
                    first_element_ = false;
                }
                else {
                    os_ << ",\n";
                }

                os_ << "    { "
                    << "\"position\": " << ps[i];
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
                << "\"position\": " << substr.pos();
            print_rest(substr);
        }
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
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                std::string encoded = cvt.to_bytes(
                        boost::copy_range<std::u32string>(
                            substr | boost::adaptors::transformed(*to_unicode_char_)));
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

struct JsonLinesResultPrinter {
    JsonLinesResultPrinter(std::ostream& os, const column_set_t cs, bool show_all_pos, bool show_substr)
        : os_(os),
          to_unicode_char_(),
          column_set_(cs),
          show_all_pos_(show_all_pos),
          show_substr_(show_substr)
    {}

    template <class F>
    JsonLinesResultPrinter(std::ostream& os, F to_unicode_char, const column_set_t cs, bool show_all_pos, bool show_substr)
        : os_(os),
          to_unicode_char_(to_unicode_char),
          column_set_(cs),
          show_all_pos_(show_all_pos),
          show_substr_(show_substr)
    {}

    void print_header() {}

    void print_footer() {}

    template <class S>
    void print(const S& substr) {
        json j;
        j["position"]  = -1;
        j["length"]    = substr.length();
        j["frequency"] = substr.frequency();
        if (column_set_ & COLUMN_STRICT_PURITY)      j["strict_purity"]      = substr.spurity();
        if (column_set_ & COLUMN_LOOSE_PURITY)       j["loose_purity"]       = substr.lpurity();
        if (column_set_ & COLUMN_LEFT_UNIVERSALITY)  j["left_universality"]  = substr.luniversality();
        if (column_set_ & COLUMN_RIGHT_UNIVERSALITY) j["right_universality"] = substr.runiversality();
        if (show_substr_) {
            if (to_unicode_char_) {
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                std::string encoded = cvt.to_bytes(
                        boost::copy_range<std::u32string>(
                            substr | boost::adaptors::transformed(*to_unicode_char_)));
                j["substring"] = encoded;
            }
            else {
                j["substring"] = boost::copy_range<std::string>(substr);
            }
        }

        if (show_all_pos_) {
            const auto ps = substr.allpos();
            for (size_t i = 0; i < ps.size(); ++i) {
                j["position"] = ps[i];
                os_ << j.dump() << std::endl;
            }
        }
        else {
            j["position"] = substr.pos();
            os_ << j.dump() << std::endl;
        }
    }

private:
    using unicode_char_type = std::uint32_t;
    using largest_id_type = std::uint32_t;
    std::ostream& os_;
    boost::optional<boost::function<unicode_char_type (largest_id_type)>> to_unicode_char_;
    column_set_t column_set_;
    bool show_all_pos_;
    bool show_substr_;
};

struct NumberArrayJsonLinesResultPrinter {
    NumberArrayJsonLinesResultPrinter(std::ostream& os, const column_set_t cs, bool show_all_pos, bool show_substr)
        : os_(os),
          column_set_(cs),
          show_all_pos_(show_all_pos),
          show_substr_(show_substr)
    {}

    void print_header() {}

    void print_footer() {}

    template <class S>
    void print(const S& substr) {
        using id_type = typename std::iterator_traits<typename S::const_iterator>::value_type;

        json j;
        j["position"]  = -1;
        j["length"]    = substr.length();
        j["frequency"] = substr.frequency();
        if (column_set_ & COLUMN_STRICT_PURITY)      j["strict_purity"]      = substr.spurity();
        if (column_set_ & COLUMN_LOOSE_PURITY)       j["loose_purity"]       = substr.lpurity();
        if (column_set_ & COLUMN_LEFT_UNIVERSALITY)  j["left_universality"]  = substr.luniversality();
        if (column_set_ & COLUMN_RIGHT_UNIVERSALITY) j["right_universality"] = substr.runiversality();
        if (show_substr_) {
            j["substring"] = boost::copy_range<std::vector<id_type>>(substr);
        }

        if (show_all_pos_) {
            const auto ps = substr.allpos();
            for (size_t i = 0; i < ps.size(); ++i) {
                j["position"] = ps[i];
                os_ << j.dump() << std::endl;
            }
        }
        else {
            j["position"] = substr.pos();
            os_ << j.dump() << std::endl;
        }
    }

private:
    std::ostream& os_;
    column_set_t column_set_;
    bool show_all_pos_;
    bool show_substr_;
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

enum class Enumeration {
    BlumerStrings,
    PurityMaximalStrings,
    BranchingStrings,
    NonLeafStrings,
    AllSubStrings,
    AllBlockwiseSubStrings,
    ChunkedSlidingWindow,
    SlidingWindow,
    BlockwiseSlidingWindow,
    Words,
    SingleSubstring
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
void do_rest_of_binary_mode(const std::size_t& alphabet_size, const std::vector<std::uint8_t> input, ResultPrinter& printer, Enumeration enum_type, const int resolution, const int block, const int window, std::pair<int, int> range, const substring_constraint& constraint);

template<class ID, class ResultPrinter>
void do_rest_of_text_mode(const std::size_t& alphabet_size, const std::vector<std::uint32_t>& id2char, const std::vector<ID> input, ResultPrinter& printer, Enumeration enum_type, const int resolution, const int block, const int window, std::pair<int, int> range, const substring_constraint& constraint);

template<class ID, class ResultPrinter>
void do_rest_of_json_number_array_mode(const std::size_t& alphabet_size, const std::vector<ID> input, ResultPrinter& printer, Enumeration enum_type, const int resolution, const int block, const int window, std::pair<int, int> range, const substring_constraint& constraint);

int main(int argc, char* argv[]) {
    using namespace std;

    // command line
    cmdline::parser p;
    p.add("help", 'h', "");
    p.add("version", 'V', "");
    p.add<string>("number-format", 'F',
                  "one of: fixed, scientific",
                  false, "fixed",
                  cmdline::oneof<string>("fixed", "scientific"));
    p.add<string>("format", 0,
                  "one of: tsv, json, json-lines, benchmark",
                  false, "tsv",
                  cmdline::oneof<string>("tsv", "json", "json-lines", "benchmark"));
    p.add<string>("mode", 'm', "one of: binary, text, json-number-array",
                  false, "binary",
                  cmdline::oneof<string>("binary", "text", "json-number-array"));
    p.add("header", 0, "");
    p.add("strict-purity", 0, "");
    p.add("loose-purity", 0, "");
    p.add("left-universality", 0, "");
    p.add("right-universality", 0, "");
    p.add("show-all-positions", 'a', "");
    p.add("show-substring", 's', "");
    p.add("escape", 'E', "");
    p.add<string>("purity", 'p', "one of: strict, loose",
                  false, "strict",
                  cmdline::oneof<string>("strict", "loose"));
    p.add<string>("enum", 0, "one of: blumer, purity-maximal, branching, non-leaf, all, all-blockwise, chunked-sliding-window, sliding-window, blockwise-sliding-window, words, single-substring",
                  false, "non-leaf",
                  cmdline::oneof<string>("blumer", "purity-maximal", "branching", "non-leaf", "all", "all-blockwise", "chunked-sliding-window", "sliding-window", "blockwise-sliding-window", "words", "single-substring"));
    p.add<int>("resolution", 'r', "", false, 1);
    p.add<int>("block", 0, "", false, 1);
    p.add<int>("window", 'w', "", false, 1);
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
    const Enumeration enum_type = p.get<string>("enum") == "blumer"                    ? Enumeration::BlumerStrings
                                : p.get<string>("enum") == "purity-maximal"            ? Enumeration::PurityMaximalStrings
                                : p.get<string>("enum") == "branching"                 ? Enumeration::BranchingStrings
                                : p.get<string>("enum") == "non-leaf"                  ? Enumeration::NonLeafStrings
                                : p.get<string>("enum") == "all"                       ? Enumeration::AllSubStrings
                                : p.get<string>("enum") == "all-blockwise"             ? Enumeration::AllBlockwiseSubStrings
                                : p.get<string>("enum") == "chunked-sliding-window"    ? Enumeration::ChunkedSlidingWindow
                                : p.get<string>("enum") == "sliding-window"            ? Enumeration::SlidingWindow
                                : p.get<string>("enum") == "blockwise-sliding-window"  ? Enumeration::BlockwiseSlidingWindow
                                : p.get<string>("enum") == "words"                     ? Enumeration::Words
                                : p.get<string>("enum") == "single-substring"          ? Enumeration::SingleSubstring
                                : throw runtime_error("Invalid enumeration type was specified.");
    const int resolution = p.get<int>("resolution");
    const int block = p.get<int>("block");
    const int window = p.get<int>("window");
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
        // input
        is.seekg(0);
        const auto input = vector<std::uint8_t>(std::istreambuf_iterator<char>(is), {});

        // alphabets
        const size_t alphabet_size = 0x100;

        if (p.get<string>("format") == "tsv") {
            TsvResultPrinter printer(std::cout, cs, p.exist("header"),p.exist("show-all-positions"), p.exist("show-substring"), p.exist("escape"));
            do_rest_of_binary_mode(alphabet_size, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
        }
        else if (p.get<string>("format") == "json") {
            JsonResultPrinter printer(std::cout, cs, p.exist("show-all-positions"), p.exist("show-substring"));
            do_rest_of_binary_mode(alphabet_size, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
        }
        else if (p.get<string>("format") == "json-lines") {
            JsonResultPrinter printer(std::cout, cs, p.exist("show-all-positions"), p.exist("show-substring"));
            do_rest_of_binary_mode(alphabet_size, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
        }
        else if (p.get<string>("format") == "benchmark") {
            BenchmarkPrinter printer(std::cout, cs, p.exist("show-all-positions"), p.exist("show-substring"));
            do_rest_of_binary_mode(alphabet_size, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
        }
        else {
            throw runtime_error("Unsupported output format is specified.");
        }
    }
    else if (p.get<string>("mode") == "text") {
        using char_type = std::uint32_t;

        // alphabets
        is.seekg(0);
        std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
        const auto alphabets = boost::copy_range<set<char_type>>(
                cvt.from_bytes(
                    std::string(
                        std::istreambuf_iterator<char>(is),
                        {})));
        const size_t alphabet_size = alphabets.size();

        // map: id -> char
        const auto id2char = boost::copy_range<vector<char_type>>(alphabets);

        if (p.get<string>("format") == "tsv") {
            TsvResultPrinter printer(std::cout, tr_by(id2char), cs, p.exist("header"), p.exist("show-all-positions"), p.exist("show-substring"), p.exist("escape"));
            if (alphabet_size <= 0x100) {
                using id_type = std::uint8_t;
                // map: char -> id
                map<char_type, id_type> char2id;
                for (size_t id = 0; id < alphabet_size; ++id) { char2id[id2char[id]] = id; }
                // input
                is.seekg(0);
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                const auto input = boost::copy_range<vector<id_type>>(
                        cvt.from_bytes(std::string(std::istreambuf_iterator<char>(is), {})) | boost::adaptors::transformed(tr_by(char2id)));
                // Do the rest
                do_rest_of_text_mode<std::uint8_t>(alphabet_size, id2char, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
            }
            else if (alphabet_size <= 0x10000) {
                using id_type = std::uint16_t;
                // map: char -> id
                map<char_type, id_type> char2id;
                for (size_t id = 0; id < alphabet_size; ++id) { char2id[id2char[id]] = id; }
                // input
                is.seekg(0);
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                const auto input = boost::copy_range<vector<id_type>>(
                        cvt.from_bytes(std::string(std::istreambuf_iterator<char>(is), {})) | boost::adaptors::transformed(tr_by(char2id)));
                // Do the rest
                do_rest_of_text_mode<std::uint16_t>(alphabet_size, id2char, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
            }
            else {
                using id_type = std::uint32_t;
                // map: char -> id
                map<char_type, id_type> char2id;
                for (size_t id = 0; id < alphabet_size; ++id) { char2id[id2char[id]] = id; }
                // input
                is.seekg(0);
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                const auto input = boost::copy_range<vector<id_type>>(
                        cvt.from_bytes(std::string(std::istreambuf_iterator<char>(is), {})) | boost::adaptors::transformed(tr_by(char2id)));
                // Do the rest
                do_rest_of_text_mode<std::uint32_t>(alphabet_size, id2char, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
            }
        }
        else if (p.get<string>("format") == "json") {
            JsonResultPrinter printer(std::cout, tr_by(id2char), cs, p.exist("show-all-positions"), p.exist("show-substring"));
            if (alphabet_size <= 0x100) {
                using id_type = std::uint8_t;
                // map: char -> id
                map<char_type, id_type> char2id;
                for (size_t id = 0; id < alphabet_size; ++id) { char2id[id2char[id]] = id; }
                // input
                is.seekg(0);
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                const auto input = boost::copy_range<vector<id_type>>(
                        cvt.from_bytes(std::string(std::istreambuf_iterator<char>(is), {})) | boost::adaptors::transformed(tr_by(char2id)));
                // Do the rest
                do_rest_of_text_mode<std::uint8_t>(alphabet_size, id2char, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
            }
            else if (alphabet_size <= 0x10000) {
                using id_type = std::uint16_t;
                // map: char -> id
                map<char_type, id_type> char2id;
                for (size_t id = 0; id < alphabet_size; ++id) { char2id[id2char[id]] = id; }
                // input
                is.seekg(0);
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                const auto input = boost::copy_range<vector<id_type>>(
                        cvt.from_bytes(std::string(std::istreambuf_iterator<char>(is), {})) | boost::adaptors::transformed(tr_by(char2id)));
                // Do the rest
                do_rest_of_text_mode<std::uint16_t>(alphabet_size, id2char, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
            }
            else {
                using id_type = std::uint32_t;
                // map: char -> id
                map<char_type, id_type> char2id;
                for (size_t id = 0; id < alphabet_size; ++id) { char2id[id2char[id]] = id; }
                // input
                is.seekg(0);
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                const auto input = boost::copy_range<vector<id_type>>(
                        cvt.from_bytes(std::string(std::istreambuf_iterator<char>(is), {})) | boost::adaptors::transformed(tr_by(char2id)));
                // Do the rest
                do_rest_of_text_mode<std::uint32_t>(alphabet_size, id2char, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
            }
        }
        else if (p.get<string>("format") == "json-lines") {
            JsonLinesResultPrinter printer(std::cout, tr_by(id2char), cs, p.exist("show-all-positions"), p.exist("show-substring"));
            if (alphabet_size <= 0x100) {
                using id_type = std::uint8_t;
                // map: char -> id
                map<char_type, id_type> char2id;
                for (size_t id = 0; id < alphabet_size; ++id) { char2id[id2char[id]] = id; }
                // input
                is.seekg(0);
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                const auto input = boost::copy_range<vector<id_type>>(
                        cvt.from_bytes(std::string(std::istreambuf_iterator<char>(is), {})) | boost::adaptors::transformed(tr_by(char2id)));
                // Do the rest
                do_rest_of_text_mode<std::uint8_t>(alphabet_size, id2char, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
            }
            else if (alphabet_size <= 0x10000) {
                using id_type = std::uint16_t;
                // map: char -> id
                map<char_type, id_type> char2id;
                for (size_t id = 0; id < alphabet_size; ++id) { char2id[id2char[id]] = id; }
                // input
                is.seekg(0);
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                const auto input = boost::copy_range<vector<id_type>>(
                        cvt.from_bytes(std::string(std::istreambuf_iterator<char>(is), {})) | boost::adaptors::transformed(tr_by(char2id)));
                // Do the rest
                do_rest_of_text_mode<std::uint16_t>(alphabet_size, id2char, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
            }
            else {
                using id_type = std::uint32_t;
                // map: char -> id
                map<char_type, id_type> char2id;
                for (size_t id = 0; id < alphabet_size; ++id) { char2id[id2char[id]] = id; }
                // input
                is.seekg(0);
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                const auto input = boost::copy_range<vector<id_type>>(
                        cvt.from_bytes(std::string(std::istreambuf_iterator<char>(is), {})) | boost::adaptors::transformed(tr_by(char2id)));
                // Do the rest
                do_rest_of_text_mode<std::uint32_t>(alphabet_size, id2char, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
            }
        }
        else if (p.get<string>("format") == "benchmark") {
            BenchmarkPrinter printer(std::cout, tr_by(id2char), cs, p.exist("show-all-positions"), p.exist("show-substring"));
            if (alphabet_size <= 0x100) {
                using id_type = std::uint8_t;
                // map: char -> id
                map<char_type, id_type> char2id;
                for (size_t id = 0; id < alphabet_size; ++id) { char2id[id2char[id]] = id; }
                // input
                is.seekg(0);
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                const auto input = boost::copy_range<vector<id_type>>(
                        cvt.from_bytes(std::string(std::istreambuf_iterator<char>(is), {})) | boost::adaptors::transformed(tr_by(char2id)));
                // Do the rest
                do_rest_of_text_mode<std::uint8_t>(alphabet_size, id2char, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
            }
            else if (alphabet_size <= 0x10000) {
                using id_type = std::uint16_t;
                // map: char -> id
                map<char_type, id_type> char2id;
                for (size_t id = 0; id < alphabet_size; ++id) { char2id[id2char[id]] = id; }
                // input
                is.seekg(0);
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                const auto input = boost::copy_range<vector<id_type>>(
                        cvt.from_bytes(std::string(std::istreambuf_iterator<char>(is), {})) | boost::adaptors::transformed(tr_by(char2id)));
                // Do the rest
                do_rest_of_text_mode<std::uint16_t>(alphabet_size, id2char, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
            }
            else {
                using id_type = std::uint32_t;
                // map: char -> id
                map<char_type, id_type> char2id;
                for (size_t id = 0; id < alphabet_size; ++id) { char2id[id2char[id]] = id; }
                // input
                is.seekg(0);
                std::wstring_convert<std::codecvt_utf8<char32_t>, char32_t> cvt;
                const auto input = boost::copy_range<vector<id_type>>(
                        cvt.from_bytes(std::string(std::istreambuf_iterator<char>(is), {})) | boost::adaptors::transformed(tr_by(char2id)));
                // Do the rest
                do_rest_of_text_mode<std::uint32_t>(alphabet_size, id2char, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
            }
        }
        else {
            throw runtime_error("Unsupported output format is specified.");
        }
    }
    else if (p.get<string>("mode") == "json-number-array") {
        using char_type = std::uint32_t;

        // Load
        is.seekg(0);
        json input_json;
        is >> input_json;

        // Type check
        bool format_is_ok = true;
        if (!input_json.is_array()) {
            format_is_ok = false;
        }
        else {
            for (auto& el : input_json.items()) {
                if (!el.value().is_number()) {
                    format_is_ok = false;
                    break;
                }
            }
        }
        if (!format_is_ok) {
            throw runtime_error("Input JSON must represent an array of numbers");
        }

        // alphabet size
        size_t alphabet_size = 0;
        for (auto& el : input_json.items()) {
            if (el.value().get<char_type>() >= alphabet_size) {
                alphabet_size = el.value().get<char_type>() + 1;
            }
        }

        NumberArrayJsonLinesResultPrinter printer(std::cout, cs, p.exist("show-all-positions"), p.exist("show-substring"));
        if (alphabet_size <= 0x100) {
            // Read into vector
            const auto input = boost::copy_range<vector<std::uint8_t>>(input_json.get<vector<std::uint8_t>>());
            // Do the rest
            do_rest_of_json_number_array_mode<std::uint8_t>(alphabet_size, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
        }
        else if (alphabet_size <= 0x10000) {
            // Read into vector
            const auto input = boost::copy_range<vector<std::uint16_t>>(input_json.get<vector<std::uint16_t>>());
            // Do the rest
            do_rest_of_json_number_array_mode<std::uint16_t>(alphabet_size, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
        }
        else {
            // Read into vector
            const auto input = boost::copy_range<vector<std::uint32_t>>(input_json.get<vector<std::uint32_t>>());
            // Do the rest
            do_rest_of_json_number_array_mode<std::uint32_t>(alphabet_size, std::move(input), printer, enum_type, resolution, block, window, range, constraint);
        }
    }
}

template <class ResultPrinter>
void do_rest_of_binary_mode(const std::size_t& alphabet_size, const std::vector<std::uint8_t> input, ResultPrinter& printer, Enumeration enum_type, const int resolution, const int block, const int window, std::pair<int, int> range, const substring_constraint& constraint)
{
    using namespace std;

    using id_type = std::uint8_t;

    // sast
    sast::sast<decltype(input), index_type> sast(input, alphabet_size);

    // enumerate substrings
    printer.print_header();
    switch (enum_type) {
    case Enumeration::BlumerStrings: {
        using substr_type = typename atm::blumer_substrings<decltype(input), index_type>::substr;

        atm::blumer_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::PurityMaximalStrings: {
        using substr_type = typename atm::purity_maximal_substrings<decltype(input), index_type>::substr;

        atm::purity_maximal_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::BranchingStrings: {
        using substr_type = typename atm::branching_substrings<decltype(input), index_type>::substr;

        atm::branching_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::NonLeafStrings: {
        using substr_type = typename atm::substrings<decltype(input), index_type>::substr;

        atm::substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::AllSubStrings: {
        using substr_type = typename atm::substrings_from_longest<decltype(input), index_type>::substr;

        atm::substrings_from_longest<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::AllBlockwiseSubStrings: {
        using substr_type = typename atm::coarse_substrings<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }

        atm::coarse_substrings<decltype(input), index_type> substrs(sast, resolution);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::ChunkedSlidingWindow: {
        using substr_type = typename atm::chunked_sliding_window<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }

        atm::chunked_sliding_window<decltype(input), index_type> substrs(sast, block, window);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::SlidingWindow: {
        using substr_type = typename atm::ngrams<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(window) > input.size()) {
            throw runtime_error("Specified N for N-grams is out of the size of the input.");
        }

        atm::ngrams<decltype(input), index_type> substrs(sast, window);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::BlockwiseSlidingWindow: {
        using substr_type = typename atm::coarse_ngrams<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }
        if (static_cast<std::size_t>(window) > input.size()) {
            throw runtime_error("Specified N for N-grams is out of the size of the input.");
        }

        atm::coarse_ngrams<decltype(input), index_type> substrs(sast, resolution, window);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::Words: {
        using substr_type = typename atm::words<decltype(input), index_type>::substr;

        atm::words<decltype(input), index_type> substrs(sast, [](const id_type id) { return std::isalpha(id); });
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::SingleSubstring: {
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
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    }
    printer.print_footer();
}

template<class ID, class ResultPrinter>
void do_rest_of_text_mode(const std::size_t& alphabet_size, const std::vector<std::uint32_t>& id2char, const std::vector<ID> input, ResultPrinter& printer, Enumeration enum_type, const int resolution, const int block, const int window, std::pair<int, int> range, const substring_constraint& constraint)
{
    using namespace std;

    using char_type = std::uint32_t;
    using id_type = ID;

    // sast
    sast::sast<decltype(input), index_type> sast(input, alphabet_size);

    // enumerate substrings
    printer.print_header();
    switch (enum_type) {
    case Enumeration::BlumerStrings: {
        using substr_type = typename atm::blumer_substrings<decltype(input), index_type>::substr;

        atm::blumer_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::PurityMaximalStrings: {
        using substr_type = typename atm::purity_maximal_substrings<decltype(input), index_type>::substr;

        atm::purity_maximal_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::BranchingStrings: {
        using substr_type = typename atm::branching_substrings<decltype(input), index_type>::substr;

        atm::branching_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::NonLeafStrings: {
        using substr_type = typename atm::substrings<decltype(input), index_type>::substr;

        atm::substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::AllSubStrings: {
        using substr_type = typename atm::substrings_from_longest<decltype(input), index_type>::substr;

        atm::substrings_from_longest<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::AllBlockwiseSubStrings: {
        using substr_type = typename atm::coarse_substrings<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }

        atm::coarse_substrings<decltype(input), index_type> substrs(sast, resolution);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::ChunkedSlidingWindow: {
        using substr_type = typename atm::chunked_sliding_window<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }

        atm::chunked_sliding_window<decltype(input), index_type> substrs(sast, block, window);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::SlidingWindow: {
        using substr_type = typename atm::ngrams<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(window) > input.size()) {
            throw runtime_error("Specified N for N-grams is out of the size of the input.");
        }

        atm::ngrams<decltype(input), index_type> substrs(sast, window);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::BlockwiseSlidingWindow: {
        using substr_type = typename atm::coarse_ngrams<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }
        if (static_cast<std::size_t>(window) > input.size()) {
            throw runtime_error("Specified N for N-grams is out of the size of the input.");
        }

        atm::coarse_ngrams<decltype(input), index_type> substrs(sast, resolution, window);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::Words: {
        using substr_type = typename atm::words<decltype(input), index_type>::substr;

        atm::words<decltype(input), index_type> substrs(sast, [&](const id_type id) { return std::isalpha(id2char[id]); });
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::SingleSubstring: {
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
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    }
    printer.print_footer();
}

template<class ID, class ResultPrinter>
void do_rest_of_json_number_array_mode(const std::size_t& alphabet_size, const std::vector<ID> input, ResultPrinter& printer, Enumeration enum_type, const int resolution, const int block, const int window, std::pair<int, int> range, const substring_constraint& constraint)
{
    using namespace std;

    using id_type = ID;

    // sast
    sast::sast<decltype(input), index_type> sast(input, alphabet_size);

    // enumerate substrings
    printer.print_header();
    switch (enum_type) {
    case Enumeration::BlumerStrings: {
        using substr_type = typename atm::blumer_substrings<decltype(input), index_type>::substr;

        atm::blumer_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::PurityMaximalStrings: {
        using substr_type = typename atm::purity_maximal_substrings<decltype(input), index_type>::substr;

        atm::purity_maximal_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::BranchingStrings: {
        using substr_type = typename atm::branching_substrings<decltype(input), index_type>::substr;

        atm::branching_substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::NonLeafStrings: {
        using substr_type = typename atm::substrings<decltype(input), index_type>::substr;

        atm::substrings<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::AllSubStrings: {
        using substr_type = typename atm::substrings_from_longest<decltype(input), index_type>::substr;

        atm::substrings_from_longest<decltype(input), index_type> substrs(sast);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::AllBlockwiseSubStrings: {
        using substr_type = typename atm::coarse_substrings<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }

        atm::coarse_substrings<decltype(input), index_type> substrs(sast, resolution);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::ChunkedSlidingWindow: {
        using substr_type = typename atm::chunked_sliding_window<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }

        atm::chunked_sliding_window<decltype(input), index_type> substrs(sast, block, window);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::SlidingWindow: {
        using substr_type = typename atm::ngrams<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(window) > input.size()) {
            throw runtime_error("Specified N for N-grams is out of the size of the input.");
        }

        atm::ngrams<decltype(input), index_type> substrs(sast, window);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::BlockwiseSlidingWindow: {
        using substr_type = typename atm::coarse_ngrams<decltype(input), index_type>::substr;

        if (static_cast<std::size_t>(resolution) > input.size()) {
            throw runtime_error("Specified resolution is out of the size of the input.");
        }
        if (static_cast<std::size_t>(window) > input.size()) {
            throw runtime_error("Specified N for N-grams is out of the size of the input.");
        }

        atm::coarse_ngrams<decltype(input), index_type> substrs(sast, resolution, window);
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    case Enumeration::Words: {
        throw runtime_error("Word enumeration is not supported when JSON array is given");
        break;
    }
    case Enumeration::SingleSubstring: {
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
        for (auto substr : boost::adaptors::filter(substrs, satisfy<substr_type>(constraint))) {
            printer.print(substr);
        }
        break;
    }
    }
    printer.print_footer();
}
