#include "CExplicitXY.h"
#include "CMendelianAnalyzer.h"
#include "CMixedPloidyReplay.h"
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include <zlib.h>
#include <dirent.h>
#include <algorithm>
#include <array>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace mendelian {
namespace xy {
struct Interval { std::string chr; int begin; int end; std::string region; };
struct Options {
    std::string sex = "unknown", par, reference, output, prefix = "out", pedigree;
    std::array<std::string, 3> files, samples;
    bool masked = false;
};
std::string Chromosome(std::string chr) {
    if (chr.compare(0, 3, "chr") == 0) chr.erase(0, 3);
    return chr;
}
bool XY(const std::string& chr) { return Chromosome(chr) == "X" || Chromosome(chr) == "Y"; }
std::vector<std::string> Lines(const std::string& path) {
    gzFile file = gzopen(path.c_str(), "rb");
    if (!file) throw std::runtime_error("Cannot open resource: " + path);
    std::vector<std::string> lines;
    std::string line;
    char buffer[8192];
    while (gzgets(file, buffer, sizeof(buffer))) {
        line += buffer;
        if (!line.empty() && line.back() == '\n') { lines.push_back(line); line.clear(); }
    }
    int error = Z_OK;
    gzerror(file, &error);
    const int closed = gzclose(file);
    if (error != Z_OK && error != Z_STREAM_END) throw std::runtime_error("Invalid gzip resource: " + path);
    if (closed != Z_OK) throw std::runtime_error("Truncated resource: " + path);
    if (!line.empty()) lines.push_back(line);
    return lines;
}
int Coordinate(const std::string& text) {
    if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
        throw std::runtime_error("Invalid BED coordinate: " + text);
    const unsigned long long value = std::stoull(text);
    if (value > 2147483647ULL) throw std::runtime_error("BED coordinate exceeds supported reference length");
    return static_cast<int>(value);
}
std::vector<Interval> ReadIntervals(const std::string& file, const std::string& expected,
                                    const std::string& region) {
    std::vector<Interval> intervals;
    for (const auto& line : Lines(file)) {
        if (line.empty() || line[0] == '#' || line == "\n") continue;
        std::istringstream input(line);
        std::string chr, begin, end, extra;
        if (!(input >> chr >> begin >> end)) throw std::runtime_error("Malformed BED row: " + file);
        if (Chromosome(chr) != expected) throw std::runtime_error("Wrong chromosome in BED: " + file);
        Interval interval{expected, Coordinate(begin), Coordinate(end), region};
        if (interval.begin >= interval.end) throw std::runtime_error("Empty/reversed BED interval: " + file);
        intervals.push_back(interval);
    }
    if (intervals.empty()) throw std::runtime_error("Empty required BED: " + file);
    return intervals;
}
struct Policy {
    std::vector<Interval> intervals;
    std::vector<Interval> names;
    explicit Policy(const std::string& directory) {
        if (directory.empty()) return;
        DIR* dir = opendir(directory.c_str());
        if (!dir) throw std::runtime_error("Cannot open --par-dir");
        std::vector<std::string> files;
        while (dirent* entry = readdir(dir)) files.push_back(entry->d_name);
        closedir(dir);
        for (const std::string chr : {"X", "Y"}) for (const std::string type : {"PAR", "nonPAR"}) {
            const std::string suffix = "chr" + chr + "_" + type + ".bed";
            std::vector<std::string> matches;
            for (std::string file : files) {
                std::string plain = file;
                if (plain.size() >= 3 && plain.substr(plain.size() - 3) == ".gz") plain.resize(plain.size() - 3);
                if (plain.size() >= suffix.size() && plain.substr(plain.size() - suffix.size()) == suffix) matches.push_back(file);
            }
            if (matches.size() != 1) throw std::runtime_error("Expected exactly one *" + suffix + "[.gz]");
            const auto more = ReadIntervals(directory + "/" + matches[0], chr, chr + "_" + type);
            intervals.insert(intervals.end(), more.begin(), more.end());
        }
        std::sort(intervals.begin(), intervals.end(), [](const Interval& a, const Interval& b) {
            return a.chr < b.chr || (a.chr == b.chr && a.begin < b.begin);
        });
        for (size_t i = 1; i < intervals.size(); ++i)
            if (intervals[i].chr == intervals[i-1].chr && intervals[i].begin < intervals[i-1].end)
                throw std::runtime_error("Overlapping PAR/nonPAR BED intervals");
        // Optional assembly-specific naming only. No ordinal/fragment inference.
        if (std::find(files.begin(), files.end(), "par_regions.tsv") != files.end()) {
            for (const auto& line : Lines(directory + "/par_regions.tsv")) {
                if (line.empty() || line[0] == '#' || line == "\n") continue;
                std::istringstream input(line);
                std::string chr, begin, end, label;
                if (!(input >> chr >> begin >> end >> label)) throw std::runtime_error("Malformed par_regions.tsv");
                if (chr == "chrom" || chr == "chromosome" || chr == "contig") continue;
                chr = Chromosome(chr);
                if ((chr != "X" && chr != "Y") || (label != "PAR1" && label != "PAR2" && label != chr + "_PAR1" && label != chr + "_PAR2"))
                    throw std::runtime_error("Invalid PAR name in par_regions.tsv");
                if (label == "PAR1" || label == "PAR2") label = chr + "_" + label;
                Interval name{chr, Coordinate(begin), Coordinate(end), label};
                if (name.begin >= name.end) throw std::runtime_error("Invalid PAR naming interval");
                names.push_back(name);
            }
            for (size_t i = 0; i < names.size(); ++i) for (size_t j = i+1; j < names.size(); ++j)
                if (names[i].chr == names[j].chr && names[i].begin < names[j].end && names[j].begin < names[i].end)
                    throw std::runtime_error("Overlapping PAR naming intervals");

        }
    }
    int Locate(const std::string& chromosome, int begin, int end, std::string& region, std::string& reason) const {
        const std::string chr = Chromosome(chromosome);
        if (!XY(chr)) { region = "autosome"; return -2; }
        region = "unclassified";
        if (intervals.empty()) { reason = "missing_region_resources"; return -1; }
        for (size_t i = 0; i < intervals.size(); ++i) {
            const auto& interval = intervals[i];
            if (interval.chr == chr && interval.begin <= begin && interval.end >= end) {
                region = interval.region;
                if (region == chr + "_PAR") for (const auto& name : names)
                    if (name.chr == chr && name.begin <= begin && name.end >= end) region = name.region;
                return static_cast<int>(i);
            }
        }
        reason = "unclassified_region";
        for (const auto& interval : intervals)
            if (interval.chr == chr && begin < interval.end && end > interval.begin) reason = "boundary_crossing";
        return -1;
    }
};

struct Record {
    bcf1_t* raw = nullptr;
    int source = 0;
    int ordinal = 0;
    std::string chr;
    core::SMixedVariant variant;
    std::string region, reason, status = "not_evaluable";
    int domain = -1;
    int spanEnd = 0;
    bool complete = true;
    bool nonref = false;
    bool phased = false;
};
struct Input {
    bcf_hdr_t* header = nullptr;
    std::vector<Record> records;
    ~Input() { for (auto& record : records) bcf_destroy(record.raw); if (header) bcf_hdr_destroy(header); }
    Input() = default;
    Input(const Input&) = delete;
    Input& operator=(const Input&) = delete;
};
Options Parse(int argc, char** argv) {
    Options options;
    const std::map<std::string, int> roles{{"-mother", 0}, {"-father", 1}, {"-child", 2}};
    const std::map<std::string, int> samples{{"-sample-mother", 0}, {"-sample-father", 1}, {"-sample-child", 2}};
    for (int i = 2; i < argc; ++i) {
        const std::string key = argv[i];
        if (key == "--ypar-masked") { options.masked = true; continue; }
        if (key == "--help") {
            std::cout << "Explicit native XY mode: --child-sex male|female|unknown --par-dir DIR [--ypar-masked]\n"
                      << "Use -mother/-father/-child VCF, optional -sample-* NAME, -ref FASTA, -outDir DIR and -out-prefix NAME.\n"
                      << "--ypar-masked is a flag; =false is invalid. Original GTs are preserved.\n"
                      << "Outputs PREFIX_{mother,father,child}_xy.vcf and PREFIX_xy_summary.tsv.\n";
            return Options{};
        }
        if (key == "--autosome-only" || key == "--disable-ref-overlap") continue;
        if (i + 1 >= argc) throw std::runtime_error("Missing value for " + key);
        const std::string value = argv[++i];
        if (roles.count(key)) options.files[roles.at(key)] = value;
        else if (samples.count(key)) options.samples[samples.at(key)] = value;
        else if (key == "--child-sex") options.sex = value;
        else if (key == "--par-dir") options.par = value;
        else if (key == "-ref") options.reference = value;
        else if (key == "-outDir") options.output = value;
        else if (key == "-out-prefix") options.prefix = value;
        else if (key == "-pedigree") options.pedigree = value;
        else if (key == "-thread-count" || key == "-no-call" || key == "-output-info-tags") {}
        else if (key == "-filter" && (value == "PASS" || value == "none")) {} 
        else throw std::runtime_error("Unsupported option in explicit XY mode: " + key);
    }
    if (options.sex != "male" && options.sex != "female" && options.sex != "unknown") throw std::runtime_error("Invalid --child-sex");
    if (options.masked && options.par.empty()) throw std::runtime_error("--ypar-masked requires valid --par-dir");
    if (options.reference.empty() || options.output.empty() || options.files[0].empty() || options.files[1].empty() || options.files[2].empty())
        throw std::runtime_error("Explicit XY requires -mother, -father, -child, -ref and -outDir");
    if (!options.pedigree.empty()) {
        struct Person { std::string id, father, mother, sex; };
        std::vector<Person> children;
        for (const auto& line : Lines(options.pedigree)) {
            if (line.empty() || line[0] == '#' || line == "\n") continue;
            std::istringstream input(line);
            std::string family;
            Person person;
            if (!(input >> family >> person.id >> person.father >> person.mother >> person.sex)) throw std::runtime_error("Malformed pedigree");
            if (person.father != "0" && person.mother != "0") children.push_back(person);
        }
        if (children.size() != 1) throw std::runtime_error("Explicit XY requires one unambiguous child in pedigree");
        const auto& child = children[0];
        const std::array<std::string, 3> ids{{child.mother, child.father, child.id}};
        for (int role = 0; role < 3; ++role) {
            if (!options.samples[role].empty() && options.samples[role] != ids[role]) throw std::runtime_error("Sample and pedigree disagree");
            options.samples[role] = ids[role];
        }
        if ((child.sex == "1" && options.sex == "female") || (child.sex == "2" && options.sex == "male"))
            throw std::runtime_error("--child-sex conflicts with pedigree");
    }
    return options;
}
void Read(Input& input, const Options& options, int role) {
    htsFile* file = bcf_open(options.files[role].c_str(), "r");
    if (!file) throw std::runtime_error("Cannot open input VCF");
    input.header = bcf_hdr_read(file);
    if (!input.header) { bcf_close(file); throw std::runtime_error("Cannot read VCF header"); }
    int sample = -1;
    if (!options.samples[role].empty()) sample = bcf_hdr_id2int(input.header, BCF_DT_SAMPLE, options.samples[role].c_str());
    else if (bcf_hdr_nsamples(input.header) == 1) sample = 0;
    if (sample < 0) { bcf_close(file); throw std::runtime_error("Select a unique sample with -sample-* or pedigree"); }
    bcf1_t* raw = bcf_init();
    int32_t* gt = nullptr;
    int capacity = 0;
    int32_t* ps = nullptr;
    int psCapacity = 0;
    int result;
    int ordinal = 0;
    while ((result = bcf_read(file, input.header, raw)) == 0) {
        ++ordinal;
        bcf_unpack(raw, BCF_UN_ALL);
        if (!XY(bcf_hdr_id2name(input.header, raw->rid))) continue;
        Record record;
        record.raw = bcf_dup(raw);
        record.source = role;
        record.ordinal = ordinal;
        record.chr = bcf_hdr_id2name(input.header, raw->rid);
        record.variant.sample = role;
        record.variant.position = raw->pos;
        record.spanEnd = raw->pos + std::max<int>(1, raw->rlen);
        for (int a = 0; a < raw->n_allele; ++a) record.variant.alleles.push_back(raw->d.allele[a]);
        const int count = bcf_get_genotypes(input.header, raw, &gt, &capacity);
        if (count <= 0) record.complete = false;
        else {
            const int stride = count / bcf_hdr_nsamples(input.header);
            for (int a = 0; a < stride; ++a) {
                const int32_t value = gt[sample * stride + a];
                if (value == bcf_int32_vector_end) break;
                if (bcf_gt_is_missing(value)) record.complete = false;
                const int allele = bcf_gt_is_missing(value) ? -1 : bcf_gt_allele(value);
                if (allele < 0 || allele >= raw->n_allele) record.complete = false;
                record.variant.genotype.push_back(allele);
                record.nonref = record.nonref || allele > 0;
                if (a > 0 && bcf_gt_is_phased(value)) record.phased = true;
            }
        }
        if (record.variant.genotype.empty()) record.complete = false;
        if (record.phased) {
            const int countPS = bcf_get_format_int32(input.header, raw, "PS", &ps, &psCapacity);
            record.variant.phaseSet = countPS > sample && ps[sample] != bcf_int32_missing
                                    ? std::to_string(ps[sample]) : "implicit";
        }
        input.records.push_back(std::move(record));
    }
    free(gt); free(ps); bcf_destroy(raw); bcf_close(file);
    if (result < -1) throw std::runtime_error("Malformed/truncated VCF");
}
core::SMixedModel Model(const Record& record, const Options& options) {
    if (record.region == "autosome") return {2, 2, 2};
    if (record.region == "X_nonPAR") return options.sex == "male" ? core::SMixedModel{2, 0, 1} : core::SMixedModel{2, 1, 2};
    if (record.region == "Y_nonPAR") return {0, 1, 1};
    return {2, 2, 2};
}
void Eligibility(Record& record, const Options& options, const Policy& policy) {
    int end = record.spanEnd;
    bool insertion = record.variant.alleles.size() > 1;
    for (size_t i = 1; i < record.variant.alleles.size(); ++i) {
        const auto& alt = record.variant.alleles[i];
        const auto& ref = record.variant.alleles[0];
        insertion = insertion && alt.size() > ref.size() && alt.compare(0, ref.size(), ref) == 0;
    }
    if (insertion) end = record.variant.position + 1;
    record.domain = policy.Locate(record.chr, record.variant.position, end, record.region, record.reason);
    const std::string chr = Chromosome(record.chr);
    if (XY(chr) && options.sex == "unknown") { record.reason = "unknown_child_sex"; return; }
    if (!record.reason.empty()) return;
    if (record.region.compare(0, 5, "Y_PAR") == 0) record.reason = "y_par_excluded";
    else if (record.region.compare(0, 5, "X_PAR") == 0 && !options.masked) record.reason = "par_requires_ypar_masked";
    else if (record.region == "Y_nonPAR" && options.sex == "female") record.reason = "female_y_nonpar";
    if (!record.reason.empty()) return;
    const core::SMixedModel model = Model(record, options);
    const int expected = record.source == 0 ? model.motherPloidy : record.source == 1 ? model.fatherPloidy : model.childPloidy;
    // Non-participating calls are preserved but never constrain eligibility.
    if (!expected) return;
    if (!record.complete) { record.reason = "required_gt_missing"; return; }
    auto& genotype = record.variant.genotype;
    if (expected == 1 && genotype.size() == 2) {
        if (genotype[0] != genotype[1]) { record.reason = "haploid_heterozygous"; return; }
        genotype.resize(1); // computation only; raw output remains unchanged
    }
    if (static_cast<int>(genotype.size()) != expected) { record.reason = "unexpected_ploidy"; return; }
    for (const auto& allele : record.variant.alleles)
        if (allele.empty() || allele.find_first_not_of("ACGTNacgtn") != std::string::npos) {
            record.reason = "unsupported_variant"; return;
        }
}
}
bool HasExplicitXY(int argc, char** argv) {
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg.compare(0, 11, "--child-sex") == 0 || arg.compare(0, 9, "--par-dir") == 0 || arg.compare(0, 13, "--ypar-masked") == 0) return true;
    }
    return false;
}
namespace xy {
std::string Status(const std::string& reason) {
    return reason == "female_y_nonpar" || reason == "unsupported_variant" || reason == "unsupported_replay"
         ? "unsupported" : "not_evaluable";
}
struct Domain {
    std::string chr;
    int interval;
    std::vector<Record*> records;
};
void Analyze(std::array<Input, 3>& inputs, const Options& options, const Policy& policy) {
    std::map<std::pair<std::string, int>, Domain> domains;
    int id = 0;
    for (auto& input : inputs) for (auto& record : input.records) {
        record.variant.id = id++;
        Eligibility(record, options, policy);
    }
    // Missing required calls outrank ploidy errors at the same observation.
    for (auto& input : inputs) for (auto& record : input.records) {
        if (record.reason.empty() || record.reason == "haploid_heterozygous" || record.reason == "unexpected_ploidy") {
            const core::SMixedModel model = Model(record, options);
            const std::array<int, 3> required{{model.motherPloidy, model.fatherPloidy, model.childPloidy}};
            for (const auto& otherInput : inputs) for (const auto& other : otherInput.records)
                if (required[other.source] && !other.complete && other.reason == "required_gt_missing"
                    && Chromosome(other.chr) == Chromosome(record.chr) && other.domain == record.domain
                    && record.variant.position < other.spanEnd && other.variant.position < record.spanEnd)
                    record.reason = "required_gt_missing";
        }
        if (!record.reason.empty()) { record.status = Status(record.reason); continue; }
        auto& domain = domains[{Chromosome(record.chr), record.domain}];
        domain.chr = record.chr;
        domain.interval = record.domain;
        domain.records.push_back(&record);
    }
    for (auto& item : domains) {
        Domain& domain = item.second;
        std::stable_sort(domain.records.begin(), domain.records.end(), [](const Record* a, const Record* b) {
            return a->variant.position < b->variant.position;
        });
        const core::SMixedModel model = Model(*domain.records.front(), options);
        const std::array<int, 3> ploidies{{model.motherPloidy, model.fatherPloidy, model.childPloidy}};
        // Invalid REQUIRED calls are replay barriers; irrelevant parents cannot
        // create gaps or remove a chromosome from the required-member union.
        std::vector<std::pair<int, int>> barriers;
        for (const auto& input : inputs) for (const auto& record : input.records)
            if (Chromosome(record.chr) == Chromosome(domain.chr) && !record.reason.empty() && ploidies[record.source]
                && (record.reason == "required_gt_missing" || record.reason == "unexpected_ploidy"
                    || record.reason == "haploid_heterozygous" || record.reason == "unsupported_variant"
                    || record.reason == "boundary_crossing"))
                barriers.emplace_back(record.variant.position, record.spanEnd);
        std::sort(barriers.begin(), barriers.end());
        std::vector<std::pair<int, int>> merged;
        for (const auto& barrier : barriers) {
            if (!merged.empty() && barrier.first <= merged.back().second) merged.back().second = std::max(merged.back().second, barrier.second);
            else merged.push_back(barrier);
        }
        int domainBegin = domain.interval >= 0 ? policy.intervals[domain.interval].begin : 0;
        faidx_t* fasta = fai_load(options.reference.c_str());
        if (!fasta) throw std::runtime_error("Cannot load FASTA/index");
        int length = 0;
        char* sequence = faidx_fetch_seq(fasta, domain.chr.c_str(), 0, 2147483646, &length);
        if (!sequence || length < 0) {
            free(sequence);
            const std::string alternate = domain.chr.compare(0, 3, "chr") == 0 ? domain.chr.substr(3) : "chr" + domain.chr;
            sequence = faidx_fetch_seq(fasta, alternate.c_str(), 0, 2147483646, &length);
        }
        fai_destroy(fasta);
        if (!sequence || length < 0) { free(sequence); throw std::runtime_error("Missing FASTA contig: " + domain.chr); }
        std::string reference(sequence, length);
        free(sequence);
        std::transform(reference.begin(), reference.end(), reference.begin(), [](unsigned char c) { return std::toupper(c); });
        const int domainEnd = domain.interval >= 0 ? policy.intervals[domain.interval].end : length;
        if (domainEnd > length) throw std::runtime_error("BED exceeds FASTA contig length");
        std::vector<std::pair<int, int>> segments;
        for (const auto& barrier : merged) {
            if (barrier.second <= domainBegin || barrier.first >= domainEnd) continue;
            if (barrier.first > domainBegin) segments.emplace_back(domainBegin, barrier.first);
            domainBegin = std::max(domainBegin, barrier.second);
        }
        if (domainBegin < domainEnd) segments.emplace_back(domainBegin, domainEnd);
        std::set<int> processed;
        for (const auto& segment : segments) {
            std::vector<Record*> records;
            std::array<bool, 3> seen{{false, false, false}};
            bool nonref = false;
            for (Record* record : domain.records)
                if (record->variant.position >= segment.first && record->spanEnd <= segment.second) {
                    records.push_back(record);
                    processed.insert(record->variant.id);
                    if (ploidies[record->source]) { seen[record->source] = true; nonref = nonref || record->nonref; }
                }
            if (records.empty()) continue;
            bool complete = true;
            // A variant-only VCF replays reference between records, exactly as
            // legacy VBT; this is NOT evidence of a complete explicit GT.
            // Explicit missing records were retained above as barriers.
            (void)seen;
            if (!complete || !nonref) {
                for (Record* record : records) {
                    record->reason = complete ? "all_required_reference" : "required_gt_missing";
                    record->status = "not_evaluable";
                }
                continue;
            }
            std::vector<core::SMixedVariant> variants;
            std::map<int, Record*> byID;
            for (Record* record : records) {
                byID[record->variant.id] = record;
                if (!ploidies[record->source]) continue;
                core::SMixedVariant variant = record->variant;
                for (auto& allele : variant.alleles)
                    std::transform(allele.begin(), allele.end(), allele.begin(), [](unsigned char c) { return std::toupper(c); });
                variants.push_back(std::move(variant));
            }
            const core::SMixedResult result = core::ReplayMixed(reference, segment.first, segment.second, model, variants);
            for (Record* record : records) { record->status = "consistent"; record->reason = "pass"; }
            for (int excluded : result.excluded) {
                byID.at(excluded)->status = "violation";
                byID.at(excluded)->reason = "mendelian_violation";
            }
            // Reference-only records remain visible but are not denominator
            // observations. Irrelevant parental ALT cannot open that gate.
            for (Record* record : records) {
                bool eventNonref = false;
                for (const Record* other : records)
                    if (ploidies[other->source] && other->nonref
                        && record->variant.position < other->spanEnd
                        && other->variant.position < record->spanEnd) eventNonref = true;
                if (!eventNonref) { record->status = "not_evaluable"; record->reason = "all_required_reference"; }
            }
            // Attribution is representation dependent, not de-novo localization.
            for (int excluded : result.excluded) {
                const Record* event = byID.at(excluded);
                for (Record* record : records)
                    if (record->variant.position < event->spanEnd && event->variant.position < record->spanEnd) {
                        record->status = "violation"; record->reason = "mendelian_violation";
                    }
            }
        }
        for (Record* record : domain.records) if (!processed.count(record->variant.id)) {
            record->reason = "required_gt_missing"; record->status = "not_evaluable";
            for (const auto& input : inputs) for (const auto& other : input.records)
                if (Chromosome(other.chr) == Chromosome(record->chr) && ploidies[other.source] && !other.reason.empty()
                    && record->variant.position < other.spanEnd && other.variant.position < record->spanEnd) {
                    record->reason = other.reason; record->status = Status(other.reason); break;
                }
        }
    }
}
void Write(std::array<Input, 3>& inputs, const Options& options) {
    std::map<std::string, unsigned long> counts;
    const std::array<std::string, 3> roles{{"mother", "father", "child"}};
    for (int role = 0; role < 3; ++role) {
        Input& input = inputs[role];
        for (const std::string tag : {"XY_STATUS", "XY_REGION", "XY_REASON", "XY_SOURCE_ID"}) {
            const int id = bcf_hdr_id2int(input.header, BCF_DT_ID, tag.c_str());
            if (id >= 0 && bcf_hdr_idinfo_exists(input.header, BCF_HL_INFO, id)) {
                if (bcf_hdr_id2type(input.header, BCF_HL_INFO, id) != BCF_HT_STR
                    || bcf_hdr_id2length(input.header, BCF_HL_INFO, id) != BCF_VL_FIXED
                    || bcf_hdr_id2number(input.header, BCF_HL_INFO, id) != 1)
                    throw std::runtime_error("Incompatible existing INFO/" + tag);
            } else {
                const std::string line = "##INFO=<ID=" + tag + ",Number=1,Type=String,Description=\"Native explicit XY " + tag + "\">";
                if (bcf_hdr_append(input.header, line.c_str()) < 0) throw std::runtime_error("Cannot append XY header");
            }
        }
        bcf_hdr_append(input.header, "##XY_POLICY=native-joint-replay-v1");
        bcf_hdr_append(input.header, ("##XY_CHILD_SEX=" + options.sex).c_str());
        bcf_hdr_append(input.header, options.masked ? "##XY_YPAR_MASKED_DECLARATION=true" : "##XY_YPAR_MASKED_DECLARATION=false");
        if (bcf_hdr_sync(input.header) < 0) throw std::runtime_error("Cannot synchronize output header");
        const std::string path = options.output + "/" + options.prefix + "_" + roles[role] + "_xy.vcf";
        htsFile* file = bcf_open(path.c_str(), "w");
        if (!file) throw std::runtime_error("Cannot write " + path);
        if (bcf_hdr_write(file, input.header) < 0) { bcf_close(file); throw std::runtime_error("Cannot write VCF header"); }
        for (auto& record : input.records) {
            const std::string sourceID = roles[role] + ":" + std::to_string(record.ordinal);
            if (bcf_update_info_string(input.header, record.raw, "XY_SOURCE_ID", sourceID.c_str()) < 0) throw std::runtime_error("Cannot write source record ID");
            if (record.reason.empty()) throw std::runtime_error("Unclassified output record");
            if (bcf_update_info_string(input.header, record.raw, "XY_STATUS", record.status.c_str()) < 0
                || bcf_update_info_string(input.header, record.raw, "XY_REGION", record.region.c_str()) < 0
                || bcf_update_info_string(input.header, record.raw, "XY_REASON", record.reason.c_str()) < 0
                || bcf_write(file, input.header, record.raw) < 0) {
                bcf_close(file); throw std::runtime_error("Cannot write annotated record");
            }
            ++counts[roles[role] + "\t" + record.region + "\t" + record.status + "\t" + record.reason];
        }
        if (bcf_close(file) < 0) throw std::runtime_error("Output VCF close failed");
    }
    std::ofstream summary(options.output + "/" + options.prefix + "_xy_summary.tsv");
    if (!summary) throw std::runtime_error("Cannot write XY summary");
    summary << "source\tregion\tstatus\treason\trecords\n";
    for (const auto& count : counts) summary << count.first << '\t' << count.second << '\n';
    summary.flush();
    if (!summary) throw std::runtime_error("Cannot finish XY summary");
}
}
int RunExplicitXY(int argc, char** argv) {
    try {
        const xy::Options options = xy::Parse(argc, argv);
        if (options.reference.empty()) return 0;
        const xy::Policy policy(options.par);
        std::array<xy::Input, 3> inputs;
        for (int role = 0; role < 3; ++role) xy::Read(inputs[role], options, role);
        xy::Analyze(inputs, options, policy);
        xy::Write(inputs, options);
        // Remove only the new policy arguments. The original autosome analyzer,
        // filters, matching, decisions, writer and counters remain unchanged.
        std::vector<std::string> legacy;
        for (int i = 0; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--child-sex" || arg == "--par-dir") { ++i; continue; }
            if (arg == "--ypar-masked") continue;
            legacy.push_back(arg);
        }
        legacy.push_back("--autosome-only");
        std::vector<char*> args;
        for (auto& arg : legacy) args.push_back(&arg[0]);
        args.push_back(nullptr);
        CMendelianAnalyzer autosomes;
        return autosomes.run(static_cast<int>(legacy.size()), args.data());
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Explicit XY error: " << error.what() << '\n';
        return 1;
    }
}
}
