#include "CMixedPloidyReplay.h"
#include "CHaplotypeSequence.h"
#include "CBaseVariantProvider.h"
#include <algorithm>
#include <array>
#include <map>
#include <sstream>
#include <stdexcept>
#include <utility>
namespace core {
namespace {
struct Track { int cursor; std::string pending; };
struct State {
    std::array<Track, 4> tracks;
    std::map<std::string, int> phase;
    std::vector<int> included;
    std::vector<int> excluded;
};
void Advance(Track& track, int position, const std::string& reference) {
    if (position > track.cursor) {
        track.pending.append(reference, track.cursor, position - track.cursor);
        track.cursor = position;
    }
}
bool Match(Track& first, Track& second) {
    const size_t length = std::min(first.pending.size(), second.pending.size());
    if (first.pending.compare(0, length, second.pending, 0, length) != 0)
        return false;
    first.pending.erase(0, length);
    second.pending.erase(0, length);
    return true;
}
bool Prune(State& state, int position, const std::string& reference, int pairs) {
    for (int i = 0; i < pairs; ++i) {
        Track& a = state.tracks[2 * i];
        Track& b = state.tracks[2 * i + 1];
        if (a.pending.empty() && b.pending.empty() && a.cursor == b.cursor && a.cursor < position) {
            a.cursor = position; b.cursor = position;
        } else { Advance(a, position, reference); Advance(b, position, reference); }
    }
    for (int i = 0; i < pairs; ++i)
        if (!Match(state.tracks[2 * i], state.tracks[2 * i + 1]))
            return false;
    return true;
}
// The same provider context trimmer is used by the legacy diploid path.
// Its implementation is stateless; the adapter adds no normalization rules.
struct ContextTrimmer : CBaseVariantProvider {
    using CBaseVariantProvider::FindOptimalTrimmings;
};
std::map<int, CVariant> Normalize(const std::vector<SMixedVariant>& variants) {
    std::map<int, CVariant> result;
    for (int sample = 0; sample < 3; ++sample) {
        std::vector<CVariant> ambiguous;
        std::vector<std::vector<CVariant>> ordinary(1);
        for (const auto& variant : variants) {
            if (variant.sample != sample) continue;
            CVariant native;
            native.m_nId = variant.id;
            native.m_nVcfId = sample;
            native.m_nChrId = 0;
            native.m_nOriginalPos = variant.position;
            native.m_nStartPos = variant.position;
            native.m_nEndPos = variant.position + static_cast<int>(variant.alleles.at(0).size());
            native.m_refSequence = variant.alleles.at(0);
            native.m_nZygotCount = static_cast<int>(variant.genotype.size());
            native.m_bIsHeterozygous = native.m_nZygotCount == 2 && variant.genotype[0] != variant.genotype[1];
            native.m_nAlleleCount = native.m_bIsHeterozygous ? 2 : 1;
            for (int slot = 0; slot < native.m_nZygotCount; ++slot) {
                const int allele = variant.genotype.at(slot);
                native.m_genotype[slot] = allele;
                SAllele& selected = native.m_alleles[slot];
                selected.m_sequence = variant.alleles.at(allele);
                selected.m_nStartPos = native.m_nStartPos;
                selected.m_nEndPos = native.m_nEndPos;
                selected.m_bIsIgnored = allele == 0;
                if (selected.m_bIsIgnored) continue;
                unsigned int first = 0, last = 0;
                native.GetMaxTrimStartEnd(slot, first, last);
                const bool multiple = selected.m_sequence.size() >= 2
                    && native.m_refSequence.size() >= selected.m_sequence.size()
                    && first + last > std::min(native.m_refSequence.size(), selected.m_sequence.size());
                if (multiple) native.m_bHaveMultipleTrimOption = true;
                else native.TrimVariant(slot, false); // legacy mendelian default
            }
            if (native.m_bHaveMultipleTrimOption) ambiguous.push_back(native);
            else {
                int first = 2147483647, last = -1;
                for (int slot = 0; slot < native.m_nAlleleCount; ++slot)
                    if (!native.m_alleles[slot].m_bIsIgnored) {
                        first = std::min(first, native.m_alleles[slot].m_nStartPos);
                        last = std::max(last, native.m_alleles[slot].m_nEndPos);
                    }
                native.m_nStartPos = first == 2147483647 ? variant.position : first;
                native.m_nEndPos = last == -1 ? variant.position : last;
                ordinary[0].push_back(native);
            }
        }
        // CVariant has a legacy partial copy constructor. Sort indices rather
        // than invoking its implicit assignment on partially initialized fields.
        std::vector<size_t> order;
        for (size_t i = 0; i < ordinary[0].size(); ++i) order.push_back(i);
        std::stable_sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            return ordinary[0][a].GetStart() < ordinary[0][b].GetStart();
        });
        std::vector<CVariant> sorted;
        for (size_t i : order) sorted.push_back(ordinary[0][i]);
        ordinary[0].swap(sorted);
        SConfig config;
        config.m_bTrimBeginningFirst = false;
        ContextTrimmer::FindOptimalTrimmings(ambiguous, &ordinary, config);
        for (const auto& native : ordinary[0]) result.emplace(native.m_nId, native);
        for (const auto& native : ambiguous) result.emplace(native.m_nId, native);
    }
    return result;
}
// Replay one selected actual chromosome through VBT single-haplotype playback.
bool Apply(Track& track, const CVariant& native, int slot,
           const std::string& reference) {
    const SAllele& normalized = native.m_alleles[slot];
    if (normalized.m_bIsIgnored) return true;
    const int start = normalized.m_nStartPos;
    const int stop = normalized.m_nEndPos;
    if (start == stop && normalized.m_sequence.empty()) return true;
    if (start < track.cursor) return false;
    Advance(track, start, reference);
    CHaplotypeSequence sequence(reference.data(), static_cast<int>(reference.size()));
    sequence.AddVariant(COrientedVariant(native, slot));
    sequence.MoveForward(start);
    while (!sequence.IsOnTemplate() || sequence.GetTemplatePosition() < stop) {
        track.pending += sequence.NextBase();
        sequence.Next();
    }
    track.cursor = stop;
    return true;
}
std::string Key(const State& state, int pairs) {
    std::ostringstream out;
    for (int i = 0; i < 2 * pairs; ++i)
        out << state.tracks[i].cursor << ':' << state.tracks[i].pending.size()
            << ':' << state.tracks[i].pending << ';';
    for (const auto& phase : state.phase)
        out << phase.first.size() << ':' << phase.first << ':' << phase.second << ';';
    return out.str();
}
void Keep(std::map<std::string, State>& states, State state, int pairs) {
    const std::string key = Key(state, pairs);
    auto found = states.find(key);
    if (found == states.end() || state.included.size() > found->second.included.size())
        states[key] = std::move(state);
}
}
SMixedResult ReplayMixed(const std::string& reference, int begin, int end,
                        const SMixedModel& model, std::vector<SMixedVariant> variants) {
    if (begin < 0 || end < begin || end > static_cast<int>(reference.size()))
        throw std::runtime_error("Invalid mixed-ploidy replay domain");
    const bool maternal = model.motherPloidy != 0;
    const bool paternal = model.fatherPloidy != 0;
    const int pairs = maternal && paternal ? 2 : 1;
    if ((!maternal && !paternal) || model.childPloidy != pairs)
        throw std::runtime_error("Invalid mixed-ploidy transmission model");
    std::stable_sort(variants.begin(), variants.end(), [](const SMixedVariant& a, const SMixedVariant& b) {
        return a.position < b.position;
    });
    // Validate before the legacy two-slot representation is constructed.
    for (const auto& variant : variants) {
        if (variant.genotype.empty() || variant.genotype.size() > 2) throw std::runtime_error("Invalid replay ploidy");
        for (int allele : variant.genotype)
            if (allele < 0 || allele >= static_cast<int>(variant.alleles.size())) throw std::runtime_error("Invalid replay allele");
    }
    const auto normalized = Normalize(variants);
    std::stable_sort(variants.begin(), variants.end(), [&](const SMixedVariant& a, const SMixedVariant& b) {
        return normalized.at(a.id).GetStart() < normalized.at(b.id).GetStart();
    });
    State initial;
    for (auto& track : initial.tracks) track.cursor = begin;
    std::map<std::string, State> states;
    Keep(states, initial, pairs);
    unsigned long visited = 0;
    for (const auto& variant : variants) {
        const int expected = variant.sample == 0 ? model.motherPloidy
                           : variant.sample == 1 ? model.fatherPloidy : model.childPloidy;
        if (!expected) continue;
        if (variant.position < begin || variant.position + static_cast<int>(variant.alleles.at(0).size()) > end
            || static_cast<int>(variant.genotype.size()) != expected)
            throw std::runtime_error("Variant does not belong to its replay domain/model");
        for (int allele : variant.genotype)
            if (allele < 0 || allele >= static_cast<int>(variant.alleles.size()))
                throw std::runtime_error("Incomplete or invalid GT in replay");
        if (reference.compare(variant.position, variant.alleles[0].size(), variant.alleles[0]) != 0)
            throw std::runtime_error("VCF REF does not match reference");
        std::map<std::string, State> next;
        for (const auto& entry : states) {
            State base = entry.second;
            ++visited;
            if (!Prune(base, normalized.at(variant.id).GetStart(), reference, pairs)) continue;
            State excluded = base;
            excluded.excluded.push_back(variant.id);
            Keep(next, std::move(excluded), pairs);
            const bool heterozygous = expected == 2 && variant.genotype[0] != variant.genotype[1];
            for (int orientation = 0; orientation < (heterozygous ? 2 : 1); ++orientation) {
                State included = base;
                if (heterozygous && !variant.phaseSet.empty()) {
                    const std::string phase = std::to_string(variant.sample) + ":" + variant.phaseSet;
                    const auto found = included.phase.find(phase);
                    if (found != included.phase.end() && found->second != orientation) continue;
                    included.phase[phase] = orientation;
                }
                bool valid;
                if (variant.sample == 2) {
                    valid = Apply(included.tracks[1], normalized.at(variant.id), orientation, reference);
                    if (pairs == 2)
                        valid = valid && Apply(included.tracks[3], normalized.at(variant.id), 1 - orientation, reference);
                } else {
                    const int track = variant.sample == 0 || !maternal ? 0 : 2;
                    valid = Apply(included.tracks[track], normalized.at(variant.id), orientation, reference);
                }
                if (valid) {
                    included.included.push_back(variant.id);
                    Keep(next, std::move(included), pairs);
                }
            }
        }
        states.swap(next);
    }
    const State* best = nullptr;
    for (auto& entry : states) {
        State& state = entry.second;
        if (!Prune(state, end, reference, pairs)) continue;
        bool finished = true;
        for (int i = 0; i < 2 * pairs; ++i) finished = finished && state.tracks[i].pending.empty();
        if (finished && (!best || state.included.size() > best->included.size())) best = &state;
    }
    if (!best) throw std::runtime_error("No complete mixed-ploidy replay path");
    return SMixedResult{best->included, best->excluded, visited};
}
}
