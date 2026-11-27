#include "CMixedPloidyReplay.h"
#include <cassert>
#include <iostream>
using core::SMixedVariant;
int main() {
    const std::string ref = "AACCGGTTAAAAACGT";
    auto v = [](int id, int sample, int pos, std::vector<std::string> alleles,
                std::vector<int> gt, std::string ps = "") {
        return SMixedVariant{id, sample, pos, alleles, gt, ps};
    };
    auto r = core::ReplayMixed(ref, 0, ref.size(), {2, 0, 1}, {
        v(0, 0, 2, {"C", "T"}, {0, 1}), v(1, 2, 2, {"C", "T"}, {1})});
    assert(r.excluded.empty() && r.included.size() == 2);
    // Independent parental paths accept this invalid daughter.
    r = core::ReplayMixed(ref, 0, ref.size(), {2, 1, 2}, {
        v(0, 0, 2, {"C", "T"}, {1, 1}), v(1, 1, 2, {"C", "T"}, {1}),
        v(2, 2, 2, {"C", "T"}, {0, 1})});
    assert(!r.excluded.empty());
    r = core::ReplayMixed(ref, 0, ref.size(), {2, 1, 2}, {
        v(0, 0, 2, {"C", "T"}, {0, 1}), v(1, 1, 2, {"C", "G"}, {1}),
        v(2, 2, 2, {"C", "T", "G"}, {1, 2})});
    assert(r.excluded.empty());
    // Equivalent MNP and split SNVs.
    r = core::ReplayMixed(ref, 0, ref.size(), {0, 1, 1}, {
        v(0, 1, 2, {"CC", "TA"}, {1}), v(1, 2, 2, {"C", "T"}, {1}),
        v(2, 2, 3, {"C", "A"}, {1})});
    assert(r.excluded.empty());
    // Shifted deletion.
    r = core::ReplayMixed(ref, 0, ref.size(), {0, 1, 1}, {
        v(0, 1, 8, {"AA", "A"}, {1}), v(1, 2, 10, {"AA", "A"}, {1})});
    assert(r.excluded.empty());
    // Phased cis cannot be transmitted as trans.
    r = core::ReplayMixed(ref, 0, ref.size(), {2, 0, 1}, {
        v(0, 0, 2, {"C", "T"}, {0, 1}, "p"),
        v(1, 0, 3, {"C", "A"}, {0, 1}, "p"),
        v(2, 2, 2, {"C", "T"}, {1}), v(3, 2, 3, {"C", "A"}, {0})});
    assert(!r.excluded.empty());
    std::cout << "mixed replay: six native tests passed\n";
}
