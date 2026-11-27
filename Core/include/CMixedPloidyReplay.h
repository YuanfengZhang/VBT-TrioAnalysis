#ifndef VBT_MIXED_PLOIDY_REPLAY_H
#define VBT_MIXED_PLOIDY_REPLAY_H
#include <string>
#include <vector>
namespace core {
// Haploid GTs contain exactly one allele, never a duplicated pseudo-diploid.
struct SMixedVariant {
    int id;
    int sample; // mother=0, father=1, child=2
    int position;
    std::vector<std::string> alleles;
    std::vector<int> genotype;
    std::string phaseSet;
};
struct SMixedModel {
    int motherPloidy;
    int fatherPloidy;
    int childPloidy;
};
struct SMixedResult {
    std::vector<int> included;
    std::vector<int> excluded;
    unsigned long statesVisited;
};
// Joint sequence replay: maternal and paternal constraints share both child
// paths. Caller-provided domains cannot cross model boundaries or gaps.
SMixedResult ReplayMixed(const std::string& reference, int begin, int end,
                        const SMixedModel& model,
                        std::vector<SMixedVariant> variants);
}
#endif
