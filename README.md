# VARIANT BENCHMARKING TOOLS (VBT)

**Author:** Berke Cagkan Toptas (berke.toptas@sbgdinc.com) <br />
**Project Language:** C++11 <br />
<br />
VBT provides a set of tools that is used for aligner/variant calling benchmarking. VBT uses [vcfeval](https://github.com/RealTimeGenomics/rtg-tools) as its core variant comparison library and [htslib](https://htslib.org) to read/write VCF and FASTA files.
<br />
<br />
Preprint of VBT is published at Biorxiv: [Comparing complex variants in family trios](https://doi.org/10.1101/253492)  <br/>
<br />
VBT contains 2 tool:

   1. [VARIANT COMPARISON](DuoComparison/README.md)
   2. [MENDELIAN VIOLATION DETECTOR](MendelianViolation/README.md)

## Installing VBT:


#### Using makefile:

VBT uses htslib as an external dependency and it needs to be installed before VBT. Please download the latest htslib version from [htslib official page](https://htslib.org). After you download and extract htslib, you can install it using the following steps: (following steps are tested with htslib 1.6)

```
cd htslib-1.x
./configure
make
make install
```

You can download the latest version of VBT using the following git command:

```
git clone https://github.com/sbg/VBT-TrioAnalysis.git
```

Once you have successfully installed htslib and download VBT, Then you can use following commands to compile VBT:

```
cp -R (HTSLIB_PATH)/htslib (VBT_PATH)          //Copy htslib folder to VBT directory
cp (HTSLIB_PATH)/libhts.a (VBT_PATH)/lib       //Copy htslib binaries to VBT/lib
cp (HTSLIB_PATH)/libhts.so (VBT_PATH)/lib      //Copy htslib binaries to VBT/lib
cp (HTSLIB_PATH)/libhts.so.2 (VBT_PATH)/lib    //Copy htslib binaries to VBT/lib
cd (VBT_PATH)                                  //Go inside the VBT folder where makefile is
make all                                       //Compile VBT
```
#### Using dockerfile:

You can run VBT as a docker image. A dockerfile is added to the project that is running on Ubuntu 14.04 instance.

```
cd (VBT_PATH)
docker build -t vbtapp:v1.0 .
docker run -ti vbtapp:v1.0
```

**Note:** If you are unable to compile program or getting runtime error, please contact me. **(berke.toptas@sbgdinc.com)**

## Parameter format:

```
./vbt varComp [PARAMETERS]
./vbt mendelian [PARAMETERS]
```
For additional help, please use the following commands:

```
./vbt --help
./vbt [module_name] --help
```

## Sample Execution (mendelian):

**Example** folder contains a minimal example for testing VBT mendelian. We placed chr21 of merged BWA+UG CEPH trio as an input VCF. For the FASTA file, we extract chr21 from GRCh37 reference. In order to run the sample test, following command can be executed :

```
./vbt mendelian -mother Example/UG_CEU_merged_cleaned_chr21.vcf -father Example/UG_CEU_merged_cleaned_chr21.vcf -child Example/UG_CEU_merged_cleaned_chr21.vcf -pedigree Example/ceu.ped -outDir Example/output -out-prefix ceu_21_sample
```
This sample execution will produce:

	1. ceu_21_sample_BestPathLogs.txt              //Logs of Best Path Algorithm
	2. ceu_21_sample_ChildReportLog.txt            //Mendelian decision for Non-Ref child variants only
	3. ceu_21_sample_DetailedLogs.txt              //Detailed human-readable logs of trio concordance analysis
	4. ceu_21_sample_tab_delim_detailed_log.tsv    //Tab delimited version of detailed logs
	5. ceu_21_sample_trio.vcf                      //Merged output trio (Mendelian decisions are annotated for each record)

For additional parameters, please check mendelian paramater specifications.

## Source Code Documentation:

A [Doxygen html documentation](Doxygen) is available for VBT under Doxygen folder.



## Explicit native XY mode

### 1. XY Mendelian inheritance model

The model assumes an XX mother and XY father. Child sex determines which parental GTs are compared with the child and how many alleles each GT must contain outside the pseudoautosomal regions (PARs). Inherited haplotypes must reproduce the child sequence; equivalent VCF representations can describe the same transmission.

| Region and condition | Required ploidies |
| --- | --- |
| X PAR, Y-PAR-masked calling declared | Mother 2, father 2, child 2 |
| Male X nonPAR | Mother 2, child 1; father not required |
| Male Y nonPAR | Father 1, child 1; mother not required |
| Female X nonPAR | Mother 2, father 1, child 2 |
| Female Y nonPAR | Excluded |
| Y PAR | Excluded for either child sex |

Female X nonPAR requires maternal and paternal contributions on different child haplotypes. X PAR is evaluated only when the caller declares that calls were generated with Y PAR masked. Y PAR is excluded to avoid duplicate assessment. VBT does not infer this calling policy from genotypes.

Required GTs must be complete and have the expected ploidy. A homozygous two-slot representation of a haploid call can be folded internally; a heterozygous two-slot call cannot supply one haplotype. An explicit missing GT in the child or a parent needed for the comparison is not reference evidence. GTs from other parents do not affect whether the comparison has complete calls or a nonreference allele. Unknown child sex, unclassified intervals and events crossing boundaries between inheritance rules are not evaluated. Pure insertions use their reference anchor for classification.

### 2. Upstream implementation at the pinned baseline

These excerpts describe upstream-derived commit `0aece17`, before this branch's changes. No equivalence between this commit and a dated container tag is assumed.

The [VCF reader](https://github.com/YuanfengZhang/VBT-TrioAnalysis/blob/0aece17bf617c3920745c732d407830253eebe8f/VcfIO/src/CVcfReader.cpp) recognizes one- and two-slot GTs:

```cpp
            bcf_get_genotypes(m_pHeader, m_pRecord, &gt_arr, &ngt_arr);
            zygotCount = ngt_arr / samplenumber;
            a_pVariant->m_nAlleleCount = zygotCount;
            if(zygotCount == 2)
                a_pVariant->m_bIsPhased = bcf_gt_is_phased(gt_arr[0]) || bcf_gt_is_phased(gt_arr[1]);
            else if(zygotCount == 1)
                a_pVariant->m_bIsPhased = bcf_gt_is_phased(gt_arr[0]);
```

The [Mendelian provider](https://github.com/YuanfengZhang/VBT-TrioAnalysis/blob/0aece17bf617c3920745c732d407830253eebe8f/MendelianViolation/src/CMendelianVariantProvider.cpp) then removes called hom-reference and non-diploid records:

```cpp
        if(!variant.m_bIsNoCall && CUtils::IsHomRef(variant))
            continue;

        //Eliminate variants rather than diploid
        if(variant.m_nZygotCount != 2)
            continue;
```

Its `SetCommonChromosomes` requires matching contig names in all three inputs and sufficient retained records in every role:

```cpp
                if(childItr->first == motherItr->first && childItr->first == fatherItr->first)
                {
                    if(m_aChildVariantList[childItr->second].size() > LEAST_VARIANT_THRESHOLD
                       &&
                       m_aFatherVariantList[fatherItr->second].size() > LEAST_VARIANT_THRESHOLD
                       &&
                       m_aMotherVariantList[motherItr->second].size() > LEAST_VARIANT_THRESHOLD)
```

[Constants.h](https://github.com/YuanfengZhang/VBT-TrioAnalysis/blob/0aece17bf617c3920745c732d407830253eebe8f/Constants.h) defines:

```cpp
const int LEAST_VARIANT_THRESHOLD = 2;
```

Each role therefore needs at least three retained records on a shared contig. For example, a father with only hom-reference calls can remove a relevant male-X contig from the legacy analysis. A successful exit with a header-only VCF does not prove that replay processed any records.

Each [CSemiPath](https://github.com/YuanfengZhang/VBT-TrioAnalysis/blob/0aece17bf617c3920745c732d407830253eebe8f/Core/include/CSemiPath.h) has two haplotype tracks:

```cpp
    CHaplotypeSequence m_haplotypeA;
    CHaplotypeSequence m_haplotypeB;
```

The [analyzer](https://github.com/YuanfengZhang/VBT-TrioAnalysis/blob/0aece17bf617c3920745c732d407830253eebe8f/MendelianViolation/src/CMendelianAnalyzer.cpp) runs father-child genotype matching (GT), father-child allele matching (AM), mother-child GT, and mother-child AM. AM processes records excluded by the preceding GT pass. These assignments occur at separate points in `ProcessChromosome`:

```cpp
        m_aBestPathsFatherChildGT[triplet.m_nTripleIndex] = replayFatherChildGT.FindBestPath(ctg, true);
        m_aBestPathsFatherChildAM[triplet.m_nTripleIndex] = replayFatherChildAM.FindBestPath(ctg, false);
        m_aBestPathsMotherChildGT[triplet.m_nTripleIndex] = replayMotherChildGT.FindBestPath(ctg, true);
        m_aBestPathsMotherChildAM[triplet.m_nTripleIndex] = replayMotherChildAM.FindBestPath(ctg, false);
```

The [Mendelian decider](https://github.com/YuanfengZhang/VBT-TrioAnalysis/blob/0aece17bf617c3920745c732d407830253eebe8f/MendelianViolation/src/CMendelianDecider.cpp) combines those pairwise results. The provider filter and two-track replay do not implement the required haploid models. Independently matching each parent also cannot by itself enforce different maternal and paternal child paths in female X nonPAR. Removing the haploid filter alone leaves these assumptions in place.

### 3. Changes in this XY branch and usage

The accepted implementation is commit `60c1a3d`. Its changed and added files have the following responsibilities:

| Files | Responsibility |
| --- | --- |
| `MendelianViolation/src/CMendelianAnalyzer.cpp` | Dispatch explicit XY arguments before legacy analysis and expose the new CLI in help. |
| `VcfIO/src/CVcfReader.cpp` | Skip X/Y before the legacy fixed-slot decoder when running the autosome-only path. |
| `MendelianViolation/include/CExplicitXY.h`, `MendelianViolation/src/CExplicitXY.cpp` | Read raw source records and actual GT slots; choose which parental GTs to compare and which inheritance rule to use from child sex and the BED intervals; check whether those GTs can be evaluated; run replay separately in each continuous interval, splitting at missing or invalid calls; write source VCFs and counts; invoke the legacy autosome analyzer. |
| `Core/include/CMixedPloidyReplay.h`, `Core/src/CMixedPloidyReplay.cpp` | Keep parental and child sequences in one replay state, read the alleles present in each GT, enforce phase and separate maternal and paternal child haplotypes, and score paths by the number of included source records. |
| `Base/CBaseVariantProvider.h` | Make the existing stateless `FindOptimalTrimmings` callable through the adapter as a protected static method; its implementation remains unchanged. |
| `Makefile` | Select matching HTSlib headers and libraries through `HTSLIB_PREFIX` and add the core test target. |
| `tests/mixed_replay.cpp` | Test actual haploid replay, female child-path separation, equivalent representations and phase constraints. |
| `tests/native_xy.py` | Exercise policy, ploidy, missing-GT priority, source IDs and preservation of original GTs through the CLI. |
| `tests/native_representations.py` | Exercise MNP/split-SNV, shifted indel, overlapping trimming and phase cases; compare nonempty autosome outputs against the baseline. |

`CMendelianVariantProvider.cpp` remains unchanged. Explicit XY input uses the provider adapter in `CExplicitXY.cpp`; the original provider still supplies autosome analysis. The code excerpts below are from accepted commit `60c1a3d`, separate from the upstream excerpts in subsection 2.

The explicit wrapper removes only its three policy options and adds `--autosome-only` before calling a fresh legacy analyzer. This preserves the original autosome provider, comparisons, decider and output path:

[MendelianViolation/src/CExplicitXY.cpp](https://github.com/YuanfengZhang/VBT-TrioAnalysis/blob/60c1a3d2dec608d343c2c2a94b5f5becfed42615/MendelianViolation/src/CExplicitXY.cpp):

```cpp
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
```

Normalization separates records with one trimming choice from those with multiple choices for each sample, then calls the existing trimmer to choose using nearby records. The adapter does not replace its implementation:

[Core/src/CMixedPloidyReplay.cpp](https://github.com/YuanfengZhang/VBT-TrioAnalysis/blob/60c1a3d2dec608d343c2c2a94b5f5becfed42615/Core/src/CMixedPloidyReplay.cpp):

```cpp
        SConfig config;
        config.m_bTrimBeginningFirst = false;
        ContextTrimmer::FindOptimalTrimmings(ambiguous, &ordinary, config);
```

Selected alleles are then replayed through the original single-haplotype sequence component. `slot` refers to an actual genotype slot, including the sole slot of a haploid GT:

[Core/src/CMixedPloidyReplay.cpp](https://github.com/YuanfengZhang/VBT-TrioAnalysis/blob/60c1a3d2dec608d343c2c2a94b5f5becfed42615/Core/src/CMixedPloidyReplay.cpp):

```cpp
    CHaplotypeSequence sequence(reference.data(), static_cast<int>(reference.size()));
    sequence.AddVariant(COrientedVariant(native, slot));
    sequence.MoveForward(start);
    while (!sequence.IsOnTemplate() || sequence.GetTemplatePosition() < stop) {
        track.pending += sequence.NextBase();
        sequence.Next();
    }
    track.cursor = stop;
```

Replay compares one parent-child sequence pair for male X or Y nonPAR, and two pairs when both parents contribute. The child must have exactly that many actual genotype slots:

[Core/src/CMixedPloidyReplay.cpp](https://github.com/YuanfengZhang/VBT-TrioAnalysis/blob/60c1a3d2dec608d343c2c2a94b5f5becfed42615/Core/src/CMixedPloidyReplay.cpp):

```cpp
    const bool maternal = model.motherPloidy != 0;
    const bool paternal = model.fatherPloidy != 0;
    const int pairs = maternal && paternal ? 2 : 1;
    if ((!maternal && !paternal) || model.childPloidy != pairs)
        throw std::runtime_error("Invalid mixed-ploidy transmission model");
```

When both parents contribute, child slot `orientation` goes to one sequence pair and slot `1 - orientation` goes to the other. Both assignments occur in the same state; a female-X father and mother cannot each claim the same child path:

[Core/src/CMixedPloidyReplay.cpp](https://github.com/YuanfengZhang/VBT-TrioAnalysis/blob/60c1a3d2dec608d343c2c2a94b5f5becfed42615/Core/src/CMixedPloidyReplay.cpp):

```cpp
                if (variant.sample == 2) {
                    valid = Apply(included.tracks[1], normalized.at(variant.id), orientation, reference);
                    if (pairs == 2)
                        valid = valid && Apply(included.tracks[3], normalized.at(variant.id), 1 - orientation, reference);
                } else {
                    const int track = variant.sample == 0 || !maternal ? 0 : 2;
                    valid = Apply(included.tracks[track], normalized.at(variant.id), orientation, reference);
                }
```

Build with matching HTSlib headers and library (tested with HTSlib 1.22.1):

```sh
make -B -j4 HTSLIB_PREFIX=/path/to/htslib-install
```

Do not link the bundled macOS libraries on Linux. The executable is `./vbt`.

```sh
./vbt mendelian -mother mother.vcf -father father.vcf -child child.vcf \
  -ref reference.fa -outDir output -out-prefix trio \
  --child-sex male --par-dir assembly-bed-directory

# Add --ypar-masked to the command only for calls generated with Y PAR masked.
```

#### Options and regions

`--child-sex` accepts `male`, `female` or `unknown`; the model assumes an XX mother and XY father. Add `--ypar-masked` only for calls generated with Y PAR masked. Otherwise omit it. `--ypar-masked=false` is invalid. The flag does not verify calling history or modify the FASTA. Without the new options, VBT keeps its legacy behavior.

The BED directory needs one file per suffix: `chrX_PAR.bed`, `chrX_nonPAR.bed`, `chrY_PAR.bed` and `chrY_nonPAR.bed`, optionally gzip-compressed. Coordinates are zero-based, half-open, nonnegative and nonoverlapping. X/chrX and Y/chrY are accepted aliases. Gaps remain unclassified. Each variant's REF span must fit within one interval; insertions use their anchor.

Optional `par_regions.tsv` columns are `chrom`, `start`, `end`, `region`, with PAR1 or PAR2 as the region name. This file only labels existing BED coverage. Without it, labels are X_PAR/Y_PAR, never numbered by BED row order.

#### Inheritance and replay

- Masked X PAR uses two allele copies in each parent and child. Unmasked X PAR is skipped; Y PAR is always skipped.
- Male X nonPAR checks the child's allele against the mother. Male Y nonPAR checks it against the father.
- Female X nonPAR requires one maternal and one paternal contribution to different child haplotypes.
- Only participating people need complete, valid calls. Homozygous diploid calls can be treated internally as haploid where appropriate, without changing output GTs. Haploid heterozygotes and missing required calls are excluded.

VBT compares the resulting haplotype sequences. It reuses `CVariant` trimming, `CBaseVariantProvider::FindOptimalTrimmings` and `CHaplotypeSequence` playback. Nearby records help choose trimming, but only within the same continuous BED interval and inheritance rule, without crossing missing or invalid calls. Legacy trimming can use records across the contig. No new sitewise matcher or repeat-window normalization is added.

Replay keeps actual haploid alleles, without duplicating them into `1/1`. Female X contributions are checked together, so both parents cannot match the same child haplotype. Phased calls and phase sets constrain the paths within each replay interval, including phased calls without PS.

When paths reach the same sequence state, VBT keeps the one containing the most source records; ties keep the first deterministic state. Which records are excluded can depend on their representation. A violation does not establish a de novo mutation or identify a unique causal allele.

A region with no record in a variant-only VCF replays reference sequence, as in legacy VBT; it is not evidence of a complete reference GT. An explicit missing GT stops replay, as do BED gaps and changes in inheritance rule.

#### Outputs and counting

Autosomes retain the original provider, comparisons, decider and reports. In the new mode, XY records bypass the legacy two-allele reader, so `PREFIX_trio.vcf` and legacy logs contain only autosomes. Read XY results from:

- `PREFIX_mother_xy.vcf`
- `PREFIX_father_xy.vcf`
- `PREFIX_child_xy.vcf`
- `PREFIX_xy_summary.tsv`

Source VCFs preserve the original records, alleles, sample columns and GTs. `XY_SOURCE_ID` links each result to its input role and one-based original record number, including non-XY records in that numbering.

`XY_STATUS` is `consistent` or `violation` for replay results, and `not_evaluable` or `unsupported` for excluded records. `XY_REGION` gives the region and `XY_REASON` explains the decision. Female Y nonPAR and unsupported variants or replay use `unsupported`; other exclusions use `not_evaluable`.

The summary columns are `source`, `region`, `status`, `reason`, `records`. Counts refer to each input source, not unique trio events. Use the child source for a child-record denominator, and report parent-only records separately. Do not add all three sources into one denominator. An empty evaluated denominator has no rate. The model does not cover aneuploidy or infer parental karyotypes.

#### Sample selection and legacy options

For multisample VCFs, select each person with `-sample-mother`, `-sample-father`, `-sample-child`, or a valid `-pedigree` file. Only single-sample inputs can omit selection without a PED. The PED must identify one unambiguous child and agree with the declared sex.

Legacy `-no-call`, `-filter` and `--disable-ref-overlap` options still apply to autosomes. They do not change XY behavior: XY records are not selected by FILTER, trimming stays enabled, and missing XY calls remain excluded even with `-no-call none`.

Run portable tests with an isolated scratch directory supplied explicitly:

```sh
make test-core HTSLIB_PREFIX=/path/to/htslib-install
python3 -I tests/native_xy.py ./vbt /path/to/scratch/xy
python3 -I tests/native_representations.py ./vbt /path/to/baseline-vbt /path/to/scratch/representations
```

## Usage

Clone this fork's XY branch and build against an installed HTSlib prefix containing matching headers and libraries. A C++11 compiler, make, zlib and pthread support are required.

```sh
git clone --branch XY https://github.com/YuanfengZhang/VBT-TrioAnalysis.git
cd VBT-TrioAnalysis
make -B -j4 HTSLIB_PREFIX=/path/to/htslib-install
./vbt mendelian --help
mkdir -p inputs resources/GRCh38 resources/CHM13 output
```

Provide a coordinate-sorted VCF for each role as `inputs/mother.vcf`, `inputs/father.vcf` and `inputs/child.vcf`. These may be independently called files with different records and variant representations. VBT does not require a merged trio VCF or a PED file. Use the same reference assembly and compatible contig names for all three files. For a multisample VCF, explicitly select the intended sample with `-sample-mother`, `-sample-father` or `-sample-child`, or supply a valid PED that selects it. Omit sample selection only for a single-sample VCF; explicit XY mode rejects an unselected multisample input. The sample names in the example below are placeholders for the actual VCF sample names.

An optional PED can select samples through `-pedigree family.ped`. In explicit XY mode it must identify exactly one child with both parental IDs, those IDs must match the VCF samples, and its sex must agree with the declared child sex. Do not add a PED merely to run three single-sample VCFs.

Use `bedtools maskfasta` with the assembly-matched Y PAR BED to generate a YPARmasked FASTA, then index it with `samtools faidx`. Keep the reference assembly and contig names consistent with the VCF and BED files. Use the exact FASTA used to generate the calls, including its masking choices.

Obtain the four XY BED files (X PAR, X nonPAR, Y PAR and Y nonPAR) from the assembly-specific XY directory in the official [GIAB genome stratifications v3.6 release](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.6/). Use the directory matching the VCF and reference assembly, and keep assemblies in separate directories. VBT reads the gzip BED files directly.

`par_regions.tsv` is optional and is not one of these four GIAB BED downloads. If canonical PAR1/PAR2 labels are needed, supply an assembly-specific file derived from that assembly's authoritative PAR coordinates, with four tab-separated header fields (literal tab characters between `chrom`, `start`, `end` and `region`) and zero-based, half-open coordinates; accepted region labels are `PAR1` and `PAR2`. It only names intervals already allowed by the BED files. Without it, VBT uses generic X_PAR/Y_PAR labels. Do not label BED fragments by their order or pair X and Y fragments by row number.

Run the three role VCFs with the appropriate child sex, reference and BED directory. Create the output directory before running. This example uses GRCh38, a male child and explicit sample selection:

```sh
./vbt mendelian \
  -mother inputs/mother.vcf -sample-mother MOTHER \
  -father inputs/father.vcf -sample-father FATHER \
  -child inputs/child.vcf -sample-child CHILD \
  -ref resources/GRCh38/reference.fa \
  -outDir output -out-prefix trio \
  --child-sex male --par-dir resources/GRCh38
```

For a female child, use `--child-sex female`. For CHM13 calls, use the matching CHM13 FASTA and `--par-dir resources/CHM13`. If calls were generated with Y PAR masked, add `--ypar-masked` to the command. For example, the following command declares masked calling; its reference must be the matching reference used for those calls, with the same masking choices:

```sh
./vbt mendelian \
  -mother inputs/mother.vcf -father inputs/father.vcf -child inputs/child.vcf \
  -ref resources/calling-reference.fa \
  -outDir output -out-prefix trio_masked \
  --child-sex female --par-dir resources/GRCh38 --ypar-masked
```

Omit the sample options in the second command only when each role has a single-sample VCF, or a valid PED selects the samples. Multisample inputs require explicit sample selection or that PED. The flag declares how calls were generated; it does not mask the FASTA or repair existing calls. For unmasked calls, omit the flag entirely. `--ypar-masked=false` is invalid. Unmasked X PAR and all Y PAR records are excluded from evaluation.

Read `output/trio_mother_xy.vcf`, `output/trio_father_xy.vcf` and `output/trio_child_xy.vcf` for the original XY records annotated with `XY_STATUS`, `XY_REGION`, `XY_REASON` and `XY_SOURCE_ID`. `output/trio_xy_summary.tsv` counts records separately for each source. `output/trio_trio.vcf` and the original report files contain only autosome results in explicit XY mode; they do not include the XY records. Use the child source for a child-record denominator, report parent-only records separately, and do not add the three source counts as though they were unique trio events.
