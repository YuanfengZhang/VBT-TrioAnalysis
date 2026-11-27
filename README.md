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
./vbt mendelian -mother Example/UG_CEU_merged_cleaned_chr21.vcf -father Example/UG_CEU_merged_cleaned_chr21.vcf -child Example/UG_CEU_merged_cleaned_chr21.vcf -ped Example/ceu.ped --outDir Example/output -out-prefix ceu_21_sample
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

Build with matching HTSlib headers and library (tested with HTSlib 1.22.1):

```sh
make -B -j4 HTSLIB_PREFIX=/path/to/htslib-install
```

Do not link the bundled macOS libraries on Linux. The executable is `./vbt`.

```sh
./vbt mendelian -mother mother.vcf -father father.vcf -child child.vcf \
  -ref reference.fa -outDir output -out-prefix trio \
  --child-sex male --par-dir assembly-bed-directory --ypar-masked
```

`--child-sex` accepts `male`, `female` or `unknown`. It describes only the
child; the model assumes an XX mother and XY father. `--ypar-masked` is a
pure flag, default false; `--ypar-masked=false` is an error. It declares how
calls were generated and does not verify their provenance or modify a FASTA.
Without any new policy option, legacy behavior is unchanged.

The resource directory must contain exactly one file for each suffix
`chrX_PAR.bed`, `chrX_nonPAR.bed`, `chrY_PAR.bed`, `chrY_nonPAR.bed`, optionally
with `.gz`. Coordinates are zero-based, half-open, nonnegative and nonoverlapping.
Gaps remain unclassified. An optional `par_regions.tsv` with columns
`chrom`, `start`, `end`, `region` labels canonical intervals `PAR1` or `PAR2`.
It names regions only and never expands BED eligibility. Without the manifest,
PAR fragments retain generic X_PAR/Y_PAR names rather than ordinal labels.
Variant reference spans must fit in one BED interval; pure insertions use their
anchor. XY aliases X/chrX and Y/chrY are accepted.

Masked X PAR uses mother2/father2/child2. Unmasked X PAR is not evaluated;
Y PAR is always excluded. Male X nonPAR uses mother2/child1, male Y nonPAR
uses father1/child1, and female X nonPAR uses mother2/father1/child2 with
separate child paths for the two parental contributions. Homozygous pseudo2
calls are folded only in the computation; original GTs remain unchanged.
Haploid heterozygotes and explicit missing required GTs are not evaluated.
Nonparticipating parents do not determine completeness or the nonreference gate.

The mixed-ploidy engine reuses VBT CVariant trimming, the original
CBaseVariantProvider::FindOptimalTrimmings ordinary/ambiguous context
selection and CHaplotypeSequence playback. The default end-first trim and
neighbor-driven opposite trim remain native components. Joint state dominance retains the path with most included source
records; ties retain the first deterministic state. Excluded source records
are representation-dependent matching attribution, not de novo localization.
A record-free interval in a variant-only VCF replays reference as in legacy VBT;
this does not establish a known, complete reference GT. Explicit missing GTs
are retained and form replay barriers. BED gaps and model boundaries also
separate replay domains. Phase sets, including phased GTs without PS, constrain
orientation jointly within a domain. No new sitewise transmission matcher or
repeat-window normalization is applied.

Autosomes use the original provider, four GT/AM comparisons, decider, writer
and reports. XY input is routed away before the legacy two-slot GT decoder.
In explicit mode, `PREFIX_trio.vcf` and legacy logs therefore contain only
autosomes. Consumers must also read all three new XY products:

- `PREFIX_mother_xy.vcf`
- `PREFIX_father_xy.vcf`
- `PREFIX_child_xy.vcf`
- `PREFIX_xy_summary.tsv`

Each source VCF preserves original records, alleles, all sample columns and GTs.
INFO `XY_SOURCE_ID` is `mother:N`, `father:N` or `child:N`, where N is the
one-based record ordinal in the original input, including non-XY records.
INFO `XY_STATUS` is consistent, violation, not_evaluable or unsupported;
`XY_REGION` and `XY_REASON` identify native policy decisions. Only
female_y_nonpar, unsupported_variant and unsupported_replay map to unsupported;
other policy exclusions are not_evaluable. Consistent/violation are actual
replay decisions. The summary header is
`source\tregion\tstatus\treason\trecords` and counts source records, not
unique trio events. Never sum three source counts as a sitewise denominator.
A child-record denominator uses only the child source; parent-only events
must be reported separately. An empty evaluated denominator has no rate.

Run portable tests with an isolated scratch directory supplied explicitly:

```sh
make test-core HTSLIB_PREFIX=/path/to/htslib-install
python3 -I tests/native_xy.py ./vbt /path/to/scratch/xy
python3 -I tests/native_representations.py ./vbt /path/to/baseline-vbt /path/to/scratch/representations
```
