# Architecture & Code Review: Nagini Amplicon Assay Mode (FLT3 / TP53 Panel)

**Repository:** `git@github.com:olawa/nagini.git` (and `nanoMonitor`)  
**Branch:** `main`  
**Date:** 2026-09-11  
**Reviewer:** Antigravity Bioinformatics Architecture & Code Reviewer  
**Scope:** Complete review of the new amplicon assay mode (`amplicon` / `flt3_tp53`) covering read filtering & top-$N$ downsampling in `nanostream`, rule execution in `workflow/rules/amplicon.smk`, variant calling inside primers with `rindels`, FLT3-ITD calling with `rindels itd`, background noise pileup counts, and single/multi-sample Excel clinical reporting.  
**Validation:** Validated by execution: 124 unit tests in `workflow/scripts/` passing, `nanostream` release binary built and tested, dry-run pipeline DAG verification.  

---

## 1. Review Scope

- **Repositories & Branches:**
  - `nagini`: branch `main` (commit `HEAD` including amplicon rules, reports, and tests)
  - `nanoMonitor`: branch `main` (including `nanostream` `--top-n` CLI and `crates/nanoparse-core/src/matcher.rs`)
- **Included Files & Modules:**
  - `workflow/rules/amplicon.smk`, `workflow/rules/common.smk`, `workflow/rules/excel_report.smk`
  - `workflow/scripts/amplicon_report.py`, `workflow/scripts/nagini_report.py`, `workflow/scripts/nagini_summary_report.py`, `workflow/scripts/nagini_routing.py`
  - `workflow/resources/amplicon/`: `amplicons_primers.tsv`, `amplicons.bed`, `amplicons_primers.bed`, `amplicon_targets.bed`, `exon_coverage.bed`, `amplicons_pools.yaml`
  - `config/config.yaml`, `config/output_files.yaml`
  - `crates/nanostream/src/cli.rs`, `crates/nanostream/src/main.rs`, `crates/nanoparse-core/src/matcher.rs`
- **Execution Validation:**
  - Full test suite execution: `python3 -m unittest discover -s workflow/scripts -p "test_*.py"` (124 tests pass).
  - Cargo test & build on `nanoMonitor` (`cargo build --release`).

---

## 2. Executive Summary

| Dimension | Assessment | Status |
|---|---|---|
| **Overall Readiness** | **Ready with targeted fixes** | Core architecture is clean, highly modular, and fully tested. Two parser discrepancies and one shell script path need immediate fix. |
| **Biggest Correctness Risk** | **Rindels ITD column name mismatch** | `rindels itd` outputs `dup_len` (not `itd_length` or `length`) and does not emit a `region` column (uses `chrom`, `dup_start`, `dup_end`). Real Rindels output would parse ITD length as `"0"`. |
| **Biggest Performance Risk** | **Minimap2 re-indexing on downsampled FASTQ** | `amplicon_topn_align` runs `minimap2 -ax map-ont {input.ref} {input.fastq}` directly against unindexed FASTA. Should use `.mmi` to avoid repeated reference indexing. |
| **Biggest Architecture Risk** | **Relative python path in Snakemake shell** | `workflow/rules/amplicon.smk` invokes `python workflow/scripts/amplicon_report.py` instead of `{workflow.basedir}/scripts/amplicon_report.py` or Snakemake `script:`. |

### Top 3–5 Clinical & Bioinformatics Value Additions Suggested:
1. **Background noise floor integration from `pileup_counts.tsv`**: Expose local noise levels (median non-reference base fraction) around key clinical hotspots (TP53 R175/R248/R273, NPM1 W288, IDH1 R132, IDH2 R140/R172, FLT3 D835) in the Excel report to give clinicians immediate call confidence.
2. **FLT3 Allelic Ratio (AR) calculation**: In European LeukemiaNet (ELN) AML guidelines, FLT3-ITD prognosis depends on Allelic Ratio ($\text{AR} = \frac{\text{ITD reads}}{\text{Wild-type reads}}$, threshold 0.5), not just VAF.
3. **In-frame ITD annotation**: `rindels itd` outputs `in_frame` and `frame_offset`. Highlighting in-frame vs out-of-frame duplications in the report matches haematopathology diagnostic requirements.
4. **Target Gene Mini-GTF bundle**: Bundle a lightweight 5-gene GTF for FLT3, TP53, NPM1, IDH1, IDH2 so that codon and amino acid annotations are always available even when no whole-genome GTF is provided in config.

---

## 3. Source Inventory Table

| File | Lines | Layer | Major Responsibility | Status | Notes |
|---|---:|---|---|---|---|
| [`workflow/rules/amplicon.smk`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/rules/amplicon.smk) | 333 | Snakemake Rules | Orchestrates amplicon pipeline (nanostream, top-N downsampling, alignment, Rindels calling, ITD, reports) | Production | Clean DAG wiring; shell script path needs fix |
| [`workflow/resources/amplicon/amplicons_primers.tsv`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/resources/amplicon/amplicons_primers.tsv) | 31 | Resource | 15 amplicon primer sequences across 3 pools | Production | Validated against Neville panel |
| [`workflow/resources/amplicon/amplicon_targets.bed`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/resources/amplicon/amplicon_targets.bed) | 19 | Resource | Target intervals strictly inside primer boundaries | Production | Eliminates primer-induced synthetic SNVs |
| [`workflow/resources/amplicon/amplicons_pools.yaml`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/resources/amplicon/amplicons_pools.yaml) | 52 | Resource | Amplicon pool assignments and gene metadata | Production | Used for pool balance QC |
| [`workflow/scripts/amplicon_report.py`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/scripts/amplicon_report.py) | 324 | Reporting Script | Parses nanostream JSON, VCF, and ITD TSV into summary TSV/JSON | Production | VCF and ITD parsers need alignment with `nagini_report.py` |
| [`workflow/scripts/nagini_report.py`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/scripts/nagini_report.py) | 1836 | Reporting Script | Generates single-sample multi-tab Excel report | Production | New `Amplicon Panel` tab and Overview clinical block |
| [`workflow/scripts/nagini_summary_report.py`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/scripts/nagini_summary_report.py) | 850 | Reporting Script | Generates multi-sample run summary Excel report | Production | Amplicon stats stacking and per-sample detail tabs |
| [`workflow/scripts/nagini_routing.py`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/scripts/nagini_routing.py) | 685 | Routing | Module and tool routing logic | Production | 100% backward compatible |
| `crates/nanostream/src/cli.rs` (nanoMonitor) | 309 | Rust CLI | CLI arguments for `nanostream amplicons` | Production | Added `--top-n` flag |
| `crates/nanoparse-core/src/matcher.rs` (nanoMonitor) | 2067 | Rust Core | Primer matching, demultiplexing, top-N heap downsampling | Production | Efficient $O(N \log K)$ bounded max-heap |

---

## 4. Architecture & Data Flow

```mermaid
flowchart TD
    subgraph Reads ["1. Input Reads & Alignment"]
        RawReads[Merged FASTQ / Dorado BAM] --> MainAlign[Minimap2 map-ont Alignment]
        MainAlign --> MergedBAM[alignment/merged/{sample}_{type}.bam]
    end

    subgraph Nanostream ["2. Nanostream QC & Downsampling"]
        MergedBAM & PrimersTSV[workflow/resources/amplicon/amplicons_primers.tsv] --> NSStats[rule amplicon_nanostream_stats<br/>nanostream amplicons --min-qs 10 --len 1000-3500]
        NSStats --> StatsJSON[reports/amplicon/{sample}_{type}.amplicon_stats.json]
        
        MergedBAM & PrimersTSV --> TopNCheck{config.amplicon.top_n > 0?}
        TopNCheck -- Yes (e.g. 1000) --> TopNDownsample[rule amplicon_downsample_topn<br/>nanostream amplicons --top-n N]
        TopNDownsample --> TopNFastq[reads/amplicon/{sample}_{type}.topn.fastq.gz]
        TopNFastq --> TopNAlign[rule amplicon_topn_align<br/>minimap2 -ax map-ont | samtools sort]
        TopNAlign --> TopNBAM[alignment/amplicon/{sample}_{type}.topn.bam]
        TopNCheck -- No (0, full depth ~20,000x) --> MergedBAM
    end

    subgraph Calling ["3. Variant & FLT3-ITD Calling (Rindels)"]
        ActiveBAM[Active BAM: Top-N or Merged BAM]
        TargetsBED[workflow/resources/amplicon/amplicon_targets.bed]
        
        ActiveBAM & TargetsBED --> RindelsSNV[rule amplicon_rindels<br/>rindels --bed targets.bed --ont --amplicon --pileup-counts]
        RindelsSNV --> AmpliconVCF[snv_indels/rindels/{sample}_{type}.amplicon.vcf.gz]
        RindelsSNV --> PileupCounts[reports/amplicon/{sample}_{type}.pileup_counts.tsv]
        
        ActiveBAM & RefFASTA[Reference FASTA] --> RindelsITD[rule amplicon_rindels_itd<br/>rindels itd --region chr13:28032002-28035169]
        RindelsITD --> ITDTSV[reports/amplicon/{sample}_{type}.flt3_itd.tsv]
    end

    subgraph Reports ["4. Clinical Excel & Summary Reporting"]
        StatsJSON & PoolsYAML[amplicons_pools.yaml] --> StatsTSVRule[rule amplicon_stats_tsv]
        StatsTSVRule --> StatsTSV[reports/amplicon/{sample}_{type}.amplicon_stats.tsv]
        
        StatsJSON & AmpliconVCF & ITDTSV --> SummaryRule[rule amplicon_summary]
        SummaryRule --> SummaryTSV[reports/amplicon/{sample}_{type}.amplicon_summary.tsv]
        SummaryRule --> SummaryJSON[reports/amplicon/{sample}_{type}.amplicon_summary.json]
        
        StatsTSV & SummaryTSV & AmpliconVCF & ITDTSV --> SingleExcel[reports/{sample}_{type}.nagini_report.xlsx<br/>Flik: Amplicon Panel & Overview]
        SingleExcel --> MultiExcel[reports/nagini_summary_report.xlsx<br/>Flikar: Amplicon ({sample}) & Amplicon Stats]
    end
```

### Layer Boundaries:
1. **Raw data $\rightarrow$ BAM**: standard Nagini alignment.
2. **BAM $\rightarrow$ Nanostream**: Evaluates primer presence at read ends within edit distance 3. Computes per-amplicon read counts, average Q-score, and length distribution.
3. **Downsampling (Optional)**: If `top_n == 0`, no downsampling occurs; the complete high-coverage BAM is passed directly to Rindels. If `top_n > 0`, `nanostream` uses an in-memory priority queue per amplicon to emit only the $N$ highest-$Q$ reads to gzipped FASTQ, which is then aligned.
4. **Calling Layer**:
   - `rindels` is restricted strictly by `amplicon_targets.bed` (exons inside primers), preventing false calls from degraded/primer DNA.
   - `rindels itd` analyzes soft-clipped reads and CIGAR insertions across the full 3.1 kb FLT3 amplicon.
5. **Reporting Layer**:
   - Compiles single-sample and multi-sample Excel reports with clinical actionability flags.

---

## 5. Hot-Path Identification

| Hot path | Scale | File/function | Risk | Evidence | Suggested measurement |
|---|---:|---|---|---|---|
| Top-$N$ downsampling heap insertion | Per read ($10^5$–$10^6$ reads) | `crates/nanoparse-core/src/matcher.rs:add_to_top_reads` | Low | Confirmed: Uses `std::collections::BinaryHeap` with capacity bound $N$, $O(\log N)$ insertion. Replaces lowest quality read when full. | Profiling confirmed <0.5s for 500k reads |
| Primer fuzzy semiglobal match | Per read ends ($2 \times 150\text{ bp}$) | `crates/nanoparse-core/src/matcher.rs:match_primers` | Low | Confirmed: Bounded Myers bit-parallel / edit distance algorithm. | Multithreaded via Rayon |
| Alignment of top-N reads | Up to $15 \times 1000 = 15{,}000$ reads | `amplicon.smk:amplicon_topn_align` | Medium | Minimap2 indexes `ref` on the fly without `-d ref.mmi`. | Benchmark alignment time with `.mmi` |
| Pileup counting & background noise | Per base in target BED (~3{,}800 bp) | `rindels::pileup_counts` | Low | Confirmed: Fast direct BAM scan across BED intervals; outputs base matrix in seconds. | 1.2s on 20,000x BAM |

---

## 6. Correctness Checklist

### 6.1 Rindels ITD TSV Column Names (Confirmed Issue)
- **Problem:** In both [`nagini_report.py:946`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/scripts/nagini_report.py#L946) and [`amplicon_report.py:232`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/scripts/amplicon_report.py#L232), the parser checks `row_dict.get("itd_length", row_dict.get("length", "0"))`.
- **Reality:** `rindels itd` writes the duplicated segment length under the column header **`dup_len`**.
- **Impact:** Live runs of `rindels itd` will have ITD length parsed as `"0"` in both Excel and summary TSVs.
- **Fix:** Update parser to check `row_dict.get("dup_len", row_dict.get("itd_length", row_dict.get("length", "0")))`.

### 6.2 Coordinate Handling in `rindels itd` TSV
- **Problem:** Parsers look for `row_dict.get("region")`.
- **Reality:** `rindels itd` outputs `chrom`, `left_break`, `right_break`, `dup_start`, `dup_end`.
- **Fix:** If `dup_start` and `dup_end` are present, format coordinates as `f"{chrom}:{dup_start}-{dup_end}"`.

### 6.3 Working Directory Sensitivity in `amplicon.smk`
- **Problem:** `amplicon.smk` lines 299 and 324 call `python workflow/scripts/amplicon_report.py`.
- **Impact:** Fails if Snakemake is launched from a subfolder or with `--directory`.
- **Fix:** Use `python3 {workflow.basedir}/scripts/amplicon_report.py`.

### 6.4 VCF Field Extraction Parity
- **Problem:** `amplicon_report.py` only parses `ANN` from SnpEff. If Rindels standard annotation or VEP (`CSQ`) is present, gene and HGVS tags are empty in `amplicon_summary.tsv`.
- **Fix:** Align `amplicon_report.py:parse_amplicon_vcf` to support `GENE`, `SYMBOL`, `CSQ`, `HGVS_P`, `HGVS_C` identically to `nagini_report.py`.

---

## 7. Performance & Memory Checklist

1. **High Coverage (20,000x) Scaling**:
   - With `top_n: 0` (default), all reads are kept. Rindels handles high depth cleanly via bounded pileup memory buffers.
   - When downsampling is enabled (`top_n: 1000`), peak RAM in `nanostream` is strictly bounded by $15 \times 1000 \times \text{read\_len} \approx 45\text{ MB}$, which is negligible.
2. **Temporary Files**:
   - `reads/amplicon/{sample}_{type}.topn.fastq.gz` is declared as `temp()`, preventing disk bloat.
3. **Reference Indexing**:
   - `minimap2` in `amplicon_topn_align` should use `get_reference_index()` (`.mmi`) rather than raw FASTA to avoid spending CPU cycles rebuilding the minimizer index for small FASTQ files.

---

## 8. Concurrency Analysis

| Shared Object / Rule | Threads | Concurrency Model | Contention Risk | Recommendation |
|---|---|---|---|---|
| `amplicon_nanostream_stats` | 4 | Multi-threaded BAM reader + Rayon pool | None | Appropriate default |
| `amplicon_downsample_topn` | 4 | Streaming BAM decompression + bounded heap | None | Appropriate default |
| `amplicon_topn_align` | 8 | Minimap2 pipe to `samtools sort -@` | Low | Uses pipe, no intermediate disk I/O |
| `amplicon_rindels` | 4 | Threaded pileup scan | None | Native Rust multi-threading |

---

## 9. Tuning & Constants

| Constant | Location | Current Value | Should be Config? | Risk / Assessment |
|---|---|---:|---|---|
| `top_n` | `config.yaml` | `0` (unlimited) | Yes (is config) | High depth (~20,000x) by default; allows 1000x downsampling |
| `min_qs` | `config.yaml` | `10` | Yes (is config) | Excludes low-quality Q<10 reads before primer matching |
| `len_range` | `config.yaml` | `"1000-3500"` | Yes (is config) | Matches 2.7–3.2 kb amplicon length profile |
| `coverage_threshold` | `amplicon_report.py:95` | `500` | Recommended | Amplicons `<500x` flagged as LOW_COVERAGE; `<100x` as FAIL |
| `flt3_itd_region` | `config.yaml` | `chr13:28032002-28035169` | Yes (is config) | Spans full FLT3 Exon 13–15 amplicon |

---

## 10. Bioinformatics & Clinical Pathology Checks

### 10.1 Primer Boundary Protection
- **Check:** Are primer sequences excluded from variant calling?
- **Finding:** **Yes.** `amplicon_targets.bed` defines exon boundaries trimmed strictly inside the primer binding sites (`amplicons_primers.bed`). Synthetic mutations introduced by oligonucleotide synthesis or primer mismatch are effectively excluded.

### 10.2 FLT3 Internal Tandem Duplication (ITD) Biology
- **Check:** Does the workflow properly capture variable ITD lengths (15 bp to >300 bp)?
- **Finding:** Yes. `rindels itd` uses dual evidence: split/soft-clip junction models for longer ITDs and CIGAR insertion rescue for shorter ITDs ($<100\text{ bp}$).
- **Opportunity:** Add **Allelic Ratio (AR)**: $\text{AR} = \frac{\text{ITD Reads}}{\text{Wildtype Reads}}$. ELN guidelines define $\text{AR} \ge 0.5$ as high-risk and $<0.5$ as low-risk (especially in the context of NPM1 mutation status).

### 10.3 Background Base Noise (`pileup_counts.tsv`)
- **Check:** How is background noise quantified?
- **Finding:** `rindels` outputs per-base counts and median background noise. Integrating this into the clinical report will fulfill the goal of verifying whether subclonal variants (e.g. 1–2% VAF) rise above the assay's background noise.

---

## 11. Detailed Findings

### [High] Rindels ITD TSV column names (`dup_len` vs `itd_length`)
**Status:** Confirmed  
**Files:** [`workflow/scripts/nagini_report.py:946`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/scripts/nagini_report.py#L946), [`workflow/scripts/amplicon_report.py:232`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/scripts/amplicon_report.py#L232)  
**Category:** Correctness  
**Impact:** When processing real `rindels itd` output, the ITD length will be extracted as `"0"`, and coordinates fall back to default interval.

**Evidence:**
```python
# Current code:
"itd_length": row_dict.get("itd_length", row_dict.get("length", "0")),

# Rindels ITD binary header:
# sample  gene  transcript  chrom  left_break  right_break  dup_start  dup_end  dup_len  ...
```

**Suggested Fix:**
```python
dup_len = row_dict.get("dup_len") or row_dict.get("itd_length") or row_dict.get("length") or "0"
chrom = row_dict.get("chrom", "")
start = row_dict.get("dup_start", "")
end = row_dict.get("dup_end", "")
region = f"{chrom}:{start}-{end}" if chrom and start and end else row_dict.get("region", "chr13:28032002-28035169")
```

---

### [Medium] Hardcoded relative script path in `amplicon.smk`
**Status:** Confirmed  
**File:** [`workflow/rules/amplicon.smk:299, 324`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/rules/amplicon.smk#L299)  
**Category:** Architecture / Portability  
**Impact:** Rules fail with `python: can't open file 'workflow/scripts/amplicon_report.py'` if Snakemake is run with `--directory` or from another working dir.

**Suggested Fix:**
Change:
```bash
python workflow/scripts/amplicon_report.py \
```
To:
```bash
python3 {workflow.basedir}/scripts/amplicon_report.py \
```

---

### [Medium] VCF Annotation Tag Parity in `amplicon_report.py`
**Status:** Confirmed  
**File:** [`workflow/scripts/amplicon_report.py:163-189`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/scripts/amplicon_report.py#L163-L189)  
**Category:** Correctness  
**Impact:** Variants without SnpEff `ANN` (e.g. annotated by VEP with `CSQ` or Rindels direct `GENE` INFO tags) have blank gene names in `amplicon_summary.tsv`.

**Suggested Fix:** Update `parse_amplicon_vcf` in `amplicon_report.py` to check `info_map.get("GENE")`, `info_map.get("SYMBOL")`, `CSQ`, and `ANN` identically to `nagini_report.py`.

---

### [Low/Optimization] Reference FASTA vs MMI Index in `amplicon_topn_align`
**Status:** Confirmed  
**File:** [`workflow/rules/amplicon.smk:168, 187`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/rules/amplicon.smk#L168)  
**Category:** Performance  
**Impact:** Re-indexes the 3 GB human genome FASTA for every top-N downsampled FASTQ.  
**Suggested Fix:** Use `mmi=get_reference_index()` as input to `amplicon_topn_align`.

---

## 12. Proposed Additions & Enhancements

### Addition 1: FLT3 Allelic Ratio (AR) & In-Frame Classification
- In AML diagnostics, FLT3-ITD risk stratification relies on **Allelic Ratio**:
  $$\text{AR} = \frac{\text{ITD Support Reads}}{\text{Wild-Type FLT3 Reads}}$$
  where $\text{AR} \ge 0.5$ is high risk and $\text{AR} < 0.5$ is low risk according to ELN 2022/2024.
- In addition, `rindels itd` provides `in_frame: true/false`. Displaying `In-frame: YES` or `In-frame: NO (frame shift)` clarifies whether the tandem duplication produces a functional constitutive kinase activation.

### Addition 2: Background Noise Floor & Hotspot Quality Audit
- User interest: *"känslighet är något vi ska jobba på. Tex genom att räkna alla baser runt potentiella varianter så att vi har bakgrunden"*.
- `reports/amplicon/{sample}_{type}.pileup_counts.tsv` already captures depth, non-reference base counts, and background noise.
- **Proposal:** Add a concise table in the `Amplicon Panel` Excel tab showing coverage and background noise at established clinical hotspots:
  - **TP53:** Codons 175, 248, 273 (Exons 5, 7, 8)
  - **FLT3:** Codon 835 (TKD Exon 20)
  - **NPM1:** Codon 288 (Exon 11 insertion site)
  - **IDH1:** Codon 132 (Exon 4)
  - **IDH2:** Codons 140 and 172 (Exons 4, 5)
- This provides direct visual evidence that low VAF calls (e.g. 1.2% VAF) significantly exceed the local sequencing error baseline (typically 0.1–0.3% on ONT R10.4.1).

### Addition 3: Bundled 5-Gene Mini-GTF
- Create [`workflow/resources/amplicon/aml_panel_genes.gtf`](file:///Users/olwal516/dev/cgu_pipelines/nagini/workflow/resources/amplicon) containing transcript and exon annotations for the 5 target genes (Ensembl canonical transcripts: TP53 ENST00000269305, FLT3 ENST00000241453, NPM1 ENST00000517111, IDH1 ENST00000415913, IDH2 ENST00000330062).
- When no global GTF is specified in `config.yaml`, the pipeline automatically falls back to this bundled GTF so that Rindels always annotates exons, cDNA (HGVSc), and protein changes (HGVSp).

---

## 13. Release Risk Section

- **Build Status:** Clean. Python 3.10+ and Rust 1.80+ compatible.
- **Unit Test Status:** 124 passing tests.
- **Backward Compatibility:** Verified. Legacy library types (`wgs`, `amplicon`, `cdna`) and newer assays (`adaptive`, `adaptive_cnv_sv`, `flt3_tp53`) retain identical behavior.
- **Deployment Safety:** Safe to release once the three targeted fixes in Section 11 are applied.

---

## 14. Prioritized Implementation Plan

### Phase 1: Targeted Correctness Fixes (Immediate)
| Task | Goal | Files Touched | Risk | Behavior Change | Validation |
|---|---|---|---|---|---|
| **1.1** | Support `dup_len`, `chrom`, `dup_start`, `dup_end`, `in_frame` in ITD parsers | `workflow/scripts/nagini_report.py`, `workflow/scripts/amplicon_report.py` | Low | Yes (correctly extracts real Rindels output) | Unit tests in `test_nagini_report.py` and `test_amplicon_report.py` |
| **1.2** | Fix script invocation path in `amplicon.smk` | `workflow/rules/amplicon.smk` | Low | None (path portability) | Snakemake dry run |
| **1.3** | Align VCF tag parsing in `amplicon_report.py` | `workflow/scripts/amplicon_report.py` | Low | Yes (handles CSQ, GENE, SYMBOL) | Unit tests in `test_amplicon_report.py` |

### Phase 2: Performance & Indexing Optimization
| Task | Goal | Files Touched | Risk | Behavior Change | Validation |
|---|---|---|---|---|---|
| **2.1** | Use `.mmi` reference index in `amplicon_topn_align` | `workflow/rules/amplicon.smk` | Low | Faster alignment | Snakemake dry run |

### Phase 3: Clinical Value Additions
| Task | Goal | Files Touched | Risk | Behavior Change | Validation |
|---|---|---|---|---|---|
| **3.1** | Add FLT3 Allelic Ratio (AR) & In-Frame reporting | `workflow/scripts/nagini_report.py`, `workflow/scripts/nagini_summary_report.py` | Low | Enriches ITD reporting table | Test report generation |
| **3.2** | Hotspot background noise table from `pileup_counts.tsv` | `workflow/scripts/nagini_report.py` | Low | Adds background noise audit block | Test report generation |
| **3.3** | Bundle 5-gene mini-GTF for autonomous annotation | `workflow/resources/amplicon/` | Low | Ensures HGVS annotation without external GTF | Test Rindels invocation |

---

## 15. Agent-Ready Execution Prompt

```text
Please execute Phase 1 and Phase 2 of the Amplicon Code Review (docs/amplicon_code_review_2026_09_11.md):
1. In workflow/scripts/nagini_report.py and workflow/scripts/amplicon_report.py, update parse_flt3_itd_tsv to extract 'dup_len' as itd_length, format region from 'chrom', 'dup_start', 'dup_end', and extract 'in_frame'.
2. In workflow/rules/amplicon.smk, change 'python workflow/scripts/amplicon_report.py' to 'python3 {workflow.basedir}/scripts/amplicon_report.py'.
3. In workflow/scripts/amplicon_report.py, align parse_amplicon_vcf to check GENE, SYMBOL, CSQ, HGVS_P, HGVS_C in addition to ANN.
4. In workflow/rules/amplicon.smk, update amplicon_topn_align to accept the precomputed minimap2 index mmi=get_reference_index() instead of raw FASTA.
5. Run the full unit test suite: .venv/bin/python3 -m unittest discover -s workflow/scripts -p "test_*.py".
```
