# Code Review: Nagini Architecture, Routing, Excel Reporting & ClairS-TO Tuning

**Repository:** `git@github.com:olawa/nagini.git`  
**Branch:** `main`  
**Commits Reviewed:** `17c14ef`, `351694b`, `dd9ab8e` (plus BCR-ABL1 extensions)  
**Reviewer:** Antigravity Code Review Agent  
**Goal:** Review recent architectural updates (sample routing, reference registry, Excel migration, ClairS-TO parameter derivation, demux QC) for correctness, maintainability, performance, and prepare for WGS SNV caller integration.  
**Execution Validation:** All unit tests executed via `uv run` (66 script tests + 41 input preparation tests passing).

---

## 1. Executive Summary

| Metric | Assessment |
|---|---|
| **Overall Status** | **Ready with minor fixes** |
| **Architecture Quality** | **Excellent** (Clean separation of `library_type` and `assay`; standalone testable routing module) |
| **Correctness & Safety** | **High** (Dynamic platform mapping, `when:` validation, reference checksums) |
| **Top Blocker for WGS** | `clairs_to_enabled()` currently hard-fails at DAG parse time if `reference.design_bed` is empty, blocking whole-genome WGS variant calling. |

---

## 2. Source Inventory

| File | Lines | Layer | Responsibility | Status |
|---|---:|---|---|---|
| `workflow/scripts/nagini_routing.py` | 595 | Routing Core | Independent registry for assays, library types, minimap2 presets, ClairS-TO platform resolution, barcode checks | **Production** |
| `workflow/scripts/test_nagini_routing.py` | 371 | Unit Tests | Tests for assay inheritance, override precedence, barcode validation, condition parsing | **Production** |
| `workflow/scripts/nagini_report.py` | 653 | Reporting | Excel workbook compilation with `xlsxwriter` (Overview, SNVs, SVs, Translocations, QC, BCR-ABL1) | **Production** |
| `workflow/scripts/test_nagini_report.py` | 321 | Unit Tests | XML/content-level assertion tests on generated Excel workbooks | **Production** |
| `workflow/scripts/dorado_demux_summary.py` | 169 | QC / Dorado | Extraction and aggregation of per-read demultiplexing stats into MultiQC format | **Production** |
| `pipeline_start_scripts/prepare_nagini_inputs.py` | 378 | Entrypoint | Site config resolution, reference registry verification, barcode validation, `samples.tsv`/`units.tsv` generation | **Production** |
| `workflow/rules/common.smk` | 845 | Snakemake Core | Global DAG setup, dynamic input helpers, output spec evaluation | **Production** |

---

## 3. Review of Recent Changes

### 3.1. Sample Routing & Reference Registry (`17c14ef`)
* **Strengths:**
  - **Decoupled Orthogonal Properties:** Separating `library_type` (physical molecule type $\rightarrow$ alignment preset: `splice` vs `map-ont`) from `assay` (clinical analysis pipeline: `wgs`, `amplicon`, `cdna`, `bcr_abl`, `adaptive`) resolves the longstanding architectural bug where cDNA amplicons could not be represented.
  - **Zero Snakemake Dependency in Core Routing:** `nagini_routing.py` imports no Snakemake globals, allowing instant deterministic unit tests.
  - **Safe `when:` Condition Parser:** Validates `when:` tags at parse time to prevent typo-induced silent dropping of outputs.
  - **Reference Registry:** `config/reference_files/*.yaml` centralizes genome paths and separates site runtime overrides from static genome resources.

### 3.2. Excel Report Migration to `xlsxwriter` (`351694b`)
* **Strengths:**
  - Standardized on `xlsxwriter`, aligning with other SciLifeLab/CGU pipelines.
  - Clean style palette (0–10) with explicit color codes.
  - Distinction between `LAYOUT` sheets (Overview, BCR-ABL1) and `TABLE` sheets (autofilter + freeze panes).
  - Sanitized sheet naming (max 31 characters, illegal character escaping).

### 3.3. ClairS-TO Model Derivation & Demux QC (`dd9ab8e`)
* **Strengths:**
  - **Dynamic Error Model Selection:** `get_clairs_to_platform()` inspects Dorado model string (`dna_r10.4.1_e8.2_400bps_sup@v5.2.0` $\rightarrow$ `ont_r10_dorado_sup_5khz`) preventing silent error model mismatches.
  - **FASTA alignment fallback for spliced cDNA:** Prevents running spliced reads against a pre-indexed `.mmi` built with map-ont window parameters.
  - **Pre-flight Barcode Kit Checks:** Catches invalid barcode indices before dispatching expensive GPU basecalling jobs.

---

## 4. Key Findings & Recommended Fixes

### [High] ClairS-TO Hard-Dependency on `design_bed` Blocks WGS SNV Calling
* **File:** `workflow/rules/common.smk:619`
* **Evidence:**
  ```python
  if clairs_to_enabled() and not has_design_bed():
      _needs_calling = [
          str(sample).strip()
          for sample in get_samples(samples)
          if sample_needs_dna_variant_calling(str(sample).strip())
      ]
      if _needs_calling:
          raise WorkflowError(
              "ClairS-TO needs reference.design_bed, which is not set or not readable..."
          )
  ```
* **Impact:** Any WGS sample (`target_restricted: false`) routed to small variant calling crashes the pipeline unless a target BED is configured.
* **Suggested Fix:**
  1. For `target_restricted: false` samples, allow ClairS-TO to run whole-genome (without `--bed_fn` or using a whole-genome callable regions BED / per-chromosome parallelization).
  2. Or route WGS germline samples to Hydra's `deepvariant` module and WGS somatic samples to whole-genome `clairs_to` / `deepsomatic`.

### [Medium] Dynamic Spliced Minimap2 Index vs Reference Registry
* **File:** `workflow/rules/minimap2_fastq_align.smk`, `workflow/rules/align_dorado_bam.smk`
* **Observation:** Spliced libraries fall back to `{input.ref_fasta}` instead of `.mmi`. While correct, generating a dedicated spliced `.mmi` index for large reference genomes saves ~2-3 minutes of indexing overhead per sample.

---

## 5. Plan: Integrating SNV Callers for WGS Mode using Hydra

Nagini can leverage existing Hydra-genetics modules in `vendor/hydra-genetics/snv_indels`:

1. **WGS Germline Calling (DeepVariant)**:
   - Hydra rule: `deepvariant` from `vendor/hydra-genetics/snv_indels/workflow/rules/deepvariant.smk`.
   - Sharded per chromosome (`{chr}`) using `threads: 8-16` with `--model_type=ONT_R104` or `WGS`.
   - Concatenated with `bcftools concat` and indexed.

2. **WGS Somatic Tumor-Only Calling (ClairS-TO / DeepSomatic)**:
   - **ClairS-TO Whole-Genome**: Run sharded by chromosome or across whole genome without `--bed_fn`.
   - **DeepSomatic Tumor-Only**: Hydra rule `deepsomatic_t_only` from `vendor/hydra-genetics/snv_indels/workflow/rules/deepsomatic.smk`.

3. **Routing & Dispatch in Nagini**:
   - In `config/config.yaml`, configure `small_variant_caller: "clairs_to"` (somatic/amplicon) or `"deepvariant"` (germline WGS) or `"auto"` (derived from sample type: Tumor $\rightarrow$ ClairS-TO, Normal $\rightarrow$ DeepVariant).
