# SLE Landmark Real-Cohort Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a three-script R workflow that compares three fixed Landmark methods across three existing SLE cohorts using Weighted-GSEA and Concise-ORA without duplicating existing algorithms.

**Architecture:** A parameter-only entry script sources one thin helper module and the repository's existing modular ORA/comparison scripts. A separate commented walkthrough explains the same data flow without running expensive work by default. Unit tests use toy matrices and tables so correctness can be checked without the real STRING graph.

**Tech Stack:** R 4.1+, geneSetEmbedding package code, Matrix, testthat, existing base-R modular ORA scripts.

## Global Constraints

- Work only in `E:/genesetembedding/geneSetEmbedding/.worktrees/weighted-degree-landmarks`.
- Reuse existing RWR, SVD, Gaussian, enrichment, ORA, Jaccard and Spearman functions.
- Add no package dependency and do not modify the public `gsemb_fit()` API.
- Keep Weighted-GSEA as the primary analysis and Concise-ORA as a sensitivity analysis.
- Do not rerun Landmark selection when existing RDS outputs are present.

---

### Task 1: Test the missing glue interfaces

**Files:**
- Create: `tests/testthat/test-sle-landmark-real-cohort.R`
- Create: `scripts/real_cohort/02_sle_landmark_helpers.R`

**Interfaces:**
- Produces: `prepare_ranked_stats()`, `build_fit_from_landmarks()`, `as_comparable_enrichment()`, `summarize_method_reproducibility()`.

- [ ] Write tests asserting duplicate-gene resolution, fixed-Landmark fit construction, enrichment-column normalization and method-level averaging.
- [ ] Run `testthat::test_file()` and verify RED because helper functions do not exist.
- [ ] Implement the minimum helper functions using existing package functions.
- [ ] Run the targeted test and verify GREEN.

### Task 2: Add the parameter-only main workflow

**Files:**
- Create: `scripts/real_cohort/01_run_sle_landmark_benchmark.R`

**Interfaces:**
- Consumes: helper functions from Task 1 and existing `step3_modular` functions.
- Produces: cached fits, per-cohort Weighted-GSEA and Concise-ORA CSVs, cross-cohort comparison CSVs, combined method summary and run metadata.

- [ ] Declare all editable parameters at the top of the script.
- [ ] Validate required package/data/Landmark files before expensive work.
- [ ] Load exactly the three selected high-network Landmark sets.
- [ ] Build or reuse three cached fits.
- [ ] Run or reuse Weighted-GSEA results and pass them to existing comparison functions.
- [ ] Generate or reuse method-specific Concise sets, run existing ORA functions, and pass results to the same comparison functions.
- [ ] Write one combined summary and metadata record.
- [ ] Parse the script without executing it and verify no syntax errors.

### Task 3: Add the beginner walkthrough

**Files:**
- Create: `scripts/real_cohort/03_sle_landmark_beginner_walkthrough.R`

**Interfaces:**
- Consumes: the main workflow path.
- Produces: a readable, parseable Chinese tutorial with `RUN_ANALYSIS <- FALSE` safety default.

- [ ] Explain Landmark, fit, Weighted-GSEA, Concise-ORA and cross-cohort metrics in execution order.
- [ ] Show which parameters a beginner may safely change and which must stay fixed for fairness.
- [ ] Explain every output directory and the allowed conclusion language.
- [ ] Parse the walkthrough and verify no syntax errors.

### Task 4: Verification

**Files:**
- Verify all files above.

- [ ] Run the targeted helper tests.
- [ ] Parse all three scripts.
- [ ] Run a toy end-to-end check for both comparison routes without the real STRING computation.
- [ ] Inspect the worktree diff and confirm no unrelated files changed.
- [ ] Record that the Windows Git view cannot commit if the existing worktree `.git` file still points to the unavailable WSL gitdir.
