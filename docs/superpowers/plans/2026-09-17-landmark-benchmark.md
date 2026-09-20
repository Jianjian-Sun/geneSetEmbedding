# Landmark Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reproducible, extensible Landmark benchmark for the four requested degree-based strategies and nine requested metrics.

**Architecture:** A sourced function library owns validation, method registration, metrics, perturbation, and output writing. A small configuration file supplies dataset-specific values, and a single command-line script runs the benchmark.

**Tech Stack:** R 4.1+, Matrix, base R, testthat for the validation script.

## Global Constraints

- Preserve a common node universe across all and high-confidence graphs.
- Use canonical unique undirected edges before any metric or perturbation.
- Use one shared dropped-edge set per replicate across all methods.
- Allow registered methods to return `score = NULL`.
- Avoid new package dependencies.

---

### Task 1: Core edge and method interfaces

**Files:**
- Create: `scripts/landmark_benchmark.R`
- Test: `scripts/tests/test_landmark_benchmark.R`

**Interfaces:**
- Consumes: edge data frame and validated configuration.
- Produces: `prepare_landmark_context()`, `default_landmark_methods()`, and `run_landmark_methods()`.

- [ ] Write tests for canonical edges, common nodes, deterministic ties, and the four default methods.
- [ ] Run the tests and confirm they fail because the functions do not exist.
- [ ] Implement the minimum edge preparation and method registry.
- [ ] Run the tests and confirm the new tests pass.

### Task 2: Metrics and edge perturbation

**Files:**
- Modify: `scripts/landmark_benchmark.R`
- Modify: `scripts/tests/test_landmark_benchmark.R`

**Interfaces:**
- Consumes: prepared context and registered method results.
- Produces: pairwise Jaccard/Spearman, node quality summaries, coverage summaries, and dropout replicate/summary tables.

- [ ] Add failing toy-network tests for all requested metrics and a scoreless custom method.
- [ ] Implement aligned score comparisons, sparse multi-source BFS, and shared multi-seed edge dropout.
- [ ] Run the focused validation script.

### Task 3: Reproducible entry point and documentation

**Files:**
- Create: `scripts/run_landmark_benchmark.R`
- Create: `scripts/example_landmark_config.R`
- Create: `scripts/README_landmark_benchmark.md`

**Interfaces:**
- Consumes: one R configuration file.
- Produces: CSV/RDS/session artifacts in the configured output directory.

- [ ] Add a temporary toy configuration and verify command-line execution.
- [ ] Implement configuration loading, artifact writing, and provenance capture.
- [ ] Run the focused tests and one end-to-end toy benchmark.
