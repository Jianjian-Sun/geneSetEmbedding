# Landmark Benchmark Design

## Goal

Provide a reproducible, method-extensible benchmark for Landmark selection on weighted PPI networks. The first release compares degree and weighted degree on all edges and high-confidence edges, while leaving a stable registration interface for later Betweenness and K-Medoids selectors.

## Inputs

The benchmark accepts a CSV, TSV, or RDS edge table through an R configuration file. The configuration identifies endpoint columns, the confidence column, threshold rule, duplicate-edge rule, Landmark count, random-edge-drop settings, evaluation graph, seed, and output directory.

Undirected edges are canonicalized, self-loops and invalid rows are rejected or counted, and duplicate undirected pairs are collapsed using an explicit rule. All selection graphs retain one common node universe, including nodes isolated by high-confidence filtering.

## Method interface

Methods are a named list of functions with signature `function(context, k, seed)`. Each returns:

```r
list(
  landmarks = character(),
  score = named_numeric_vector_or_NULL,
  score_name = character_scalar,
  selection_graph = "all" | "high"
)
```

The initial registry contains `degree_all`, `weighted_degree_all`, `degree_high`, and `weighted_degree_high`. A future selector without a meaningful global score may return `score = NULL`; score Spearman values involving that method are then recorded as `NA` with a reason.

## Metrics

The benchmark computes:

1. Pairwise Landmark Jaccard for every method pair.
2. Pairwise Spearman correlation of aligned all-node method scores.
3. Median Landmark total degree.
4. Median Landmark high-confidence degree.
5. Median per-Landmark mean incident-edge confidence.
6. Median per-Landmark low-confidence-edge fraction.
7. Coverage at two hops on one configured reference graph.
8. Landmark Jaccard stability after repeatedly deleting 5% of canonical edges using shared perturbations across methods.
9. Mean nearest-Landmark distance among reachable nodes, accompanied by the unreachable-node fraction.

Node-level, method-level, and replicate-level results are retained so summaries can be audited.

## Reproducibility

One master seed deterministically generates replicate seeds. Each replicate records the selected edge IDs, actual edge counts, and high-confidence deletion fraction. Outputs also include the complete configuration, content fingerprints with their recorded hash algorithm, R session information, and input quality-control counts.

## Implementation constraints

Use base R and Matrix only. Coverage and nearest distances use sparse multi-source breadth-first search without constructing a node-by-Landmark distance matrix. The command-line entry point accepts exactly one configuration file and writes a fixed collection of CSV/RDS/text artifacts.

## Validation

A deterministic toy network verifies edge canonicalization, method selection, Jaccard, Spearman alignment, node-quality metrics, Coverage@2, unreachable distances, shared edge dropout, and registration of a custom scoreless selector.
