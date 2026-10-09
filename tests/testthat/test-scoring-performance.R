library(testthat)

test_that("Gaussian gene-to-set scores match hand-calculated distances", {
  E <- matrix(c(0, 1, 1, -1, 2, 0.5), 3, byrow = TRUE,
              dimnames = list(c("A", "B", "C"), c("V1", "V2")))
  mu <- matrix(c(0, 0, 2, 1), 2, byrow = TRUE,
               dimnames = list(c("S1", "S2"), c("V1", "V2")))
  v <- matrix(c(1, 4, 9, 2), 2, byrow = TRUE, dimnames = dimnames(mu))
  D <- matrix(c(1/4, 4/9, 5/4, 19/9, 65/16, 1/8), 3, byrow = TRUE,
              dimnames = list(rownames(E), rownames(mu)))
  expect_equal(gsemb_gene_to_set_score(E, mu, v, score = "neg_mahalanobis"),
               -D, tolerance = 1e-10)
  want <- -0.5 * sweep(D, 2, log(c(4, 18)), "+")
  expect_equal(gsemb_gene_to_set_score(E, mu, v, score = "loglik"),
               want, tolerance = 1e-10)
})

test_that("Unrestricted softmax selects the expected genes for multiple sets", {
  E <- matrix(c(0, 1, 2), ncol = 1, dimnames = list(c("A","B","C"), "V1"))
  mu <- matrix(c(0, 2), ncol = 1, dimnames = list(c("S1","S2"), "V1"))
  v <- matrix(1, 2, 1, dimnames = dimnames(mu))
  got <- gsemb_make_concise_gene_sets(
    E, mu, v, select = "softmax_mass", score = "neg_mahalanobis",
    restrict_to_members = FALSE, mass = 0.8, min_size = 1
  )
  expect_equal(got, list(S1 = c("A","B"), S2 = c("C","B")))
  restricted <- gsemb_make_concise_gene_sets(
    E, mu, v, gene_sets = list(S1 = rownames(E), S2 = rownames(E)),
    select = "softmax_mass", score = "neg_mahalanobis",
    restrict_to_members = TRUE, mass = 0.8, min_size = 1
  )
  expect_equal(restricted, list(S1 = c("A","B"), S2 = c("C","B")))
})

test_that("Unrestricted softmax handles a single gene", {
  E <- matrix(0, 1, 1, dimnames = list("A", "V1"))
  mu <- matrix(c(0, 2), 2, 1, dimnames = list(c("S1","S2"), "V1"))
  v <- matrix(1, 2, 1, dimnames = dimnames(mu))
  got <- gsemb_make_concise_gene_sets(
    E, mu, v, select = "softmax_mass", restrict_to_members = FALSE,
    min_size = 1
  )
  expect_equal(got, list(S1 = "A", S2 = "A"))
  one_set <- gsemb_make_concise_gene_sets(
    E, mu[1,,drop = FALSE], v[1,,drop = FALSE],
    select = "softmax_mass", restrict_to_members = FALSE, min_size = 1
  )
  expect_equal(one_set, list(S1 = "A"))
})

test_that("Unrestricted softmax caps selection by available genes and max_size", {
  E <- matrix(c(0, 1, 2), 3, 1, dimnames = list(c("A","B","C"), "V1"))
  mu <- matrix(0, 1, 1, dimnames = list("S1", "V1"))
  v <- matrix(1, 1, 1, dimnames = dimnames(mu))
  got <- gsemb_make_concise_gene_sets(
    E, mu, v, select = "softmax_mass", restrict_to_members = FALSE,
    min_size = 5
  )
  expect_equal(got, list(S1 = c("A","B","C")))
  bounded <- gsemb_make_concise_gene_sets(
    E, mu, v, select = "softmax_mass", restrict_to_members = FALSE,
    min_size = 5, max_size = 1
  )
  expect_equal(bounded, list(S1 = "A"))
})

test_that("Gaussian fitting preserves missing variance and sample variance", {
  X <- matrix(c(1, NA, 3, 2, 4, 8), 3, 2,
              dimnames = list(c("A","B","C"), c("V1","V2")))
  fit <- gsemb_fit_set_gaussians_from_members(X, list(S1 = rownames(X)))
  expect_true(is.na(fit$var[1, 1]))
  expect_equal(unname(fit$var[1, 2]), 28/3, tolerance = 1e-10)

  X[2, 1] <- 2
  complete <- gsemb_fit_set_gaussians_from_members(
    X, list(S1 = rownames(X), S2 = "A", S3 = "missing")
  )
  expect_equal(unname(complete$var[1, ]), c(1, 28/3), tolerance = 1e-10)
  expect_equal(unname(complete$var[2, ]), c(1e-8, 1e-8))
  expect_true(all(is.na(complete$var[3, ])))

  # Exercise the real stdlib variance calculation without installing/removing packages.
  fallback_variance <- get(".gsemb_col_vars", envir = environment(gsemb_fit_set_gaussians_from_members))
  environment(fallback_variance) <- list2env(
    list(requireNamespace = function(...) FALSE),
    parent = environment(fallback_variance)
  )
  fallback_fit <- gsemb_fit_set_gaussians_from_members
  environment(fallback_fit) <- list2env(
    list(.gsemb_col_vars = fallback_variance),
    parent = environment(fallback_fit)
  )
  expect_equal(
    fallback_fit(X, list(S1 = rownames(X), S2 = "A", S3 = "missing")),
    complete, tolerance = 1e-10
  )
})
