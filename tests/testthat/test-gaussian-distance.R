library(testthat)
library(geneSetEmbedding)

test_that("Gaussian distances match analytical references", {
  mu <- matrix(c(0, 2, 1, -1), nrow = 2, byrow = TRUE)
  v <- matrix(c(1, 4, 9, 1), nrow = 2, byrow = TRUE)
  mu2 <- matrix(c(1, -2, 3, 0, -1, 1), nrow = 3, byrow = TRUE)
  v2 <- matrix(c(4, 1, 1, 9, 2, 3), nrow = 3, byrow = TRUE)

  rownames(mu) <- rownames(v) <- c("A", "B")
  rownames(mu2) <- rownames(v2) <- c("X", "Y", "Z")

  # Force the pure-R path without changing the package namespace.
  fallback <- gsemb_set_gaussian_distance
  environment(fallback) <- list2env(
    list(.native_routine_available = function(symbol) FALSE),
    parent = environment(fallback)
  )

  for (metric in c("w2", "sym_kl")) {
    expected <- matrix(
      0, nrow = 2, ncol = 3,
      dimnames = list(rownames(mu), rownames(mu2))
    )

    for (i in seq_len(2)) {
      for (j in seq_len(3)) {
        delta2 <- (mu[i, ] - mu2[j, ])^2
        expected[i, j] <- if (metric == "w2") {
          sum(delta2 + (sqrt(v[i, ]) - sqrt(v2[j, ]))^2)
        } else {
          0.5 * sum(
            v[i, ] / v2[j, ] + v2[j, ] / v[i, ] +
              delta2 * (1 / v[i, ] + 1 / v2[j, ]) - 2
          )
        }
      }
    }

    symbol <- if (metric == "w2") {
      "_geneSetEmbedding_w2_distance"
    } else {
      "_geneSetEmbedding_sym_kl_distance"
    }
    expect_true(is.loaded(symbol, PACKAGE = "geneSetEmbedding"))

    expect_equal(
      fallback(mu, v, mu2, v2, metric = metric),
      expected, tolerance = 1e-10
    )
    expect_equal(
      gsemb_set_gaussian_distance(mu, v, mu2, v2, metric = metric),
      expected, tolerance = 1e-10
    )
    expect_equal(
      unname(diag(fallback(mu, v, metric = metric))),
      c(0, 0), tolerance = 1e-10
    )
    expect_equal(
      unname(diag(gsemb_set_gaussian_distance(mu, v, metric = metric))),
      c(0, 0), tolerance = 1e-10
    )
  }
})
