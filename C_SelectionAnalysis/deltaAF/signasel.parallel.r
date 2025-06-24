#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(future.apply)
})

# Set parallel backend
plan(multisession, workers = availableCores())

# --------------------- Configuration ---------------------
pre_file <- "cra_pre.test.mafs.gz"
post_file <- "cra_post.test.mafs.gz"
Ne <- 1000  # use 1000 for development, 11306700 for testing
generations <- 5
max_s <- 1
out_file <- "signasel_results.tsv"
# ---------------------------------------------------------

# Read gzipped MAF file
read_maf_file <- function(path) {
  fread(cmd = paste("zcat", path))
}

# Compute allele counts from freq and sample size
compute_allele_counts <- function(freq, nInd) {
  round(freq * 2 * nInd)
}

# ---- Wright-Fisher Functions ----
WFVecProbS <- function(p, N, s) {
  w11 <- 1 + 2 * s
  w12 <- 1 + s
  w22 <- 1
  a <- w11 - w12
  b <- w12 - w22
  c <- a - b
  pprime <- (p * (p * a + w12)) / (2 * p * b + w22 + p^2 * c)
  dbinom(0:(2*N), 2*N, pprime)
}

WFMatrixS <- function(S1, S2, N, ng, s) {
  if (ng > 1) {
    matA <- sapply(0:(2*S1), function(k) WFVecProbS(k / (2*S1), N, s))
    if (ng > 2) {
      matA2 <- sapply(0:(2*N), function(k) WFVecProbS(k / (2*N), N, s))
      for (g in 2:(ng - 1)) {
        matA <- matA2 %*% matA
      }
    }
    matB <- sapply(0:(2*N), function(k) WFVecProbS(k / (2*N), S2, s))
    return(matB %*% matA)
  } else {
    matB <- sapply(0:(2*S1), function(k) WFVecProbS(k / (2*S1), S2, s))
    return(matB)
  }
}

WFLike1S <- function(i1, S1, i2, S2, N, ng, s) {
  if (s < 0) {
    s <- -s
    i1 <- 2 * S1 - i1
    i2 <- 2 * S2 - i2
  }
  mat <- WFMatrixS(S1, S2, N, ng, s)
  v1 <- numeric(2 * S1 + 1)
  v1[i1 + 1] <- 1
  v2 <- mat %*% v1
  v2[i2 + 1]
}

WFLike2S <- function(s, data) {
  g <- data[, 1]
  i <- data[, 2]
  S <- data[, 3]
  N <- data[, 4]
  p <- 1
  for (k in 1:(length(g) - 1)) {
    p <- p * WFLike1S(i[k], S[k], i[k + 1], S[k + 1], N[k + 1], g[k + 1] - g[k], s)
  }
  return(p)
}

WFMaxiLike2S <- function(data, maxs) {
  res <- optimize(WFLike2S, interval = c(-maxs, maxs), data, maximum = TRUE)
  warn <- ifelse((res$maximum - maxs)^2 < 1e-6, 1, 0)
  return(c(res$objective, res$maximum, warn))
}

checkparam <- function(data) {
  stopifnot(all(data[, c(3,4)] > 0), all(data[, c(1,2)] >= 0), all(data[,3] <= data[,4]), all(data[,2] <= 2 * data[,3]))
}

signaseltest <- function(data, maxs = 1) {
  checkparam(data)
  L0 <- WFLike2S(0, data)
  x <- WFMaxiLike2S(data, maxs)
  Lmax <- x[1]; smax <- x[2]; warn <- x[3]
  LRT <- -2 * log(L0 / Lmax)
  pval <- -log10(1 - pchisq(LRT, 1))
  matrix(c(L0, Lmax, smax, LRT, pval, warn), nrow = 1,
         dimnames = list(NULL, c("L0", "Lmax", "smax", "LRT", "-log10pvalue", "warn")))
}

# --------------------- Main Parallel SNP Test ---------------------

run_signasel_on_mafs <- function(pre_file, post_file, Ne, generations, max_s) {
  cat("Reading files...\n")
  maf_pre <- read_maf_file(pre_file)
  maf_post <- read_maf_file(post_file)

  cat("Merging on SNPs...\n")
  maf <- merge(maf_pre, maf_post, by = c("chromo", "position"), suffixes = c(".pre", ".post"))
  if (nrow(maf) == 0) stop("No overlapping SNPs found.")

  cat("Running parallel selection tests on", nrow(maf), "SNPs...\n")

  results_list <- future_lapply(seq_len(nrow(maf)), function(i) {
    i_pre <- compute_allele_counts(maf$knownEM.pre[i], maf$nInd.pre[i])
    i_post <- compute_allele_counts(maf$knownEM.post[i], maf$nInd.post[i])

    data_mat <- matrix(c(
      0, i_pre, maf$nInd.pre[i], Ne,
      generations, i_post, maf$nInd.post[i], Ne
    ), ncol = 4, byrow = TRUE)

    res <- tryCatch(
      signaseltest(data_mat, max_s),
      error = function(e) matrix(c(NA, NA, NA, NA, NA, 1), nrow = 1,
                                 dimnames = list(NULL, c("L0", "Lmax", "smax", "LRT", "-log10pvalue", "warn")))
    )

    cbind(maf[i, .(chromo, position)], as.data.frame(res))
  }, future.seed = TRUE)

  rbindlist(results_list)
}

# --------------------- Execute and Save ---------------------

results <- run_signasel_on_mafs(pre_file, post_file, Ne, generations, max_s)
fwrite(results, file = out_file, sep = "\t")
cat("Results saved to", out_file, "\n")
