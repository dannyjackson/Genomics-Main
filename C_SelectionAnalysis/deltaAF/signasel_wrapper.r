#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

# Command-line args: species_name directory_with_mafs_files
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript lrt_delta_af.R <species> <directory>")

species <- args[1]
dir <- args[2]

# Paths to input files
pre_file <- paste0(dir, "/", species, "_pre.mafs.gz")
post_file <- paste0(dir, "/", species, "_post.mafs.gz")

# Function to read .mafs.gz ANGSD file
read_maf_file <- function(path) {
  fread(cmd = paste("zcat", path))
}

# Function to compute allele count from frequency and sample size
compute_allele_counts <- function(freq, nInd) {
  round(freq * 2 * nInd)
}

# Selection test function for one SNP
run_signasel_on_mafs <- function(pre_file, post_file, Ne = 1000, generations = 5, max_s = 1) {
  cat("Reading files...\n")
  maf_pre <- read_maf_file(pre_file)
  maf_post <- read_maf_file(post_file)

  cat("Merging on SNP coordinates...\n")
  maf_merged <- merge(maf_pre, maf_post, by = c("chromo", "position"), suffixes = c(".pre", ".post"))
  if (nrow(maf_merged) == 0) stop("No overlapping SNPs found between pre and post files.")

  cat("Running selection tests...\n")
  results_list <- vector("list", nrow(maf_merged))

  for (i in seq_len(nrow(maf_merged))) {
    i_pre <- compute_allele_counts(maf_merged$knownEM.pre[i], maf_merged$nInd.pre[i])
    S_pre <- maf_merged$nInd.pre[i]

    i_post <- compute_allele_counts(maf_merged$knownEM.post[i], maf_merged$nInd.post[i])
    S_post <- maf_merged$nInd.post[i]

    data_mat <- matrix(c(
      0, i_pre, S_pre, Ne,
      generations, i_post, S_post, Ne
    ), ncol = 4, byrow = TRUE)

    res <- tryCatch({
      signaseltest(data_mat, maxs = max_s)
    }, error = function(e) {
      matrix(c(NA, NA, NA, NA, NA, 1), nrow = 1)
    })

    results_list[[i]] <- cbind(maf_merged[i, .(chromo, position)], as.data.frame(res))
  }

  rbindlist(results_list)
}

# ---- signasel functions (abbreviated for clarity here, you already included full defs above) ----
################################################################
WFVecProbS <- function(p, N, s) {
################################################################
    ## Gives the vector of probabilities of all possible allele
    ## numbers (0 to 2N) in a population of size N, when the frequency
    ## at previous generation was p and selection coefficient is s,
    ## according to Wrigt-Fisher model. This vector will serve as a
    ## column of the recursion matrix, so p(i,j) is element a(j,i) of
    ## the matrix.
    ## p: allele frequency at generations t.
    ## N: population sizes at generations (t+1)
    ## s: coefficient of selection
################################################################
    ## Fitnesses of the 3 genotypes A1A1, A1A2, A2A2
    ## (A1 is the selected allele)
    w11 <- 1 + 2 * s
    w12 <- 1 + s
    w22 <- 1
    a <- w11 - w12
    b <- w12 - w22
    c <- a - b
    ## Allele frequencies after selection at generations t.
    pprime <- 1. * (p * (p * a + w12)) / (2. * p * b + w22 + p * p * c)
    ## Probabilities at generations (t+1)
    prob <- dbinom(0:(2*N), 2*N, pprime)
    return(prob)
}
################################################################
WFMatrixS <- function(S1, S2, N, ng, s){
################################################################
  ## Create the Wright-Fisher recursion matrix over ng generations
  ## s: coefficient of selection
  ## N: population size (Ne), for generations 1 to ng-1
  ## ng: number of generations (iterations)
  ## S1: sample size in first generation
  ## S2: sample size in last generation.
################################################################
  ## The recursion matrix for generation 1 (pop size S1->N)
  ## NB: initial pop size should be N, not S1, but we are only using
  ## allele frequency p=k/2S1, so estimates remain the same regardless
  ## of the true pop size.
  if (ng > 1) {
    matA <- NULL
    for (k in 0:(2*S1)) {                  
      matA <- cbind(matA, WFVecProbS(k/(2*S1), N, s))
    }
  }
  ## The recursion matrix for generations 2 to (ng-1) (pop size N->N)
  if (ng > 2) {
    matA2 <- NULL
    for (k in 0:(2*N)) {
      matA2 <- cbind(matA2, WFVecProbS(k/(2*N), N, s))
    }
    for (g in 2:(ng-1)) {
      matA <- matA2 %*% matA
    }
  }
  ## The recursion matrix for generation g
  ## (pop size N->S2 if ng>1 or S1->S2 if ng=1)
  matB <- NULL
  if (ng > 1) {
    for (k in 0:(2*N)) {
      matB <- cbind(matB, WFVecProbS(k/(2*N), S2, s))
    }
    matrec <- matB %*% matA
  } else {
    for (k in 0:(2*S1)) {
      matB <- cbind(matB, WFVecProbS(k/(2*S1), S2, s))
    }
    matrec <- matB
  }
  return(matrec)
}
################################################################
WFLike1S <- function(i1, S1, i2, S2, N, ng, s) {
################################################################
  ## Calculates the Likelihood, ie simply multiplies the recursion
  ## matrix by the vector of initial frequencies
################################################################
  ## 'negative' selection is done by symmetry, ie positive selection
  ## for the other allele
  if (s<0) {
    s <- -s
    i1 <- 2 * S1 - i1
    i2 <- 2 * S2 - i2
  }
  stopifnot(s>=0)
  mat <- WFMatrixS(S1, S2, N, ng, s)
  v1 <- rep(0, 2*S1+1)
  v1[i1+1] <- 1  ## i+1 here because i is in 0:(2n) but the first
                 ## element of a vector is 1 in R
  v2 <- mat %*% v1
  stopifnot(length(v2) == 2*S2+1)
  p2 <- v2[i2+1]
  return(p2)
}
################################################################
  
################################################################
WFLike2S <- function(s, data) {
################################################################
  ## Calculates the Likelihood for multiple time samples
  ## Data is assumed to have the format of a matrix with nrow=number
  ## of samples and for each sample the columns: g, i, S, N:
  ## g0 i0 S0 N0
  ## g1 i1 S1 N1
  ## g2 i2 S2 N2
  ## etc.
  ## with g: generation, i: number of allele copies, S: sample size,
  ## N: (effective) population size
  ## We simply call WFLike1S from g1 to g2 (ng=g2-g1), then from
  ## g2 to g3, etc. and multiply the resulting probabilities. s is the
  ## same for all samples and is maximized over all generations.
################################################################
  g <- data[,1]
  i <- data[,2]
  S <- data[,3]
  N <- data[,4]
  p2 <- 1
  for (k in 1:(length(g)-1)) {
      p2 <- p2 *
          WFLike1S(i[k], S[k], i[k+1], S[k+1], N[k+1],g[k+1]-g[k], s)
  }
  return(p2)
}
################################################################

################################################################
WFMaxiLike2S <- function(data, maxs) {
################################################################
  ## Finds the maximum likelihood on s, returns L(smax), smax
  ## Uses the 'optimize' function of R
  ##
  ## IMPORTANT: maxs gives min and max possible values for s 
  ## must be realistic (=1?)  
################################################################
  res <- optimize(WFLike2S, interval=c(-maxs, maxs), 
                  data,
                  maximum = TRUE)
  warn <- 0
  ## Issue warning message if maxs is reached
  if ((res$maximum - maxs)**2 < 1e-6) {
    print(paste("Warning: Maximum value smax = ", maxs,
                  " was reached."))
    warn <- 1
  }
  return(c(res$objective, res$maximum, warn))
}
################################################################
checkparam <- function(data) {
################################################################
  ## Data is assumed to have the format of a matrix with nrow=number
  ## of samples and for each sample the columns: g, i, S, N:
  ## g0 i0 S0 N0
  ## g1 i1 S1 N1
  ## g2 i2 S2 N2
  ## etc.
  ## with g: generation, i: number of allele copies, S: sample size,
  ## N: (effective) population size
################################################################
  g <- data[,1]
  i <- data[,2]
  S <- data[,3]
  n <- data[,4]
  stopifnot(all(data[,c(3,4)]>0))
  stopifnot(all(data[,c(1,2)]>=0))
  stopifnot(all(S<=n))
  stopifnot(all(i<=(2*S)))
}
################################################################
signaseltest <- function(data, maxs = 1) {
################################################################
  ## Compute the test statistics to detect selection
  ## Data is a matrix with nrow = number of samples
  ## For each sample (row) the columns are: g, i, S, N:
  ## g0 i0 S0 N0
  ## g1 i1 S1 N1
  ## g2 i2 S2 N2
  ## etc.
  ## with g: generation, i: number of allele copies, S: sample size,
  ## N: (effective) population size
################################################################
  ## verify parameter values
  checkparam(data)
  ## likelihood of the null (s=0)
  L0 <- WFLike2S(s=0, data)
  ## maximum likelihood for the alternative
  x <- WFMaxiLike2S(data, maxs)
  Lmax <- x[1]
  smax <- x[2]
  warn <- x[3]
  ## likelihood ratio test
  LRT <- -2 * log(L0 / Lmax)
  ## pvalue assuming LRT follows Chi-square with 1 df
  pvalue <- -log10(1 - pchisq(LRT, 1))
  res <- matrix(c(L0, Lmax, smax, LRT, pvalue, warn), nrow=1)
  colnames(res) <- c('L0', 'Lmax', 'smax', 'LRT', '-log10pvalue', 'warn')
  rownames(res) <- ""
  return(res)
}

################################################################



# Run the analysis
results <- run_signasel_on_mafs(pre_file, post_file)

# Save results
fwrite(results, paste0(dir, "/", species, "_signasel_selection_results.tsv"), sep = "\t")
cat("Results saved.\n")
