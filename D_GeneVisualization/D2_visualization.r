#!/usr/bin/env Rscript

# Load required packages, installing if necessary
BiocManager::install("org.Gg.eg.db")

required_packages <- c("org.Gg.eg.db", "rrvgo")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    install.packages(pkg, repos = "http://cran.us.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

cat("Parsing command-line arguments...\n")
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
scriptdir <- args [2]
file <- args[3]

infile <- file.path(scriptdir, "D_CandidateGeneAnalyses/GOterms", $file, ".txt")


# common genes in both

# go_list <- c("GO:XXXXXXX","GO:XXXXXXX")



simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")


outfile <- file.path(outdir, "analyses/genelists/plots", paste0($file, ".scatter.png"))
png("scatterplot.png", width = 800, height = 600)  
scatterPlot(simMatrix, reducedTerms)
dev.off()

outfile <- file.path(outdir, "analyses/genelists/plots", paste0($file, ".treemap.png"))
treemapPlot(reducedTerms)
dev.off()

outfile <- file.path(outdir, "analyses/genelists/plots", paste0($file, ".heatmap.png"))
heatmapPlot(simMatrix,
            annotateParent=TRUE,
            fontsize=6)
dev.off()

