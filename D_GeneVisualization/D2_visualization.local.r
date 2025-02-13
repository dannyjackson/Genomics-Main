#!/usr/bin/env Rscript

# Load required packages, installing if necessary
# BiocManager::install("org.Gg.eg.db")
# BiocManager::install("rrvgo")

required_packages <- c("org.Gg.eg.db", "rrvgo")
installed_packages <- rownames(installed.packages())

cat("Checking required packages...\n")
for (pkg in required_packages) {
  if (!(pkg %in% installed_packages)) {
    BiocManager::install(pkg)
  }
  library(pkg, character.only = TRUE)
}

cat("Parsing command-line arguments...\n")
# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
file <- args[1]

cat("reading go list\n")

vector_data <- scan(paste0(file, ".txt"), what = "", sep = "\t")

# Remove NA values
go_list <- vector_data[!is.na(vector_data) & vector_data != ""]

# Check the cleaned vector
print(go_list)


# go_list <- c("GO:XXXXXXX","GO:XXXXXXX")


cat("calculating sim matrix \n")

simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

cat("reducing terms \n")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")


cat("printing scatterplot \n")

outfile <- paste0("plot/", file, ".scatter.png")
png(outfile, width = 800, height = 600)  
scatterPlot(simMatrix, reducedTerms)
dev.off()

cat("printing treemap \n")

outfile <- paste0("plot/", file, "_treemap.png")
png(outfile, width = 800, height = 600)  
treemapPlot(reducedTerms)
dev.off()

cat("printing heatmap \n")

outfile <- paste0("plot/", file, ".heatmap.png")
png(outfile, width = 800, height = 600)  
heatmapPlot(simMatrix,
            annotateParent=TRUE,
            fontsize=6)
dev.off()

