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
# file <- "intersection.noca_all.nocorrection.GO"
# file <- "intersection.pyrr_all.nocorrection.GO"
# file <- "noca.all.genenames_filtered.GO"
cat("reading go list\n")


go_list_na <- read.delim(paste0(file, ".txt"))

go_list <- go_list_na[!is.na(go_list_na) & go_list_na != ""]

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

outfile <- paste0("plot/", file, ".scatter.pdf")
pdf(outfile, width = 10, height =10)  
scatterPlot(simMatrix, reducedTerms, labelSize=10)
dev.off()

cat("printing treemap \n")

outfile <- paste0("plot/", file, "_treemap.pdf")
pdf(outfile, width = 10, height = 10)  
treemapPlot(reducedTerms, size="score")
dev.off()

cat("printing heatmap \n")

outfile <- paste0("plot/", file, ".heatmap.pdf")
pdf(outfile, width = 10, height = 10)  
heatmapPlot(simMatrix,
            annotateParent=TRUE,
            fontsize=24)
dev.off()

