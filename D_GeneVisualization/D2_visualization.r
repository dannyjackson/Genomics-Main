# BiocManager::install("org.Gg.eg.db")
library(org.Gg.eg.db)
library(limma)
library(rrvgo)

# common genes in both

# go_list <- c("GO:XXXXXXX","GO:XXXXXXX")



simMatrix <- calculateSimMatrix(go_list,
                                orgdb="org.Gg.eg.db",
                                ont="BP",
                                method="Rel")

reducedTerms <- reduceSimMatrix(simMatrix,
                                orgdb="org.Gg.eg.db")


scatterPlot(simMatrix, reducedTerms)

treemapPlot(reducedTerms)


heatmapPlot(simMatrix,
            annotateParent=TRUE,
            fontsize=6)
