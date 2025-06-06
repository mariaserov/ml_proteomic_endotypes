install.packages("BiocManager")


BiocManager::install("clusterProfiler")
BiocManager::install("KEGGREST")
BiocManager::install("pathview")
BiocManager::install("fgsea")
BiocManager::install("msigdbr")

library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGGREST)
library(pathview)
library(fgsea)
library(msigdbr)
library(enrichplot)
library(dplyr)
library(readr)
library(rstudioapi)

current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )

# LOAD IDs

clusters = c(1,2,3,4)

for (c in clusters) {
  
  print(c)
  
  # Load cluster
  path <- paste0('/Users/mashaserova/Downloads/pathway_analysis_data/POIs_cluster', c, ".txt")
  cl <- readLines(path)
  cl <- cl[1:length(cl)]
  
  print(cl[1])
  cl <- AnnotationDbi::select(org.Hs.eg.db, keys = cl, keytype = "UNIPROT", columns = c("ENTREZID", "SYMBOL"))
  
  # Pathway analysis
  kegg <- enrichKEGG(gene=cl$ENTREZID, organism='hsa', pvalueCutoff = 0.01)
  
  # Get rid of disease terms
  filtered <- kegg 
  filtered@result <- filtered@result[!grepl("hsa05|hsa06", filtered@result$ID) & filtered@result$p.adjust < 0.05,] 
  
  # Dot plot
  p <- dotplot(filtered, showCategory = 20)
  ggsave(paste0("cluster", c, "_dotplot.png"), plot = p, width = 10, height = 6, dpi = 300)
  
  # Get proteins per pathway
  
  enriched_pathways <- as.data.frame(filtered@result)
  enriched_genes <- enriched_pathways %>% select(ID, Description, geneID) %>% separate_rows(geneID, sep = "/") 
  enriched_genes <- as.data.frame(enriched_genes)
  enriched_genes$Protein <- mapIds(org.Hs.eg.db, keys = enriched_genes$geneID, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")
  
  write.csv(enriched_genes, paste0("enriched_pathways_with_proteins/cluster", c, "_enriched_pathways_prots.csv"),
            row.names = FALSE)
}



cl1 <- readLines("/Users/mashaserova/Downloads/pathway_analysis_data/POIs_cluster1.txt")
cl1 <- cl1[2:length(cl1)]
cl1 <- AnnotationDbi::select(org.Hs.eg.db, keys = cl1, keytype = "UNIPROT", columns = c("ENTREZID", "SYMBOL"))

cl2 <- readLines("/Users/mashaserova/Downloads/pathway_analysis_data/POIs_cluster2.txt")
cl2 <- cl2[2:length(cl2)]
cl2 <- AnnotationDbi::select(org.Hs.eg.db, keys = cl2, keytype = "UNIPROT", columns = c("ENTREZID", "SYMBOL"))

cl3 <- readLines("/Users/mashaserova/Downloads/pathway_analysis_data/POIs_cluster3.txt")
cl3 <- cl3[2:length(cl3)]
cl3 <- AnnotationDbi::select(org.Hs.eg.db, keys = cl3, keytype = "UNIPROT", columns = c("ENTREZID", "SYMBOL"))

cl4 <- readLines("/Users/mashaserova/Downloads/pathway_analysis_data/POIs_cluster4.txt")
cl4 <- cl4[2:length(cl4)]
cl4 <- AnnotationDbi::select(org.Hs.eg.db, keys = cl4, keytype = "UNIPROT", columns = c("ENTREZID", "SYMBOL"))




# KEGG for CL1

cl = cl1

kegg <- enrichKEGG(gene=cl$ENTREZID, organism='hsa', pvalueCutoff = 0.01)

filtered <- kegg

filtered@result <- filtered@result[
  !grepl("hsa05|hsa06", filtered@result$ID) &
    filtered@result$p.adjust < 0.05,
]

dim(filtered)
dotplot(filtered, showCategory = 20)


# KEGG for CL2 - ZERO RESULTS

cl = cl2

kegg <- enrichKEGG(gene=cl$ENTREZID, organism='hsa', pvalueCutoff = 0.01)

filtered <- kegg

filtered@result <- filtered@result[
  !grepl("hsa05|hsa06", filtered@result$ID) &
    filtered@result$p.adjust < 0.05,
]

dim(filtered)
p <- dotplot(filtered, showCategory = 20)
ggsave("dotplot.png", plot = p, width = 10, height = 6, dpi = 300)

# KEGG for CL3

cl = cl3

kegg <- enrichKEGG(gene=cl$ENTREZID, organism='hsa', pvalueCutoff = 0.01)

filtered <- kegg

filtered@result <- filtered@result[
  !grepl("hsa05|hsa06", filtered@result$ID) &
    filtered@result$p.adjust < 0.05,
]

dim(filtered)
dotplot(filtered, showCategory = 20)

enriched_pathways <- as.data.frame(filtered@result)

enriched_genes <- enriched_pathways %>%
  select(ID, Description, geneID) %>%
  separate_rows(geneID, sep = "/") 

enriched_genes <- as.data.frame(enriched_genes)

protein_names <- mapIds(org.Hs.eg.db, keys = enriched_genes$geneID, 
                        column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")

enriched_genes$Protein <- protein_names



write.csv(enriched_genes, "enriched_pathways_with_proteins/cluster4_enriched_pathways_prots.csv", row.names = FALSE)

# KEGG FOR CL4

cl = cl4

kegg <- enrichKEGG(gene=cl$ENTREZID, organism='hsa', pvalueCutoff = 0.01)

filtered <- kegg

filtered@result <- filtered@result[
  !grepl("hsa05|hsa06", filtered@result$ID) &
    filtered@result$p.adjust < 0.05,
]

dim(filtered)
dotplot(filtered, showCategory = 20)



