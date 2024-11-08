# Documentation
# https://www.bioconductor.org/packages/release/bioc/html/RCy3.html
# https://www.bioconductor.org/packages/release/bioc/html/OmnipathR.html

# Load packages ----
suppressPackageStartupMessages({ 
  library(OmnipathR)
  library(dplyr)
  library(igraph)
  library(visNetwork)
  library(RCy3)
  library(fgsea)
  library(vroom)
  library(org.Hs.eg.db)
  library(dplyr)
  library(tidyr)
  library(clusterProfiler)
  library(ggVennDiagram)
  library(ggplot2)
  })

# Download protein interactions ----
#interactions <- import_omnipath_interactions(c("SignaLink3", "SIGNOR", "PhosphoSite"))

#post_translational i.e. physical interactions of proteins, protein-protein interactions (or PPIs)
interactions_PPI <- import_post_translational_interactions( 
  organism = 9606
)

qc_PPI_interactions <- interactions_PPI %>%
  filter(curation_effort > 2) 
qc_PPI_interactions$type <- "ppi"


# Perform quality control ----
qc_interactions <- interactions_PPI %>%
  filter(curation_effort > 2) %>%
  mutate(
    type = case_when(
      is_stimulation == 1 ~ "activation",
      is_inhibition == 1 ~ "inhibition",
      consensus_stimulation == 1 ~ "activation",
      consensus_inhibition == 1 ~ "inhibition",
      TRUE ~ NA_character_
    )
  ) %>%
  drop_na(type)

write.table(qc_interactions, "../analysis/6_NetworkCons/OmnipathQC/omnipathR_qc_interactionsCE2.txt", 
            row.names = F, sep = "\t")
 
write.table(qc_PPI_interactions, "../analysis/6_NetworkCons/OmnipathQC/omnipathR_qc_interactionsPPI.txt", 
            row.names = F, sep = "\t")

# Make a platelet network ----
#the nodes that either have an average expression of 1.0293 (25% percentile, Q1) 
#in the RNAseq dataset or were detected in the 2 proteomics datasets

## load proteomics data ----
proteomicsDataSYM <- read.table("../analysis/6_NetworkCons/Annotations/plateletProteome(ackr3data).txt", header = TRUE, sep = "\t") %>%
  pull(1) %>%
  mapIds(org.Hs.eg.db,
         keys = .,
         keytype="UNIPROT",
         column="SYMBOL", multiVals = "first") %>%
  .[!is.na(.)] %>%
  as_tibble(rownames = "UNIPROT_ID")  %>%
  mutate(value = if_else(value == "CALM1", "CALM2", value)) %>% #somehow maps it wrong, so i correct manually.
  pull(value)



## get the names of genes in RNAseq ----
#threshold <- quantile(resL_tb$AveExpr)[2]
rnaData <- resL_tb$gene

# check overlap
ggVennDiagram(list(proteome = proteomicsDataSYM, transcriptome = rnaData))

plProteins <- intersect(proteomicsDataSYM, rnaData)
# so 2028 transcripts are also in proteomics data, consistent with the protein number in 22869793.

## get platelet network ----
PlNetwork <- qc_interactions %>%
  filter(source_genesymbol %in% plProteins & 
           target_genesymbol %in% plProteins)

## check if all the network nodes are in either of the data

ggVennDiagram(list("proteome" = proteomicsDataSYM, "transcriptome" = rnaData, 
          "network" = unique(c(PlNetwork$source_genesymbol, 
                               PlNetwork$target_genesymbol)))) +
  scale_fill_gradient(low = "white", high = "white")+
  theme(legend.position = "none")

ggsave(filename = "../analysis/6_NetworkCons/PlateletNetwork/NodesVenn.tiff",
       device = "tiff",
       width = 8,
       height = 7)


## Filter interactions based on expression correlations ----
# Add source and target logFC columns using mutate and match
PlNetwork <- PlNetwork %>%
  mutate(
    source_logFC = resL_sig0$logFC[match(source_genesymbol, resL_sig0$gene)],
    target_logFC = resL_sig0$logFC[match(target_genesymbol, resL_sig0$gene)],
    multiply = source_logFC * target_logFC,
    predInt = case_when(
      multiply > 0 ~ "activation",
      multiply < 0 ~ "inhibition",
      is.na(multiply) ~ type
    )
  ) %>%
  # Filter the network based on matching type and predInt
  filter(type == predInt)

## write network to file ----
write.table(PlNetwork, "../analysis/6_NetworkCons/PlateletNetwork/PlNetwork.txt", 
            row.names = F, sep = "\t")

write.table(PlNetwork, "../analysis/6_NetworkCons/PlateletNetwork/PlPPINetwork.txt", 
            row.names = F, sep = "\t")


## prepare node annotations ----
PlNetworkNodeAnnots <- bind_rows(
  data.frame(uniprot = PlNetwork$source, symbol = PlNetwork$source_genesymbol),
  data.frame(uniprot = PlNetwork$target, symbol = PlNetwork$target_genesymbol)
) %>%
  unique()


write.table(PlNetworkNodeAnnots, "../analysis/6_NetworkCons/PlateletNetwork/PlNetworkUniprotSymbol.txt", 
            row.names = F, sep = "\t")

write.table(PlNetworkNodeAnnots, "../analysis/6_NetworkCons/PlateletNetwork/PlPPINetworkUniprotSymbol.txt", 
            row.names = F, sep = "\t")

# Make a DEG network ----
DEGnetwork <- PlNetwork %>%
  # Filter interactions where either source or target is in the significant gene list
  filter(source_genesymbol %in% resL_sig$gene | 
           target_genesymbol %in% resL_sig$gene) %>%
  # Map the logFC values and calculate the interaction type
  mutate(
    source_logFC = resL_sig0$logFC[match(source_genesymbol, resL_sig0$gene)],
    target_logFC = resL_sig0$logFC[match(target_genesymbol, resL_sig0$gene)],
    multiply = source_logFC * target_logFC,
    predInt = case_when(
      multiply > 0 ~ "activation",
      multiply < 0 ~ "inhibition",
      TRUE ~ type  # Default case when multiply is NA
    )
  ) %>%
  # Filter the network to keep only interactions where type matches predInt
  filter(type == predInt)

write.table(DEGnetwork, "../analysis/6_NetworkCons/DEGNetwork/DEGNetwork.txt", 
            row.names = F, sep = "\t")

## prepare node annotations ----
DEGnetworkNodeAnnots <- bind_rows(
  data.frame(uniprot = DEGnetwork$source, symbol = DEGnetwork$source_genesymbol),
  data.frame(uniprot = DEGnetwork$target, symbol = DEGnetwork$target_genesymbol)
) %>%
  unique()

write.table(DEGnetworkNodeAnnots, "../analysis/6_NetworkCons/DEGNetwork/DEGNetworkUniprotSymbol.txt", 
            row.names = F, sep = "\t")

unique(c(DEGnetwork$source, DEGnetwork$target))







# SNIPPETS ----
#GO enrichment
##obtained in cytoscape
firstneighbors <- c("HNRNPA0","TGFBR1","PPP3CA","PARN","ROCK2","PRKCD","AGO2","ALOX5","NSMCE3","NSMCE1","WASF1","CIB1","ITGA2B","ZFP36","SYK","TAX1BP1","AATF","DAPK3","MARK3","CDK2","LPCAT2","HSPB1","BAD","LYN","ACVRL1","BMPR1B","MYBL2","CD79A","CDKN1A","SKP2","RELA","ATM","SMAD2","FCGR2A","RPS6KA1","TP53","HSPA5","SMPD1","MAPK3","ERN1","ETV1","YWHAZ","MAPK1","ATF6","NFATC1","FBLIM1","UBE2J1","FYN","MYL9","SMAD1","MDM2","TNFAIP3","CIT","RCSD1","SMURF2","ARPC5","ROCK1","RNF168","LIMK1","ELAVL1","BAG2","CCR7","CDC42BPA","BLM","FOXO3","SRF","PIM1","TCF3","BAMBI","RPA1","POLA1","PRKCA","CREB1","PIAS1","SMAD4","MYC","SIL1","CDC25B","RNF111","ZFP36L1","PAK1","TSC2","SIN3A","MAP3K5","FLNA","DDX5","LSP1","SRC","HDAC1","SKIL","IFITM3","NFASC","ANK1","KLF10","MAP2K4","MYLK","TRAF6","MAPK14","MAPKAPK2")
PlClean <- read.csv("../analysis/networkCons/omnipathR_qc_Platelet_Clean_default_node.csv")

desc = "PlClean"
go_enrich_pl <- enrichGO(gene = PlClean$HGNC,
                         universe = resL_tb$gene,
                         OrgDb = org.Hs.eg.db, 
                         keyType = 'SYMBOL',
                         readable = T,
                         ont = "ALL",
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.2,
                         pAdjustMethod = "BH")


pl.tab = go_enrich_pl@result

write.table(pl.tab, file = paste0("../analysis/networkCons/", desc, "_GOBP_CUSTOM_all.txt"), sep = "\t", quote = F, 
            row.names = F, col.names = T)


tiff(filename = paste0("../analysis/networkCons/", desc, "_GOdotplot.tiff"),
     width = 10 * 300, 
     height = 7 * 300,
     res = 300,
     compression = "lzw")

dotplot(go_enrich_pl, title = paste(desc,"GO Overrepresentation", sep=" "), font.size=19)

dev.off()


tiff(filename = paste0("../analysis/networkCons/", desc, "_GObarplot.tiff"),
     width = 10 * 300, 
     height = 6 * 300,
     res = 300,
     compression = "lzw")

barplot(go_enrich_pl, 
        drop = TRUE, 
        x= "GeneRatio",
        showCategory = 10, 
        title = paste(desc, "GO Overrepresentation", sep=" "),
        font.size = 18,
)

dev.off()


# KEGG ---

ids_pl <-bitr(resL_tb$gene, fromType = "SYMBOL",
              toType = "ENTREZID", OrgDb=org.Hs.eg.db) 
ids_pl = ids_pl[!duplicated(ids_pl[c("SYMBOL")]),]

firstneighbors_entrez = bitr(firstneighbors, fromType = "SYMBOL",
                         toType = "ENTREZID", OrgDb=org.Hs.eg.db) 
firstneighbors_entrez = firstneighbors_entrez[!duplicated(firstneighbors_entrez[c("SYMBOL")]),]



kk_pl <- enrichKEGG(gene=firstneighbors_entrez$ENTREZID, 
                    universe=ids_pl$ENTREZID,
                    organism="hsa", 
                    pvalueCutoff = 0.05, 
                    keyType = "ncbi-geneid",
                    pAdjustMethod = "BH")



pl.tab_kk = kk_pl@result

write.table(pl.tab_kk, file = paste0("../analysis/networkCons/",desc, "_KEGG_CUSTOM_all.txt"), sep = "\t", quote = F, 
            row.names = F, col.names = T)


tiff(filename = paste0("../analysis/networkCons/", desc, "_KEGGdotplot.tiff"),
     width = 10 * 300, 
     height = 7 * 300,
     res = 300,
     compression = "lzw")

dotplot(kk_pl, title = paste(desc,"KEGG Enriched Pathways", sep=" "), font.size=19)

dev.off()


tiff(filename = paste0("../analysis/networkCons/", desc, "_KEGGbarplot.tiff"),
     width = 10 * 300, 
     height = 6 * 300,
     res = 300,
     compression = "lzw")

barplot(kk_pl, 
        drop = TRUE, 
        x= "GeneRatio",
        showCategory = 10, 
        title = paste(desc, "KEGG Enriched Pathways", sep=" "),
        font.size = 18,
)

dev.off()

pathway <- pl.tab_kk[pl.tab_kk$p.adjust < 0.05,]
i = length(pathway)

IDsOI <- pl.tab_kk[pl.tab_kk$Description == "Calcium signaling pathway", "geneID"]
protsOI <- unique(unlist(strsplit(IDsOI, split = "\\/")))
symbolsOI <- firstneighbors_entrez[firstneighbors_entrez$ENTREZID %in% protsOI,"SYMBOL"]


termOI <- AllGenes.AVG[rownames(AllGenes.AVG) %in% symbolsOI,] 

pheatmap(mat = termOI[,c(3,6)], 
         color = myheatcolors,
         scale="row",
         fontsize_col = 12,
         angle_col = 45,
         border_color = "black", 
         show_colnames = TRUE, 
         show_rownames = TRUE, 
         cellwidth = 25,
         cellheight = 12,
         drop_levels = TRUE, 
         fontsize = 12, 
         cluster_cols = F,
         cluster_rows = T,
         clustering_distance_cols=as.dist(1-cor(diffGenes.AVG, method="spearman")),
         cex=1,
         annotation_legend = F,
         annotation_names_row = F,
         main = "Calcium signaling pathway")

