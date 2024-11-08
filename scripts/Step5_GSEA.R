
#load packages ----
suppressPackageStartupMessages({  
  library(clusterProfiler)
  library(enrichplot)
  library(ggplot2)
  library(pathview)
  library(ggridges)
  library(dplyr)
  })

# SET THE DESIRED ORGANISM HERE
organism = "org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)


df <- resL_tb
df <- merge(df, TxOG, by='gene', all.x = T) %>%
  distinct(., gene, .keep_all = T)
df$entrezid <- unlist(sapply(df$entrezid,"[[",1))

# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector - 
#according to the column in your dataset containing the gene symbol or id
names(original_gene_list) <- df$entrezid
#names(original_gene_list) <- df$Symbol

# omit any NA values 
gene_list<-na.omit(original_gene_list)
removeNAnames <- names(gene_list) == "NA"
gene_list <- gene_list[!removeNAnames]

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
kegg_organism = "hsa"

#After here, you can choose and run the analysis & visualizations according to your needs
#1 KEGG GSEA####

kk2 <- gseKEGG(geneList     = gene_list,
               organism     = kegg_organism,
               minGSSize    = 10,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid"
)

kk3 <- gseMKEGG(geneList     = gene_list,
               organism     = kegg_organism,
               minGSSize    = 5,
               maxGSSize    = 500,
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               keyType       = "ncbi-geneid"
)


pl.tab_kk = kk2@result

datatable(pl.tab_kk,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 2: GSEA Resultss',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(4:8), digits=2)


write.table(pl.tab_kk, file = "../analysis/GSEA/KEGG/KEGG_CUSTOM_allLimma.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

#1a. Visualizations : dotplot & ridgeplot
#choosing manually

UPs <- c("hsa04512", "hsa04510","hsa02010","hsa00430", "hsa00512",
         "hsa04020", "hsa03320", "hsa04911","hsa04024", "hsa04010")
DOWNs <- c("hsa03010", "hsa03040", "hsa00190","hsa04612", "hsa04662",
           "hsa03430", "hsa00020","hsa04064","hsa04660","hsa04210")

kk <- kk2
kk@result <- kk@result[kk@result$ID %in% UPs | kk@result$ID %in% DOWNs, ]

tiff(filename = "../analysis/GSEA/KEGG/dotplotLimma2.tiff",
     width = 10 * 300, 
     height = 6 * 300,
     res = 300,
     compression = "lzw")

options(enrichplot.colours = c("lavenderblush1","blue"))
dotplot(kk, title = "KEGG GSEA NSCLC vs. Noncancer", font.size=12,
        showCategory = 10, split=".sign",) + facet_grid(.~.sign)

dev.off()


tiff(filename = "../analysis/GSEA/KEGG/ridgeplotLimma2.tiff",
     width = 10 * 300, 
     height = 11 * 300,
     res = 300,
     compression = "lzw")
options(enrichplot.colours = c("lavenderblush1","blue"))
ridgeplot(kk) + labs(x = "enrichment distribution")

dev.off()

# Pathway specific heatmaps ----

# Define the term and extract relevant gene information
term <- c("Focal adhesion", "ECM-receptor interaction")
select <- unlist(strsplit(kk@result[kk@result$Description == term[2], "core_enrichment"], "/"))
select_genes <- df[df$entrezid %in% select, "gene"]

# Extract and rearrange module data
myModule <- as.data.frame(AllGenes.AVG[select_genes, ])
myModule <- myModule[, c(1, 2, 4, 5, 7, 8, 9, 3, 6)]
myModule <- myModule[order(row.names(myModule)), ]

venn(list(rownames(myModule2), rownames(myModule)))

# Prepare results data frame
resL_df <- as.data.frame(resL_tb)
rownames(resL_df) <- resL_df$gene
resL_select <- resL_df[select_genes, ]

# Generate p-values and logFC matrices
targets_pvals <- as.data.frame(resL_select$adj.P.Val)
rownames(targets_pvals) <- rownames(resL_select)

targets_pvals[targets_pvals <= 0.05] <- "*"
targets_pvals[targets_pvals> 0.05] <- ""
other_pvals <- matrix(nrow = 58, ncol = 8)
targets_pvals <- cbind(as.data.frame(other_pvals), targets_pvals)
targets_pvals[is.na(targets_pvals)] <- ""


targets_logfcs <-as.data.frame(resL_select$logFC)
targets_logfcs<-round(targets_logfcs, digits = 2)
other_logFCs <- matrix(nrow = 58, ncol = 7)
targets_logfcs <- cbind(as.data.frame(other_logFCs), targets_logfcs)
targets_logfcs[is.na(targets_logfcs)] <- ""

for (j in 1:(length(colnames(targets_pvals))-1)) {
  print(j)
  for (i in 1:length(rownames(targets_logfcs))) {
    if (targets_pvals[i, j+1] == "") {
      targets_logfcs[i, j] <- ""}}}

targets_logfcs$V8<- ""
targets_logfcs <- targets_logfcs[c(1:7,9,8)]
test_labels <- as.matrix(targets_logfcs) 


tiff(paste0("../analysis/GSEA/KEGG/", gsub(" ", "",term), "Heatmap.tiff"),
     width = 10 *300, 
     height = 10 * 300,
     res = 300,
     compression = "lzw")
pheatmap(mat = myModule, color = myheatcolors,
         scale="row", fontsize_col = 12,
         display_numbers = test_labels, number_color = "darkgreen", 
         angle_col = 45, border_color = "black", 
         show_colnames = TRUE, show_rownames = TRUE, 
         cellwidth = 25, cellheight = 7,
         annotation_colors = annot_colors,
         drop_levels = TRUE, fontsize = 8, 
         cluster_cols = F, cluster_rows = T,
         cex=1, annotation_legend = F, annotation_names_row = F,
         main = paste0("Expression of ", term, " genes"))
dev.off()


# 2. Reactome GSEA ----
library(ReactomePA)

gse.r <- gsePathway(
  geneList=gene_list,
  #nPerm=100000,
  organism = "human",
  minGSSize=10,
  maxGSSize =100,
  pvalueCutoff=1, #to get the table, filtering can be done after
  pAdjustMethod="BH",
  verbose=TRUE)


gse.r_genename <- setReadable(gse.r, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

write.table(gse.r_genename, file = "../analysis/GSEA/REACTOME/Reactome_all.txt", sep = "\t", quote = F, 
              row.names = F, col.names = T)

dot.r<-dotplot(gse.r_genename, showCategory=10, split=".sign") + facet_grid(.~.sign)

tiff(filename = "../analysis/GSEA/REACTOME/Reactome_dotplot.tiff",
     width = 10 * 300, 
     height = 12 * 300,
     res = 300,
     compression = "lzw")
print(dot.r)
dev.off()


 #3. Gene Ontology GSEA####

# we want the log2 fold change 
original_gene_list <- df$logFC

# name the vector - change the df$_ 
#according to the column in your dataset containing the gene symbol or id
#names(original_gene_list) <- df$GeneID
names(original_gene_list) <- df$gene


# omit any NA values 
gene_list<-na.omit(original_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

head(gene_list)

#analysis - change keytype accordingly, minGSsize is minimum number of genes forming a geneset
gse <- gseGO(geneList=gene_list, 
             ont ="CC", 
             keyType = "SYMBOL", 
             minGSSize = 10, 
             maxGSSize = 500, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "BH")

pl.tab = gse@result


# 
# 
# ##check terms of interest
# pl.tab_all = gse_all@result
# pl.tab_all[grep(terms_of_interest[1], pl.tab_all$ID),]



write.table(pl.tab, file = "../analysis/GSEA/GO/GOBP_CUSTOM_MF.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

#1a. Visualizations : dotplot & ridgeplot



tiff(filename = paste0("../analysis/GSEA/GO/dotplotCC.tiff"),
     width = 10 * 300, 
     height = 8 * 300,
     res = 300,
     compression = "lzw")

dotplot(gse, showCategory = 10, title = "GO CC Overrepresentation", font.size=11,
        split=".sign") +
  facet_grid(.~.sign)

dev.off()


tiff(filename = "../analysis/GSEA/GO/ridgeplotMF.tiff",
     width = 10 * 300, 
     height = 15 * 300,
     res = 300,
     compression = "lzw")

ridgeplot(gse) + labs(x = "enrichment distribution")

dev.off()


