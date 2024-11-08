# Description -----------
#this script creates heatmaps from differentially expressed genes or transcripts
#and selects modules of co-expressed genes based on pearson correlations

# Load packages -----
suppressPackageStartupMessages({  
  library(tidyverse)
  library(limma) #we only use limma in this script for the 'avearrays' function
  library(RColorBrewer) #need colors to make heatmaps
  library(gplots) #the heatmap2 function in this package is a primary tool for making heatmaps
  library(pheatmap)
  library(cowplot)
  library(plyr)
  library(gprofiler2)
  library(grid)
  })


# Choose Colors ----
myheatcolors <- colorRampPalette(colors=c("blue","white","red"))(100)

# Prepare Data----
vIn <- v.myDGEListIn.norm
results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=0.58)
disease <- myDGEListIn$samples$Patient_group
disease <- revalue(disease, c("Chronic_Pancreatitis" = "Chronic Pancreatitis",
                              "Epilepsy" = "Epilepsy",
                              "Healthy_Control" = "Healthy",
                              "Multiple_Sclerosis" = "MS",
                              "Non-significant_atherosclerosis" = "NS Atherosclerosis",
                              "NSCLC" = "NSCLC",
                              "Pulmonary_Hypertension" = "Pulmonary Hypertension",
                              "Stable_Angina_Pectoris" = "Stable AP",
                              "Unstable_Angina_Pectoris" =  "Unstable AP"))

colnames(vIn$E) <- disease
AllGenesHeatmap <- vIn$E
diffGenesHeatmap <- vIn$E[results[,1] !=0,]
dim(diffGenesHeatmap)

#now an old function from the limma package to average your replicates 
AllGenes.AVG <- avearrays(AllGenesHeatmap)
AllGenes.AVG <- AllGenes.AVG[,order(colnames(AllGenes.AVG))]

write.csv(AllGenes.AVG, file="../data/processed_data/LimmaCountsAveraged.csv")


diffGenes.AVG <- avearrays(diffGenesHeatmap)
diffGenes.AVG <- diffGenes.AVG[,order(colnames(diffGenes.AVG))]



tiff("../analysis/Clusters/LimmaHeatmapAllGenes.tiff",
     width = 8 *300, 
     height = 8 * 300,
     res = 300,
     compression = "lzw")

pheatmap(mat = AllGenes.AVG, 
         color = myheatcolors,,
         fontsize_col = 11,
         angle_col = 45,
         cellwidth = 30,
         cellheight = 0.02,
         scale="row",
         border_color = NA, 
         show_colnames = TRUE, 
         show_rownames = F, 
         #drop_levels = TRUE, 
         fontsize = 10, 
         cluster_cols = T,
         cluster_rows = T,
         clustering_distance_cols=as.dist(1-cor(AllGenes.AVG, method="spearman")),
         cex=1,
         main = "")

dev.off()


# Cluster the data using heatmap ----
tiff("../analysis/Clusters/LimmaClustersSum.tiff",
     width = 8 *300, 
     height = 6 * 300,
     res = 300,
     compression = "lzw")

set.seed(1508)
clustered_data <- pheatmap(diffGenes.AVG,
         myheatcolors, scale="row", na_col = "grey",
         border_color = "black", 
         fontsize_col = 11,
         cellwidth = 20,
         cellheight = 15,
         angle_col = 45,
         show_colnames = TRUE, 
         show_rownames = TRUE,
         drop_levels = TRUE, 
         fontsize = 10, 
         cluster_cols = TRUE,
         cluster_rows = F,
         cex=1,
         kmeans_k = 9, 
         clustering_distance_rows="correlation", ##note: we use Spearman, instead of Pearson, for clustering samples because it gives equal weight to highly vs lowly expressed transcripts or genes
         clustering_distance_cols=as.dist(1-cor(diffGenes.AVG, method="spearman")),
         clustering_method="complete",
         main = "")
dev.off()

# Add cluster information to your count matrix ----
data_subset_norm.clust <- cbind(diffGenes.AVG, 
                                cluster = clustered_data$kmeans$cluster)
data_subset_norm.clust <- as.data.frame(data_subset_norm.clust)

table(data_subset_norm.clust$cluster)

# Order count matrix by cluster number
data_subset_norm.clust.ordered <- data_subset_norm.clust[order(data_subset_norm.clust$cluster, decreasing = F),]
row_annots <- as.data.frame(paste0("mod ", data_subset_norm.clust.ordered$cluster))
rownames(row_annots) <- rownames(data_subset_norm.clust.ordered)
colnames(row_annots) = "module"

# Choose cluster colors ----
annot_colors <- list(module = brewer.pal(9, "PRGn"))
names(annot_colors$module) <- unique(row_annots$module)
annot_colors$module[5] <- "burlywood"

# Make heatmap showing all clusters ----

tiff("../analysis/Clusters/LimmaClusters.tiff",
     width = 16 *300, 
     height = 12 * 300,
     res = 300,
     compression = "lzw")

pheatmap(mat = data_subset_norm.clust.ordered[1:9], 
         color = myheatcolors,,
         fontsize_col = 11,
         angle_col = 45,
         cellwidth = 40,
         cellheight = 2,
         gaps_row=cumsum(as.numeric(table(row_annots$module))),
         scale="row",
         border_color = NA, 
         show_colnames = TRUE, 
         show_rownames = F, 
         annotation_colors = annot_colors,
         annotation_row = row_annots,
         annotation_legend = T,
         annotation_names_row = F,
         #drop_levels = TRUE, 
         fontsize = 10, 
         cluster_cols = T,
         cluster_rows = F,
         clustering_distance_cols=as.dist(1-cor(diffGenes.AVG, method="spearman")),
         cex=1,
         main = "")

dev.off()




# View modules of co-regulated genes ----
# view your color assignments for the different clusters

# look at these clusters in more detail later
module.assign <- data_subset_norm.clust.ordered$cluster
names(module.assign) <- rownames(data_subset_norm.clust.ordered)

#now assign a color to each module (makes it easy to identify and manipulate)
module.color <- brewer.pal(9, "PRGn")
module.color <- module.color[as.vector(module.assign)] 

names(module.color) <- names(module.assign) 

module.assign.df <- as_tibble(as.list(module.assign))

module.assign.pivot <- pivot_longer(module.assign.df, # dataframe to be pivoted
                          cols = everything(), # column names to be stored as a SINGLE variable
                          names_to = "geneID", # name of that new variable (column)
                          values_to = "module") # name of new variable (column) storing all the values (data)

module.assign.pivot <- module.assign.pivot %>%
  mutate(moduleColor = case_when(
    module == 1 ~ "#762A83",
    module == 2 ~ "#9970AB",
    module == 3 ~ "#C2A5CF",
    module == 4 ~ "#E7D4E8",
    module == 5 ~ "burlywood",
    module == 6 ~ "#D9F0D3",
    module == 7 ~ "#A6DBA0",
    module == 8 ~ "#5AAE61",
    module == 9 ~ "#1B7837"))

ggplot(module.assign.pivot) +
  aes(as.character(module)) +
  geom_bar(aes(fill=moduleColor)) +
  labs(title = "Number of genes in clusters",
       x ="cluster")+
  scale_fill_identity()+
  theme_cowplot()

ggsave("LimmaModuleSizes.tiff", 
       device = "tiff",
       path = "../analysis/Clusters/",
       width = 6, height = 4, bg = "white")

# Heatmaps of modules ----

i = 1
for (i in unique(module.assign)) {
  myModule <- diffGenes.AVG[names(module.assign[module.assign %in% i]),] 
  assign(paste0("myModule", i), myModule)
  tiff(paste0("../analysis/Clusters/Module", i, "Heatmap2.tiff"),
     width = 10 *300, 
     height = 10 * 300,
     res = 300,
     compression = "lzw")
  pheatmap(mat = myModule, 
         color = myheatcolors,
         scale="row",
         fontsize_col = 12,
         angle_col = 45,
         border_color = "black", 
         show_colnames = TRUE, 
         show_rownames = TRUE, 
         cellwidth = 25,
         cellheight = 12,
         annotation_colors = annot_colors,
         annotation_row = row_annots,
         drop_levels = TRUE, 
         fontsize = 12, 
         cluster_cols = T,
         cluster_rows = T,
         clustering_distance_cols=as.dist(1-cor(diffGenes.AVG, method="spearman")),
         cex=1,
         annotation_legend = F,
         annotation_names_row = F,
         main = paste0("Expression of module ", i, " genes"))

dev.off()

}


# GO ORA of modules ----


logScale <- function(input, input_start = 1, input_end = 50000, 
                     output_start = 2, output_end = 10) {
  m = (output_end - output_start)/(log(input_end) - log(input_start))
  b = -m * log(input_start) + output_start
  output = m * log(input) + b
  return(output)
}

# Carry out GO enrichment using gProfiler2 ----

datalist = list()

for (i in 1:9) {
  df <- get(paste0("myModule", i))
  gost.res <- gost(rownames(df), 
                   organism = "hsapiens", 
                   correction_method = "fdr",
                   sources = c("GO", "KEGG", "REAC"),
                   custom_bg =rownames(resL_tb))
  if (is.null(gost.res)) {
    next
  } 
  
  gostTable <- data.frame(apply(gost.res[["result"]],2,as.character))
  
  essential_names <- c("source_order", "term_size", "term_name", 
                       "term_id", "source", "significant", "p_value")
  gostTable <- gostTable %>% dplyr::select(essential_names)
  gostTable <- gostTable %>% mutate_at(c('source_order', 'term_size', 'p_value'), as.numeric)
  gostTable$term_size_scaled <- logScale(gostTable$term_size)
  gostTable$module <- as.character(i)
  gostTable <- gostTable[order(gostTable$p_value),]
  gostTable <- gostTable[1:5,]
  datalist[[i]] <- gostTable
}

myGO.df = do.call(rbind, datalist)

myGO.df <- myGO.df %>%
  mutate(moduleColor = case_when(
    module == 1 ~ "#762A83",
    module == 2 ~ "#9970AB",
    module == 3 ~ "#C2A5CF",
    module == 4 ~ "#E7D4E8",
    module == 5 ~ "burlywood",
    module == 6 ~ "#D9F0D3",
    module == 7 ~ "#A6DBA0",
    module == 8 ~ "#5AAE61",
    module == 9 ~ "#1B7837"))

myGO.df[12,]$term_name <- "Regulation of IGF transport and uptake by IGFBPs" 

ggplot(myGO.df, aes(x = module, y = term_name)) +
  geom_point(aes(size = term_size_scaled,
                 fill = moduleColor, alpha = -log10(p_value)),
             shape = 21) +
  scale_fill_identity() +  
  scale_alpha(range = c(0.4, 1)) +
  theme_bw()

ggsave("GOBPORA_Clusters.tiff", 
       device = "tiff",
       path = "../analysis/Clusters/",
       width = 7, height = 7, bg = "white")

gostplot(gost.res, interactive = TRUE, capped = F)

# mark the genes in GOtermsOI in the module heatmaps
library(grid)
interesting_GOs <- list("cell surface" = "GO:0009986",
                        "external side of plasma membrane" = "GO:0009897")
i = 1
description <- names(interesting_GOs[i])
results <- AnnotationDbi::select(org.Hs.eg.db, keys=c(interesting_GOs[[i]]), columns = c('SYMBOL'), keytype = "GOALL")
symbols1 <- unique(results$SYMBOL)

# Extract and rearrange 
inGO <- intersect(rownames(myModule2), symbols)
myModule <- as.data.frame(myModule2)
myModule$color <- ifelse(rownames(myModule) %in% inGO, "#D73027", "black")
mat_text_colors <- data.frame(myModule$color)
rownames(mat_text_colors) <- rownames(myModule)

e <- pheatmap(mat = myModule[1:9], 
         color = myheatcolors,
         scale="row",
         fontsize_col = 12,
         angle_col = 45,
         border_color = "black", 
         show_colnames = TRUE, 
         show_rownames = TRUE, 
         cellwidth = 25,
         cellheight = 12,
         annotation_colors = annot_colors,
         annotation_row = row_annots,
         drop_levels = TRUE, 
         fontsize = 12, 
         cluster_cols = T,
         cluster_rows = T,
         clustering_distance_cols=as.dist(1-cor(diffGenes.AVG, method="spearman")),
         cex=1,
         annotation_legend = F,
         annotation_names_row = F,
         main = paste0("Expression of module 2 "," genes"))

cols=mat_text_colors[order(match(rownames(mat_text_colors), e$gtable$grobs[[6]]$label)), ]

e$gtable$grobs[[6]]$gp=gpar(col=cols)

tiff(paste0("../analysis/Clusters/Module", gsub(" ", "",description), "Heatmap.tiff"),
     width = 10 *300, 
     height = 10 * 300,
     res = 300,
     compression = "lzw")
print(e)
dev.off()

# Tables of modules
library(ReactomeContentService4R)

PlateletDeg <- event2Ids(event.id = "R-HSA-114608")$geneSymbol
PlateletResptoCa2 <- event2Ids(event.id = "R-HSA-76005")$geneSymbol
PlateletAct <- event2Ids(event.id = "R-HSA-76002")$geneSymbol
venn(list(PD=PlateletDeg, PA=PlateletAct, PR=PlateletResptoCa2)) #PA covers all.

Insulin <- event2Ids(event.id = "R-HSA-381426")$geneSymbol


module2_df <- resL_sig2[resL_sig2$gene %in% rownames(myModule3),] %>% 
  dplyr::select(gene, description,AveExpr, logFC, P.Value, adj.P.Val) %>%
  arrange(desc(abs(logFC)))

module2_df$P.Value <- sprintf('%.2e', module2_df$P.Value)
module2_df$adj.P.Val<- sprintf('%.2e', module2_df$adj.P.Val)


module2_df %>%
gt()  %>%
  fmt_number(columns=c(3:6), decimals = 2) %>%
  cols_width(
    gene ~ px(300),
    P.Value ~ px(200),
    adj.P.Val ~ px(200),
    description ~ px(1000)) %>%
  cols_label(gene = md("**gene**"),
             AveExpr = md("**Avg. exp.**"),
             logFC = md("**log2FC**"),
             P.Value = md("**p-value**"),
             adj.P.Val = md("**adj. p-value**"),
             description = md("**description**")
  ) %>%
  cols_align("center",
             columns = 3:6) %>%
  tab_header(title = md("**Genes in Module 3**"),
             #subtitle = md("**to choose fit type***")
  )  %>%
  tab_source_note(source_note = "red: Platelet pathways, blue: Insulin Growth Factor (IGF) transport") %>%
  tab_style(
    style = cell_text(color = "red"),
    locations = cells_body(
      columns = gene,
      rows = gene %in% PlateletAct
    )
  ) %>% 
  tab_style(
    style = cell_text(color = "blue"),
    locations = cells_body(
      columns = gene,
      rows = gene %in% Insulin
    )
  ) %>%
  tab_options(table.font.size = 10,
              data_row.padding = px(1)) %>%
  gtsave("../analysis/Clusters/module2Table.pdf")




# HEATMAPS OF GSEA RESULTS ----

#get pathway of interest
interestingpathwaysUP <- c("hsa04512", "hsa04974","hsa05033", "hsa04820", "hsa04510", "hsa04814",
                         "hsa00430", "hsa00512", "hsa03320", "hsa04020","hsa04911",
                         "hsa04024")

interestingpathwaysDOWN <- c("hsa03010", "hsa05171", "hsa03040", "hsa05140", "hsa00190",
                             "hsa04672", "hsa04612", "hsa03050", "hsa04662", "hsa05320", "hsa00020",
                             "hsa04210", "hsa04064", "hsa04660","hsa04620", "hsa04061")

category <- "DOWN"
listTerms <- get(paste0("interestingpathways", category))


for (term in listTerms){
  description <- pl.tab_kk[term,]$Description
  genes <- unlist(strsplit(pl.tab_kk[term,]$core_enrichment, split = "/"))
  genesName <- df[df$entrezid %in% genes,"gene"]
  myPathway <- AllGenes.AVG[genesName,]
  
  # Extract and rearrange 
  DEGsRelaxed <- resL_sig0[abs(resL_sig0$logFC) > 0.2,]$gene
  DEGs <- resL_sig$gene
  

  
  myPathway <- as.data.frame(myPathway)
  myPathway$color <- ifelse(rownames(myPathway) %in% setdiff(DEGsRelaxed,DEGs), "black",
                            ifelse(rownames(myPathway) %in% DEGs, "#D73027","gray"))
  mat_text_colors <- data.frame(myPathway$color)
  rownames(mat_text_colors) <- rownames(myPathway)
  
  e <- pheatmap(mat = myPathway[1:9], 
                color = myheatcolors,
                scale="row",
                fontsize_col = 12,
                angle_col = 45,
                border_color = "black", 
                show_colnames = TRUE, 
                show_rownames = TRUE, 
                cellwidth = 25,
                cellheight = 10,
                #annotation_colors = annot_colors,
                #annotation_row = row_annots,
                drop_levels = TRUE, 
                fontsize = 12, 
                cluster_cols = T,
                cluster_rows = T,
                clustering_distance_cols=as.dist(1-cor(diffGenes.AVG, method="spearman")),
                cex=1,
                annotation_legend = F,
                annotation_names_row = F,
                main = paste0("Expression of ",description ," genes"))
  
  cols=mat_text_colors[order(match(rownames(mat_text_colors), e$gtable$grobs[[6]]$label)), ]
  
  e$gtable$grobs[[6]]$gp=gpar(col=cols)
  
  
  tiff(paste0("../analysis/4_GSEA/KEGG/Heatmaps/", category, "/", gsub(" ", "_", description), ".tiff"),
       width = 10 *300, 
       height = 18 * 300,
       res = 300,
       compression = "lzw")
  print(e)
  dev.off()
}

# get manual proteins heatmap



integrins <- annotations[grep("integrin", annotations$description),"gene"]
glycoproteins <- annotations[grep("glycoprotein", annotations$description),"gene"]

description <- "Glycoproteins and related"
genesName <- integrins
myPathway <- AllGenes.AVG[genesName,]

DEGsRelaxed <- resL_sig0[abs(resL_sig0$logFC) > 0.2,]$gene
DEGs <- resL_sig$gene

myPathway <- as.data.frame(myPathway)
myPathway$color <- ifelse(rownames(myPathway) %in% setdiff(DEGsRelaxed,DEGs), "black",
                          ifelse(rownames(myPathway) %in% DEGs, "#D73027","gray"))
mat_text_colors <- data.frame(myPathway$color)
rownames(mat_text_colors) <- rownames(myPathway)

e <- pheatmap(mat = myPathway[1:9], 
              color = myheatcolors,
              scale="row",
              fontsize_col = 12,
              angle_col = 45,
              border_color = "black", 
              show_colnames = TRUE, 
              show_rownames = TRUE, 
              cellwidth = 25,
              cellheight = 10,
              #annotation_colors = annot_colors,
              #annotation_row = row_annots,
              drop_levels = TRUE, 
              fontsize = 12, 
              cluster_cols = T,
              cluster_rows = T,
              clustering_distance_cols=as.dist(1-cor(diffGenes.AVG, method="spearman")),
              cex=1,
              annotation_legend = F,
              annotation_names_row = F,
              main = paste0("Expression of ",description ," genes"))

cols=mat_text_colors[order(match(rownames(mat_text_colors), e$gtable$grobs[[6]]$label)), ]

e$gtable$grobs[[6]]$gp=gpar(col=cols)


tiff(paste0("../analysis/4_GSEA/KEGG/Heatmaps/", gsub(" ", "_", description), ".tiff"),
     width = 10 *300, 
     height = 18 * 300,
     res = 300,
     compression = "lzw")
print(e)
dev.off()
