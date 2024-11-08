# Description ----

# Load packages ----
suppressPackageStartupMessages({  
  library(tidyverse)
  library(limma)
  library(DT) 
  library(GSEABase) 
  library(Biobase) 
  library(GSVA) 
  library(gprofiler2) 
  library(clusterProfiler) 
  library(msigdbr) 
  library(enrichplot)
  })
 



# Carry out GO enrichment using gProfiler2 ----
# use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
myTopHitsUp <- resL_sig2[resL_sig2$logFC > 0,]
myTopHitsDown <- resL_sig2[resL_sig2$logFC < 0,]

##up and down
gost.resUp <- gost(myTopHitsUp$gene, 
                 organism = "hsapiens", 
                 correction_method = "fdr",
                 sources = c("GO", "KEGG", "REAC"),
                 custom_bg =rownames(resL_tb))

gostTableUp <- as.data.frame(t(unlist(gost.resUp[["result"]])))
write.table(gostTableUp, file = "../analysis/FunctionalEnrichment/ORA_UpGenes.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T) 

gost.resDown <- gost(myTopHitsDown$gene, 
                   organism = "hsapiens", 
                   correction_method = "fdr",
                   sources = c("GO", "KEGG", "REAC"),
                   custom_bg =rownames(resL_tb))

gostTableDown <- data.frame(apply(gost.resDown[["result"]],2,as.character))

write.table(gostTableDown, file = "../analysis/FunctionalEnrichment/ORA_DownGenes.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

gostTableDown <- gostTableDown[order(gostTableDown$p_value, decreasing = F),]

# produce an interactive manhattan plot of enriched GO terms
gostplot(gost.resUp, interactive = TRUE, capped = F)
#set interactive=FALSE to get plot for publications
mygostplotDown <- gostplot(gost.resDown, interactive = FALSE, capped = F)
# produce a publication quality static manhattan plot with specific GO terms highlighted

tiff("../analysis/FunctionalEnrichment/ORA_DownGenesPlot.tiff",
     width = 8 * 300,
     height = 8 * 300,
     res = 300,
     compression = "lzw")
publish_gostplot(
  mygostplotDown, #your static gostplot from above
  highlight_terms = gostTableDown[1:10,]$term_id,
  filename = NULL,
  width = NA,
  height = NA)
dev.off()


# Perform GSEA using clusterProfiler ----
# there are a few ways to get msigDB collections into R

#msigdbr_species()
hs_gsea <- msigdbr(species = "Homo sapiens") #gets all collections/signatures with human gene IDs
#take a look at the categories and subcategories of signatures available to you
print(hs_gsea %>%
  dplyr::distinct(gs_cat, gs_subcat) %>%
  dplyr::arrange(gs_cat, gs_subcat), n = 25)

# choose a specific msigdb collection/subcollection
# since msigdbr returns a tibble, we'll use dplyr to do a bit of wrangling
DBinput <- msigdbr(species = "Homo sapiens", 
                      category = "C2",
                      subcategory = "CP:WIKIPATHWAYS"
                   ) %>% # GO
  dplyr::select(gs_name, gene_symbol) 

# Now that you have your msigdb collections ready, prepare your data
# grab the dataframe you made in step3 script
# Pull out just the columns corresponding to gene symbols and LogFC for at least one pairwise comparison for the enrichment analysis
mydata.df.sub <- dplyr::select(resL_tb, gene, logFC)
# construct a named vector
mydata.gsea <- mydata.df.sub$logFC
names(mydata.gsea) <- as.character(mydata.df.sub$gene)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)

# run GSEA using the 'GSEA' function from clusterProfiler
set.seed(123) #set a random seed so that we can reproducible ordering for our GSEA results below
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=DBinput, verbose=FALSE) #could replace C2CP with hs_gsea_c2 object you retrieved from msigdb above
myGSEA.df <- as_tibble(myGSEA.res@result)

# view results as an interactive table
datatable(myGSEA.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE,
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:10), digits=2)
# create enrichment plots using the enrichplot package
gseaplot2(myGSEA.res,
          geneSetID = c(6,22,11,37), #can choose multiple signatures to overlay in this plot
          pvalue_table = FALSE, #can set this to FALSE for a cleaner plot
          #title = myGSEA.res$Description[14]
          ) #can also turn off this title

# add a variable to this result that matches enrichment direction with phenotype
myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "NSCLC",
    NES < 0 ~ "healthy"))


write.table(myGSEA.df, file = "../analysis/FunctionalEnrichment/GSEAWikiPathways.txt", sep = "\t", quote = F, 
            row.names = F, col.names = T)

myGSEA.nsclc <- myGSEA.df[myGSEA.df$phenotype == "NSCLC",]
myGSEA.healthy <- myGSEA.df[myGSEA.df$phenotype == "healthy",]

myGSEA.top <- rbind(myGSEA.nsclc[1:10,], myGSEA.healthy[1:10,])


# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(myGSEA.top, aes(x=phenotype, y=ID)) +
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  #scale_color_gradient(low="#87bcbd", high="#7c4b73") +
  scale_color_gradient(low="blue", high="red") +
  theme_bw()

ggsave("GSEAWikiPathways.tiff", 
       device = "tiff",
       path = "../analysis/FunctionalEnrichment/",
       width = 8, height = 6, bg = "white")


