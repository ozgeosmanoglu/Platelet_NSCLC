
# Load packages -----
suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(matrixStats)
  library(cowplot)
  library(RUVSeq)
  library(plotly)
  library(DESeq2)
  library(DEFormats)
  library(glmGamPoi)
  library(swamp)
  library(ggpubr)
  library(MetBrewer)
  library(gt)
  library(Glimma)
  library(DT)
  library(dplyr)
  library(tibble)
  library(org.Hs.eg.db)
})


#keep only cancer samples 
ddsCancer <- ddsColl[,ddsColl$class == "NSCLC"]
ddsCancer <- ddsCancer[, !is.na(ddsCancer$Metastasis)]
design(ddsCancer) <-  ~Metastasis

#Filtering low counts! ----
y <- as.DGEList(ddsCancer)
keep <- filterByExpr(y) #design ~class 
table(keep)

y <- y[keep,] #14258 21039
ddsFilt <- ddsCancer[keep,]

# Make a DGElist from your counts ----

myDGEListFilt <- as.DGEList(ddsFilt)
myDGEListFilt$samples <- droplevels(myDGEListFilt$samples)

#Let's explore the data! ----

##PCA plots ----

###variance-stabilizing transformation----
sub = "Filt"
desc = "filtered data (variance stabilizing transformation)"

vsd <- vst(ddsFilt, blind=F)
pc <- prcomp(t(assay(vsd)))

#calculate percentages
pc.var<-pc$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
head(pc.per)

pca.res.df <- as_tibble(pc$x)
class = vsd$Metastasis

myplot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color = class, text = vsd$Run)+
  geom_point()+
  theme_bw() +
  #geom_label(nudge_y = 0, nudge_x = 0) +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")"))+
  ylab(paste0("PC2 (",pc.per[2],"%",")"))+ 
  labs(title=paste0("PCA plot of ", desc),
       caption=paste0("produced on ", Sys.time()))+
  coord_fixed()


ggsave(paste0("PCA", sub,"Data.tiff"), 
       myplot, 
       device = "tiff",
       path = "../analysis/Metastasis /",
       width = 8, height = 6, bg = "white")


ggplotly(myplot) #interactive plot

###estimate RUV factors + technical variation----
#https://github.com/mikelove/preNivolumabOnNivolumab/blob/main/preNivolumabOnNivolumab.knit.md

#First we see what happens when we run simple DESeq2 analysis comparing the two groups, 
#without attempting to control for the technical variation. Here we use glmGamPoi which 
#is an efficient method for estimating dispersion when we have many samples. This can 
#speed up the analysis by an order of magnitude. For details, see: https://doi.org/10.1093/bioinformatics/btaa1009

ddsIn <- ddsFilt
ddsIn <- DESeq(ddsIn, test="LRT", reduced=~1, fitType="glmGamPoi")

res <- results(ddsIn)
table(res$padj < .05) 
#FALSE  TRUE 
#7381 13658 
#There are thousands of DE genes. As we will see, these are the result of uncontrolled confounding.
#Ignoring the technical variation is not appropriate, if there is correlation 
#between the technical variation and the condition. 

table(abs(res$log2FoldChange) < .05 & res$padj > .05) #empiricals

set <- newSeqExpressionSet(counts(ddsIn))
set <- betweenLaneNormalization(set, which="upper")
not_sig <- rownames(res)[which(res$padj > .05 & abs(res$log2FoldChange) < .05)]
empirical <- rownames(set)[ rownames(set) %in% not_sig ]
set <- RUVg(set, empirical, k=5)

pdat <- pData(set)
pdat$class <- class

vsd$W1 <- pdat$W_1
vsd$W2 <- pdat$W_2
vsd$W3 <- pdat$W_3
vsd$W4 <- pdat$W_4
vsd$W5 <- pdat$W_5


#We can visualize how the factors of unwanted variation describe the samples 
#in the PC1 and PC2 space:

sub = "vstF"
desc = "Unwanted Variation in Filtered Data (vst, outlier removed)"

pc <- prcomp(t(assay(vsd)))

#calculate percentages
pc.var<-pc$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
head(pc.per)
pca.res.df <- as_tibble(pc$x)

##visualize the unwanted variation

hospitalplot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color = vsd$Hospital, text = vsd$Run)+
  geom_point()+
  theme_bw() +
  #geom_label(nudge_y = 0, nudge_x = 0) +
  #stat_ellipse() +
  guides(color=guide_legend(title="hospital"))+
  xlab(paste0("PC1 (",pc.per[1],"%",")"))+
  ylab(paste0("PC2 (",pc.per[2],"%",")"))+ 
  coord_fixed()+
  scale_color_brewer(palette = "Set1")  

p1 <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color = vsd$W1, text = vsd$Run)+
  geom_point()+
  theme_bw() +
  #geom_label(nudge_y = 0, nudge_x = 0) +
  #stat_ellipse() +
  guides(color=guide_colourbar(title="UV1"))+
  xlab(paste0("PC1 (",pc.per[1],"%",")"))+
  ylab(paste0("PC2 (",pc.per[2],"%",")"))+ 
  coord_fixed()+
  scale_color_met_c("VanGogh1")

p2 <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color = vsd$W2, text = vsd$Run)+
  geom_point()+
  theme_bw() +
  #geom_label(nudge_y = 0, nudge_x = 0) +
  #stat_ellipse() +
  guides(color=guide_colourbar(title="UV2"))+
  xlab(paste0("PC1 (",pc.per[1],"%",")"))+
  ylab(paste0("PC2 (",pc.per[2],"%",")"))+ 
  coord_fixed()+
  scale_color_met_c("VanGogh1")

p3 <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color = vsd$W3, text = vsd$Run)+
  geom_point()+
  theme_bw() +
  guides(color=guide_colourbar(title="UV3"))+
  #geom_label(nudge_y = 0, nudge_x = 0) +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")"))+
  ylab(paste0("PC2 (",pc.per[2],"%",")"))+ 
  coord_fixed()+
  scale_color_met_c("VanGogh1")

p4 <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color = vsd$W4, text = vsd$Run)+
  geom_point()+
  theme_bw() +
  guides(color=guide_colorbar(title="UV4"))+
  #geom_label(nudge_y = 0, nudge_x = 0) +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")"))+
  ylab(paste0("PC2 (",pc.per[2],"%",")"))+ 
  coord_fixed()+
  scale_color_met_c("VanGogh1")


p5 <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color = vsd$W5, text = vsd$Run)+
  geom_point()+
  theme_bw() +
  guides(color=guide_colorbar(title="UV5"))+
  #geom_label(nudge_y = 0, nudge_x = 0) +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")"))+
  ylab(paste0("PC2 (",pc.per[2],"%",")"))+ 
  coord_fixed()+
  scale_color_met_c("VanGogh1")



figure <- ggarrange(p1,p2,p3,p4,p5,hospitalplot,
                    #labels = c("UV1", "UV2", "UV3", "UV4", "UV5"),
                    #font.label = list(size = 10, face = "bold"),
                    ncol = 2, nrow = 3, vjust = 1.5, hjust = -1)



figureAnnot <- annotate_figure(figure,
                               top = text_grob("Unwanted Variation in Filtered Data", face = "bold", size = 12),
                               bottom = text_grob("(vst, outlier removed)",
                                                  hjust = 1, x = 1, size = 10))

ggsave(paste0("PCA", sub,"DataUVs.tiff"), 
       figureAnnot, 
       device = "tiff",
       path = "../analysis/Metastasis/",
       width = 10, height = 7, bg = "white")


##Prince plots----
sub = "vst" #Raw or Coll or Filt
mydds <- ddsFilt #ddsRaw or ddsColl
myswamp <- mydds@assays@data@listData[["counts"]]  

set.seed(19)
o1<-data.frame(
  Hospital = mydds$Hospital,
  StorageTime = mydds$Storage_time, 
  W1 = vsd$W1,
  W2 = vsd$W2,
  W3 = vsd$W3,
  W4 = vsd$W4,
  W5 = vsd$W5,
  Gender = mydds$Gender,
  Age = mydds$Age,
  Smoking = mydds$Smoking,
  Treatment = mydds$Treatment,
  Metastasis = mydds$Metastasis,
  SubmissionDate = mydds$submission_date,
  DataProvider = mydds$provider,
  Random = rnorm(ncol(myswamp)), row.names=colnames(myswamp))


# PCA analysis
res1<-prince(myswamp,o1,top=10,permute=F)
#beep(sound=2)

prince.plot(prince=res1, key = T, note = T, breaks = 16, colsep = 1:9,
            col =  rev(colorRampPalette(c("white", "lightgray", 
                                          "gold", "orange", "red"))(15)))


tiff(paste0("../analysis/Metastasis/PrincePlot", sub,"Data_UVs.tiff"),
     width = 10 * 300,
     height = 6 * 300,
     res = 300,
     compression = "lzw")
prince.plot(prince=res1, key = T, note = T, breaks = 16, colsep = 1:9,
            col =  rev(colorRampPalette(c("white", "lightgray", "gold", "orange", "red"))(15)),
            #               col =  heat.colors(15)
)
dev.off()

# Differential Expression ----
colData(ddsIn) <- cbind(colData(ddsIn), pdat[,1:5])

myDGEListIn <- as.DGEList(ddsIn)
myDGEListIn.norm <- calcNormFactors(myDGEListIn, method = "TMM")
myDGEListIn.norm$samples <- droplevels(myDGEListIn.norm$samples)

group <- myDGEListIn.norm$samples$group
W_1 <- myDGEListIn.norm$samples$W_1
W_2 <- myDGEListIn.norm$samples$W_2
W_3 <- myDGEListIn.norm$samples$W_3
W_4 <- myDGEListIn.norm$samples$W_4
W_5 <- myDGEListIn.norm$samples$W_5


design <- model.matrix(~0 + group + W_1 + W_2 + W_3 + W_4 + W_5)
colnames(design)[1:2] <- c("primary", "metastasis")


v.myDGEListIn.norm <- voom(myDGEListIn.norm, design, plot = TRUE)

#fit linear model
fit <- lmFit(v.myDGEListIn.norm, design)
head(fit$coefficients)

# Contrast matrix ----
contrast.matrix <- makeContrasts(metastasis - primary, levels = design) #limma

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=0.58)

plotMD(ebFit, status=results)

top.table <- topTable(ebFit, adjust ="BH", sort.by = "P", n = Inf)

res_tb <- top.table %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()


diffGenes <- v.myDGEListIn.norm$E[results[,1] !=0,]
#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

#get sig dfs----
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

res_sig <- res_tb %>%
  dplyr::filter(`adj.P.Val` < padj.cutoff & abs(logFC) > lfc.cutoff)

res_sig0 <- res_tb %>%
  dplyr::filter(`adj.P.Val` < padj.cutoff)

##make the top table----

#add descriptions
resL_sig2 <- merge(resL_sig, annotations, by='gene', all.x = T) %>% 
  distinct(., gene, .keep_all = T)

resL_sig2$description <- unlist(lapply(strsplit(as.character(resL_sig2$description), "\\["), '[[', 1))

# save 
# Arrange and merge in a single step
orderedRes <- res_tb %>%
  arrange(adj.P.Val) %>%
  left_join(annotations, by = 'gene') %>%
  distinct(gene, .keep_all = TRUE) %>%
  left_join(TxOG[6:8], by = 'gene') %>%
  distinct(gene, .keep_all = TRUE)

# Add UniProt IDs
orderedRes$uniprotid <- mapIds(org.Hs.eg.db, keys = orderedRes$gene, keytype = "SYMBOL", column = "UNIPROT", multiVals = "first")

# Clean description column and select relevant columns
orderedRes <- orderedRes %>%
  mutate(description = sub("\\[.*", "", description)) %>%
  dplyr::select(gene, uniprotid, description, target_id, entrezid, logFC, AveExpr, P.Value, adj.P.Val) %>%
  mutate(entrezid = sapply(entrezid, "[[", 1))


write.csv(orderedRes, file="../analysis/Metastasis/LimmaResults.csv")

##prepare publication tables----
number = 20

topnL <- resL_sig2 %>% 
  arrange(desc(abs(logFC))) %>%
  head(n=number)

topnL$P.Value <- sprintf('%.2e', topnL$P.Value)
topnL$adj.P.Val<- sprintf('%.2e', topnL$adj.P.Val)

#make the gt table
topnL %>%
  dplyr::select(c(gene,description, AveExpr,logFC, P.Value, adj.P.Val)) %>%
  gt()  %>%
  fmt_number(columns=c(2:5), decimals = 2) %>%
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
  tab_header(title = md(paste0("Top ", number,  " differentially regulated genes")),
             #subtitle = md("**to choose fit type***")
  ) %>%
  gtsave(paste0("../analysis/Metastasis/top", number, "DEGs.pdf"))

