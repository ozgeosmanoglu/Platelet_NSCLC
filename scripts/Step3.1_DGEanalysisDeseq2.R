#############################################################
###### DGE Analysis with Deseq2 and limma-voom & EdgeR ######
#############################################################
# load packages ----
suppressPackageStartupMessages({  
  library("apeglm")
  library(DESeq2)
  library(gt)
  library(biomaRt)
  library(org.Hs.eg.db)
  library(tibble)
  library(dplyr)
  })


#add variation factors to design and do DGE analysis----
colData(ddsIn) <- cbind(colData(ddsIn), pdat[,1:5])
design(ddsIn) <- ~W_1 + W_2 + W_3 + W_4 + W_5 + class


#check dispersion to figure out the fit type----
#following https://support.bioconductor.org/p/81094/

DDS <- estimateSizeFactors(ddsIn)
par <- estimateDispersions(DDS, fitType = "parametric")
loc <- estimateDispersions(DDS, fitType = "local")
glmgam <- estimateDispersions(DDS, fitType = "glmGamPoi")


tiff("../analysis/dataExploration/DeseqDispEstimates.tiff",
     width = 8 * 300,
     height = 3 * 300,
     res = 300,
     compression = "lzw")
par(mfrow =c(1,3))
plotDispEsts(glmgam, main= "dispEst: glmGamPoi")
plotDispEsts(par, main= "dispEst: parametric")
plotDispEsts(loc, main= "dispEst: local")
dev.off()

residualg <- log(mcols(glmgam)$dispGeneEst) - log(mcols(glmgam)$dispFit)
residualp <- log(mcols(par)$dispGeneEst) - log(mcols(par)$dispFit)
residuall <- log(mcols(loc)$dispGeneEst) - log(mcols(loc)$dispFit)

residualsT <- tibble("fit type" = c("glmGamPoi", "parametric", "local"),
          "median absolute residual" =c(median(abs(residualg), na.rm = T),
                                      median(abs(residualp), na.rm = T),
                                      median(abs(residuall), na.rm = T)))

residualsT %>%
  gt() %>%
  fmt_number(columns=2, decimals = 3) %>%
  tab_header(title = md("**Dispersion Estimates**"),
             subtitle = md("**to choose fit type***")) %>%
  tab_source_note(
    source_note = md("*the type of fitting of dispersions to the mean intensity."))




# g: 0.3823132 #this looks the best, but not possible to lfcshrink
# p: 0.4509929
# l: 0.3972261 #next best


#run with selected fit type----

dds <- DESeq(ddsIn, test="LRT", reduced=~W_1 + W_2 + W_3 + W_4 + W_5,
             fitType="glmGamPoi")


resultsNames(dds) 

resD<- results(dds, name ="classNSCLC")
summary(resD)
plotMA(resD,ylim=c(-2,2))
#mcols(resD, use.names=T)
table(abs(resD$log2FoldChange) > 0.58 & resD$padj < .05)


#logFC shrinkage---- optional, only works with default fittypes ("local", "parametric" etc.)
#Shrinkage of effect size (LFC estimates) is useful for visualization and ranking of genes. 
#To shrink the LFC, we pass the dds object to the function lfcShrink. Below we specify to use 
#the apeglm method for effect size shrinkage (Zhu, Ibrahim, and Love 2018), which improves 
#on the previous estimator (https://bioconductor.statistik.tu-dortmund.de/packages/3.12/bioc/vignettes/DESeq2/inst/doc/DESeq2.html).

#corrects the log2 fold changes for genes with high dispersion. 
#This does not change the p-values or the list of DE genes.

#resDlocalShrunken <- lfcShrink(dds, coef=7, res=resDlocal, type="apeglm", lfcThreshold = 0)
#summary(resLFC)
#plotMA(resDlocalShrunken, ylim=c(-2,2))



#Get results in tibbles ----

resD_tb <- resD %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

##get descriptions for the genes----

mart <- useEnsembl("ensembl","hsapiens_gene_ensembl")

annotations <- getBM(c("description",
             #"uniprotswissprot",
             #"entrezgene_id",
             "external_gene_name"), 
           "external_gene_name", 
           resD_tb$gene, mart)

annotations <- annotations %>% rename_with(~"gene", external_gene_name)

# 
# annotations2 <- merge(annotations, TxOG, by='gene', all = T) %>% 
#   distinct(., gene, .keep_all = T)
# annotations2$entrezid <- unlist(sapply(annotations2$entrezid,"[[",1))
# annotations2 <- annotations2 %>%
#   mutate(across(everything(), as.character))
# 
# annotations2 <- annotations2 %>% 
#   mutate(newcol = case_when(
#     is.na(entrezgene_id) & is.na(entrezid) ~ NA_character_,
#     is.na(entrezgene_id)  & !is.na(entrezid) ~ entrezid,
#     !is.na(entrezid) & is.na(entrezid) ~ entrezgene_id,
#     entrezgene_id == entrezid ~ entrezid,
#   ))




##get differentially expressed genes----
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

resD_sig <- resD_tb %>%
  dplyr::filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

resD_sig0 <- resD_tb %>%
  dplyr::filter(padj < padj.cutoff)

orderedRes <- resD_tb %>% arrange(padj)
write.csv(orderedRes, file="../data/processed_data/Deseq2Results.csv")

normCountsD <- counts(dds, normalized = TRUE)


##histograms ----
hist(resD$log2FoldChange, 
     col = "white",
     breaks = 100,
     #breaks=0:50/50, 
     xlab="log2FC", main="Histogram of log2FCs")
abline(v = c(-0.5, 0.5), col = "#973d21", lwd = 1.5, lty = 'dashed')


hist(resD$pvalue, 
     col = "white",
     breaks=0:50/50,
     xlab="p value", main="Histogram of p values")
abline(v = c(0.05, 0.01), col = c("#e78429","#ab3329"), lwd = 1.5, lty = 'dashed')
text(x=0.025, y=6200, '0.01', col="#ab3329", srt = 270, cex = 0.8)
text(x=0.065, y=6200, '0.05', col="#e78429", srt = 270, cex = 0.8)



#add descriptions
resD_sig2 <- merge(resD_sig, annotations, by='gene', all.x = T) %>% 
  distinct(., gene, .keep_all = T)

resD_sig2$description <- unlist(lapply(strsplit(as.character(resD_sig2$description), "\\["), '[[', 1))

##prepare publication tables----

number = 20

topnD <- resD_sig2 %>% 
  arrange(desc(abs(log2FoldChange))) %>%
  head(n=number)

topnD$pvalue <- sprintf('%.2e', topnD$pvalue)
topnD$padj<- sprintf('%.2e', topnD$padj)

#make the gt table
topnD %>%
  dplyr::select(c(gene,description, baseMean,log2FoldChange, pvalue, padj)) %>%
  gt()  %>%
  fmt_number(columns=c(2:5), decimals = 2) %>%
  cols_width(
    gene ~ px(120),
    pvalue ~ px(200),
    padj ~ px(200),
    description ~ px(1000)) %>%
  cols_label(gene = md("**gene**"),
             baseMean = md("**base mean**"),
             log2FoldChange = md("**log2FC**"),
             pvalue = md("**p-value**"),
             padj = md("**adj. p-value**"),
             description = md("**description**")
             ) %>%
  cols_align("center",
             columns = 3:6) %>%
  tab_header(title = md(paste0("Top ", number,  " differentially regulated genes")),
             #subtitle = md("**to choose fit type***")
             ) %>%
  gtsave(paste0("../analysis/diffExpressionDESeq/top", number, "DEGs.pdf"))




####make count plots----
dds$Patient_group <- gsub("_", " ", dds$Patient_group)
dds$Patient_group <- factor(dds$Patient_group)

#for one gene
gene = "PCSK5"

#for class
d <- plotCounts(dds, gene=gene, intgroup="class", 
                returnData=TRUE)

mean_team <- d %>% group_by(class) %>% 
  summarise(mean_pts=mean(count))

ggplot(d, aes(x=class, 
              y=count, 
              colour = class)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))+
  theme_cowplot() +
  ggtitle(gene)+
  theme(legend.position = "none", axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill="white"),
        strip.placement = "outside",
        strip.clip = "off",
        strip.text = element_text(size = 10, angle = 90, hjust = 1)) +
  facet_grid(~class, 
             shrink = F, scales = "free", space = "free", switch = "x") +
  geom_hline(aes(yintercept = mean_pts, color = class), 
             mean_team, linewidth = 1)

#for Patient_group
d <- plotCounts(dds, gene=gene, intgroup="Patient_group", 
                returnData=TRUE)

mean_team <- d %>% group_by(Patient_group) %>% 
  summarise(mean_pts=mean(count))


ggplot(d, aes(x=Patient_group, 
              y=count, 
              colour = Patient_group)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))+
  theme_cowplot() +
  ggtitle(gene)+
  theme(legend.position = "none", axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        strip.background = element_rect(fill="white"),
        strip.placement = "outside",
        strip.clip = "off",
        strip.text = element_text(size = 10, angle = 90, hjust = 1)) +
  facet_grid(~Patient_group, 
             shrink = F, scales = "free", space = "free", switch = "x") +
  geom_hline(aes(yintercept = mean_pts, color = Patient_group), 
             mean_team, linewidth = 1)


#for multiple #for multiple #for multiple genes

gene_set <- topnD$gene[1:6]
names(gene_set) <- gene_set

for (i in 1:length(gene_set)) {
  gene <- gene_set[i]
  
  d <- plotCounts(dds, gene=gene, intgroup="class", 
                  returnData=TRUE)
  
  mean_team <- d %>% group_by(class) %>% 
    summarise(mean_pts=mean(count))
  plot <- ggplot(d, aes(x=class, 
                        y=count, 
                        colour = class)) + 
    geom_point(position=position_jitter(w=0.1,h=0)) +
    scale_y_log10(breaks=c(25,100,400,1000))+
    theme_cowplot() +
    ggtitle(gene)+
    theme(legend.position = "none", axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          strip.background = element_rect(fill="white"),
          strip.placement = "outside",
          strip.clip = "off",
          strip.text = element_text(size = 12, angle = 90, hjust = 1)) +
    facet_grid(~class, 
               shrink = F, scales = "free", space = "free", switch = "x") +
    geom_hline(aes(yintercept = mean_pts, color = class), 
               mean_team, linewidth = 1)
  
  assign(paste0("plot", i), plot)
  
}
classplot <- ggarrange(plot1, plot2, plot3, plot4, plot5, plot6)
ggsave(paste0("top6PlotCounts_", "class" , ".tiff"), 
       classplot,
       device = "tiff",
       path = "../analysis/diffExpressionDESeq/",
       width = 10, height = 8, bg = "white")

for (i in 1:length(gene_set)) {
  gene <- gene_set[i]
  
  d <- plotCounts(dds, gene=gene, intgroup="Patient_group", 
                  returnData=TRUE)
  
  mean_team <- d %>% group_by(Patient_group) %>% 
    summarise(mean_pts=mean(count))
  plot <- ggplot(d, aes(x=Patient_group, 
                        y=count, 
                        colour = Patient_group)) + 
    geom_point(position=position_jitter(w=0.1,h=0)) +
    scale_y_log10(breaks=c(25,100,400))+
    theme_cowplot() +
    ggtitle(gene)+
    theme(legend.position = "none", axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          strip.background = element_rect(fill="white"),
          strip.placement = "outside",
          strip.clip = "off",
          strip.text = element_text(size = 12, angle = 90, hjust = 1)) +
    facet_grid(~Patient_group, labeller = label_wrap_gen(width=10), 
               shrink = F, scales = "free", space = "free", switch = "x") +
    geom_hline(aes(yintercept = mean_pts, color = Patient_group), 
               mean_team, linewidth = 1)
  
  assign(paste0("plot", i), plot)
  
}
patientplot <- ggarrange(plot1, plot2, plot3, plot4, plot5, plot6)
ggsave(paste0("top6PlotCounts_", "patient" , ".tiff"), 
       patientplot,
       device = "tiff",
       path = "../analysis/diffExpressionDESeq/",
       width = 16, height = 8, bg = "white")


#Volcano plots----

top.x = resD_tb
top.x.b <- top.x[abs(top.x$log2FoldChange) > 0.58,]
top.x.b <- top.x.b[order(top.x.b$padj),]
GeneSymbol <- top.x$gene
GeneSymbol2 <- top.x.b$gene

tiff(file="../analysis/diffExpressionDESeq/VolcanoPlot.tiff", 
     width = 16 * 300, height = 12 * 300, res = 300, compression = "lzw")

print(EnhancedVolcano(top.x,
                      lab = GeneSymbol,
                      selectLab = GeneSymbol2[1:30],
                      col = c("grey30",  "#7db0ea", "#89ab7c", "#a00e00"),
                      x = 'log2FoldChange',
                      y = 'padj',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      cutoffLineType = 'twodash',
                      cutoffLineWidth = 0.8,
                      pointSize = 6.0,
                      labSize = 8.0,
                      colAlpha = 0.3,
                      axisLabSize = 22,
                      legendLabels=c('Not sig.','Not sig. Log2FC','adj.p.value',
                                     'adj.p.value & Log2FC'),
                      legendPosition = 'right',
                      legendLabSize = 22,
                      legendIconSize = 4.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5,
                      maxoverlapsConnectors = 50,
                      boxedLabels = T,
                      xlim = c(-2, 2),
                      ylim = c(0, -log10(min(top.x$padj))*1.1),
                      title = "NSCLC vs. nonCancer"
))
dev.off()

