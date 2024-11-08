##Differential Expression Analysis
#https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/thursday/DE.html

# load packages ----
suppressPackageStartupMessages({  
  library(edgeR)
  library(gt)
  library(DEFormats)
  library(Glimma)
  library(DT)
  library(dplyr)
  library(tibble)
  library(org.Hs.eg.db)
  library(EnhancedVolcano)
  })


colData(ddsIn) <- cbind(colData(ddsIn), pdat[,1:5])
myDGEListIn <- as.DGEList(ddsIn)
myDGEListIn.norm <- calcNormFactors(myDGEListIn, method = "TMM")
myDGEListIn.norm$samples <- droplevels(myDGEListIn.norm$samples)

#d1 <- estimateCommonDisp(d, verbose=T)
#d1 <- estimateTagwiseDisp(d1)
#plotBCV(d1)


group <- myDGEListIn.norm$samples$group
disease <- myDGEListIn.norm$samples$Patient_group
#hospital <- myDGEListIn.norm$samples$Hospital
W_1 <- myDGEListIn.norm$samples$W_1
W_2 <- myDGEListIn.norm$samples$W_2
W_3 <- myDGEListIn.norm$samples$W_3
W_4 <- myDGEListIn.norm$samples$W_4
W_5 <- myDGEListIn.norm$samples$W_5
age <- as.numeric(myDGEListIn.norm$samples$Age)
storageT <- myDGEListIn.norm$samples$Storage_time
levels(storageT) <- c("long", "short")
# gender <-  d$samples$Gender


design <- model.matrix(~0 + group + W_1 + W_2 + W_3 + W_4 + W_5)
colnames(design)[1:2] <- c("noncancer", "NSCLC")

# colnames(design)[1:9] <- c("ChronicPanc", "Epilepsy",
#                            "Healthy", "MS", "NSathero",
#                            "NSCLC", "PulmonaryHyperT",
#                            "StableAP", "UnstableAP")


#tmp <- voom(myDGEListColl, mm, plot = T) #check how unfiltered data looks

v.myDGEListIn.norm <- voom(myDGEListIn.norm, design, plot = TRUE)

#fit linear model
fit <- lmFit(v.myDGEListIn.norm, design)
head(fit$coefficients)

# Contrast matrix ----
contrast.matrix <- makeContrasts(NSCLC - noncancer, levels = design) #limma

# extract the linear model fit -----
fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=0.58)

plotMD(ebFit, status=results)


#tmp2 <- limma::treat(tmp)
#top.table2 <- topTreat(tmp2, adjust ="BH", sort.by = "P", n = Inf)

top.table <- topTable(ebFit, adjust ="BH", sort.by = "P", n = Inf)

resL_tb <- top.table %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

#quantile(resL_tb$AveExpr, probs = seq(0,1,0.1))
#        0%        10%        20%        30%        40%        50%        60% 
#-0.9725104  0.1048774  0.7581599  1.2964740  1.8638633  2.4042587  2.9906460 
#70%        80%        90%       100% 
#  3.6196317  4.3745322  5.5070536 15.7778974 

orderedResL <- resL_tb %>% arrange(adj.P.Val)
orderedResL <- merge(orderedResL, annotations, by='gene', all.x = T) %>% 
  distinct(., gene, .keep_all = T)
orderedResL <- merge(orderedResL, TxOG[6:8], by='gene', all.x = T) %>% 
  distinct(., gene, .keep_all = T)
orderedResL$uniprotid <- mapIds(org.Hs.eg.db, keys = orderedResL$gene, keytype = "SYMBOL", column = "UNIPROT", multiVals = "first")
  

orderedResL$description <- unlist(lapply(strsplit(as.character(orderedResL$description), "\\["), '[[', 1))
orderedResL <- orderedResL %>% dplyr::select(gene, uniprotid,description, target_id, entrezid, logFC, AveExpr, P.Value, adj.P.Val)
orderedResL$entrezid <-  unlist(sapply(orderedResL$entrezid,"[[",1))

write.csv(orderedResL, file="../data/processed_data/LimmaResults.csv")

#histograms ----
hist(resL_tb$logFC, 
     col = "white",
     breaks = 100,
     #breaks=0:50/50, 
     xlab="log2FC", main="Histogram of log2FCs")
abline(v = c(-0.5, 0.5), col = "#973d21", lwd = 1.5, lty = 'dashed')


hist(resL_tb$P.Value, 
     col = "white",
     breaks=0:50/50,
     xlab="p value", main="Histogram of p values")
abline(v = c(0.05, 0.01), col = c("#e78429","#ab3329"), lwd = 1.5, lty = 'dashed')
text(x=0.025, y=5800, '0.01', col="#ab3329", srt = 270, cex = 0.8)
text(x=0.065, y=5800, '0.05', col="#e78429", srt = 270, cex = 0.8)



##get differentially expressed genes----

results <- decideTests(ebFit, method="global", adjust.method="BH", p.value=0.05, lfc=0.58)
#colnames(v.myDGEListIn.norm$E) <- sampleLabels
diffGenes <- v.myDGEListIn.norm$E[results[,1] !=0,]
dim(diffGenes)

#convert your DEGs to a dataframe using as_tibble
diffGenes.df <- as_tibble(diffGenes, rownames = "geneID")

# create interactive tables to display your DEGs ----
datatable(diffGenes.df,
          extensions = c('KeyTable', "FixedHeader"),
          caption = 'Table 1: DEGs in NSCLC',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:11), digits=2)

#get sig dfs----
padj.cutoff <- 0.05
lfc.cutoff <- 0.58

resL_sig <- resL_tb %>%
  dplyr::filter(`adj.P.Val` < padj.cutoff & abs(logFC) > lfc.cutoff)

resL_sig0 <- resL_tb %>%
  dplyr::filter(`adj.P.Val` < padj.cutoff)

#length(which(top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > 0.58))

##make the top table----

#add descriptions
resL_sig2 <- merge(resL_sig, annotations, by='gene', all.x = T) %>% 
  distinct(., gene, .keep_all = T)

resL_sig2$description <- unlist(lapply(strsplit(as.character(resL_sig2$description), "\\["), '[[', 1))

##prepare publication tables----
number = 20

topnL <- resL_sig2 %>% 
  arrange(desc(abs(logFC))) %>%
  head(n=number)

topnL$P.Value <- sprintf('%.2e', topnL$P.Value)
topnL$adj.P.Val<- sprintf('%.2e', topnL$adj.P.Val)

topnL <- topnL %>% dplyr::select(c(gene,description, AveExpr,logFC, P.Value, adj.P.Val))


#make the gt table
topnL %>%
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
  gtsave(paste0("../analysis/diffExpressionlimma/top", number, "DEGs.pdf"))





#compare DESeq and Limma----
library(ggVennDiagram)

ggVennDiagram(list(topnD$gene, topnL$gene),
              category.names = c("DESeq2", "Limma voom"),
              label_alpha = 0,
              set_size = 4,
              set_color = "#5F5F5F",
              edge_size = 0)+
  scale_fill_gradient(low = "#f6f2ee", high =  "#E4B09A")+
  #scale_x_continuous(expand = expansion(mult = .3))+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12, colour = "#5F5F5F", face = "bold"),
        plot.subtitle  = element_text(hjust = 0.5, colour = "#5F5F5F"))+
  labs(title = "Top 20 genes",
       subtitle = "assessed by two methods")+
  coord_flip()

v1 <- ggVennDiagram(list(resD_sig$gene, resL_sig$gene),
              category.names = c("DESeq2", "Limma voom"),
              label_alpha = 0,
              set_size = 4,
              set_color = "#5F5F5F",
              edge_size = 0) +
  scale_fill_gradient(low = "#f6f2ee", high =  "#87bcbd")+
  scale_x_continuous(expand = expansion(mult = .3))+
  #coord_flip()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12, colour = "#4E696A", face = "bold"),
        plot.caption = element_text(hjust = 0.5, colour = "#5F5F5F"))+
  labs(title = "Differentially expressed genes",
       caption = "|log2FC| > 0.58 & adj.p-value < 0.05")


v2 <- ggVennDiagram(list(resD_sig0$gene, resL_sig0$gene),
              category.names = c("DESeq2", "Limma voom"),
              label_alpha = 0,
              set_size = 4,
              set_color = "#5F5F5F",
              edge_size = 0) +
  scale_fill_gradient(low = "#f6f2ef", high =  "#7c4b73")+
  scale_x_continuous(expand = expansion(mult = .3))+
  #coord_flip()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 12, colour = "#52334C", face = "bold"),
        plot.caption = element_text(hjust = 0.5, colour = "#5F5F5F"))+
  labs(title = "Significantly expressed genes",
       caption = "adj.p-value < 0.05")

ggarrange(v1,v2)

ggsave("DESeqvsLimmaVenns.tiff", 
       device = "tiff",
       path = "../analysis/",
       width = 8, height = 6, bg = "white")

#make Volcano plots-------
top.x = resL_tb
top.x.b <- top.x[abs(top.x$logFC) > 0.58,]
#top.x.b <- top.x.b[order(top.x.b$adj.P.Val),]
top.x.b <- top.x.b[order(abs(top.x.b$logFC), decreasing = T),]
GeneSymbol <- top.x$gene
GeneSymbol2 <- top.x.b$gene

keyvals <- ifelse(
  top.x$adj.P.Val < 0.05 & top.x$logFC < -0.58, 'royalblue',
  ifelse(top.x$adj.P.Val < 0.05 & top.x$logFC > 0.58, 'orangered',
         ifelse(top.x$adj.P.Val < 0.05, 'grey', 'black'))
)

names(keyvals)[keyvals == 'black'] <- 'Not.sig.'
names(keyvals)[keyvals == 'orangered'] <- 'Upregulated'
names(keyvals)[keyvals == "grey"] <- 'Sig.'
names(keyvals)[keyvals == 'royalblue'] <- 'Downregulated'

tiff(file="../analysis/diffExpressionlimma/VolcanoPlot3.tiff", 
     width = 16 * 300, height = 12 * 300, res = 300, compression = "lzw")

print(EnhancedVolcano(top.x,
                      lab = GeneSymbol,
                      selectLab = GeneSymbol2[1:30],
                      #col = c("grey30",  "#7db0ea", "#89ab7c", "#a00e00"),
                      colCustom = keyvals,
                      x = 'logFC',
                      y = 'adj.P.Val',
                      xlab = bquote(~Log[2]~ 'fold change'),
                      pCutoff = 0.05,
                      FCcutoff = 0.58,
                      cutoffLineType = 'twodash',
                      cutoffLineWidth = 0.8,
                      pointSize = 6.0,
                      labSize = 8.0,
                      colAlpha = 0.5,
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
                      ylim = c(0, -log10(min(top.x$adj.P.Val))*1.1),
                      title = "NSCLC vs. nonCancer"
))
dev.off()


# Results of adding age and storage time to design ----
#design <- model.matrix(~0 + group + age + storageT+W_1 + W_2)
resL_tbN <- read.csv("../data/processed_data/LimmaResults2.csv")
resL_sigN <- resL_tbN %>%
  dplyr::filter(`adj.P.Val` < padj.cutoff & abs(logFC) > lfc.cutoff)


