
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
})



# Examine your data up to this point ----
myTPM <- Txi_gene$abundance
myCounts <- Txi_gene$counts

colSums(myTPM)
colSums(myCounts)

colnames(myCounts) <- phenodataF$sampleLabels
colnames(myTPM) <- phenodataF$sampleLabels


# Generate summary stats for your data ----
# 1st, calculate summary stats for each transcript or gene, and add these to your data matrix
# then use the base R function 'transform' to modify the data matrix (equivalent of Excel's '=')
# then we use the 'rowSds', 'rowMeans' and 'rowMedians' functions from the matrixStats package
myTPM.stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM))

# look at what you created
head(myTPM.stats)

TPMplot <- ggplot(myTPM.stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape=16, size=2) +
  geom_smooth(method=lm) +
  geom_hex(show.legend = T) +
  labs(y="Median", x = "Standard deviation",
       title="Transcripts per million (TPM)",
       subtitle="unfiltered, non-normalized data",
       caption="") +
  theme_cowplot()

tiff("../analysis/dataExploration/heteroscedasticity.tiff",
     width = 6*300,
     height = 5*300,
     res = 300, 
     compression = "lzw")
TPMplot
dev.off()

#Collapse technical replicates ----

#create Deseq dataset ----

#remove the one sample from PISA hospital
ddsRaw <- DESeqDataSetFromTximport(Txi_gene, 
                                colData=phenodataF, 
                                design= ~class)
colnames(ddsRaw@assays@data@listData[["counts"]]) <- ddsRaw@colData@listData[["sampleLabels"]]

#ageSubidx <- which(colData(ddsRaw)$Age < 61 & colData(ddsRaw)$Age  > 39)
#ddsRaw <- ddsRaw[, ageSubidx]

#collapse technical replicates
ddsColl <- collapseReplicates(ddsRaw, ddsRaw$geo_accession, ddsRaw$Run)

phenodataColl <- data.frame(ddsColl@colData@listData)
phenodataColl$Run <- phenodataColl$runsCollapsed
phenodataColl <- phenodataColl[1:23]

#remove the one sample from PISA hospital
rmix <- which(rownames(colData(ddsColl)) == "GSM2391029")
ddsColl <- ddsColl[,-rmix]
ddsColl$Hospital <- droplevels(ddsColl$Hospital)

#remove samples with no age assigned

# ddsColl$Age <- as.numeric(ddsColl$Age)
# ageidx <- which(is.na(ddsColl$Age))
# ddsColl <- ddsColl[, -ageidx]
# 
# 

colnames(ddsColl@assays@data@listData[["counts"]]) <-
  ddsColl@colData@listData[["sampleLabels"]]

#assign collapsed counts
myCountsColl <- ddsColl@assays@data@listData[["counts"]]
colnames(myCountsColl) <- ddsColl@colData@listData[["sampleLabels"]]


#Filtering low counts! ----
y <- as.DGEList(ddsColl)
keep <- filterByExpr(y) #design ~class 
table(keep)

y <- y[keep,] #20752 14545
ddsFilt <- ddsColl[keep,]

#plot the count distributions before and after filtering

inputCounts <-cpm(ddsColl@assays@data@listData[["counts"]], log = T)
dd <- reshape2::melt(inputCounts, variable.name = "sample")

densR <- ggplot(dd, aes(value, colour = Var2)) + geom_density() +
  theme_cowplot() +
  xlab("Log-cpm")+
  ggtitle("Raw Data")+
  theme(legend.position="none") + 
  xlim(min(dd$value), max(dd$value))

inputCounts <-cpm(ddsFilt@assays@data@listData[["counts"]], log = T)
dd <- reshape2::melt(inputCounts, variable.name = "sample")

densFilt <- ggplot(dd, aes(value, colour = Var2)) + geom_density() +
  theme_cowplot() +
  xlab("Log-cpm")+
  ggtitle("Filtered Data")+
  theme(legend.position="none") + 
  xlim(min(dd$value), max(dd$value))

ggarrange(densR, densFilt)

ggsave("DensityPlots.tiff",
       device = "tiff",
       path = "../analysis/dataExploration/",
       width = 10, height = 4, bg = "white")

# Make a DGElist from your counts ----

##make the dge list for rawdata (before collapsing) ----

myDGEListRaw = as.DGEList(ddsRaw)
myDGEListColl <- as.DGEList(ddsColl)
myDGEListFilt <- as.DGEList(ddsFilt)

# The names argument tells the barplot to use the sample names on the x-axis
# The las argument rotates the axis names
LibSizeDF <- data.frame(Sample = myDGEListColl$samples$Run, 
                        Libsize = myDGEListColl$samples$lib.size,
                        Class = myDGEListColl$samples$group, 
                        Hospital =myDGEListColl$samples$Hospital,
                        Patient = myDGEListColl$samples$Patient_group,
                        StorageT = myDGEListColl$samples$submission_date)  %>%
  arrange(Libsize)  # Order by LibSize


#LibSizeDF$Sample <- factor(LibSizeDF$Sample, levels = LibSizeDF$Sample)
libsizeplot <- ggplot(LibSizeDF, aes(x = Sample, y = Libsize, fill = Patient)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Library Size", title = "Library Size per Sample")+
  theme(axis.text.x = element_blank())+
  scale_fill_met_d("Signac")

ggsave("LibSizes.tiff", 
       libsizeplot, 
       device = "tiff",
       path = "../analysis/dataExploration/",
       width = 7, height = 5, bg = "white")

# take a look at the DGEList object 
myDGEListColl

save(myDGEListColl, file = "../data/raw_data/myDGEListColl")

#load(file = "../data/raw_data/myDGEList")


#Let's explore the data! ----

##Prince plots----
sub = "vstF" #Raw or Coll or Filt
mydds <- ddsIn #ddsRaw or ddsColl
myswamp <- mydds@assays@data@listData[["counts"]]  

#myswamp <- assay(vsd)

set.seed(19)
o1<-data.frame(Class = mydds$class,
               Disease = mydds$Patient_group,
               Hospital = mydds$Hospital,
               StorageTime = mydds$Storage_time, 
               #W1 = vsd$W1,
               #W2 = vsd$W2,
               #W3 = vsd$W3,
               #W4 = vsd$W4,
               #W5 = vsd$W5,
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


tiff(paste0("../analysis/dataExploration/PrincePlot", sub,"Data_UVs.tiff"),
     width = 10 * 300,
     height = 6 * 300,
     res = 300,
     compression = "lzw")
prince.plot(prince=res1, key = T, note = T, breaks = 16, colsep = 1:9,
            col =  rev(colorRampPalette(c("white", "lightgray", "gold", "orange", "red"))(15)),
            #               col =  heat.colors(15)
)
dev.off()

##PCA plots ----

##Following https://github.com/mikelove/preNivolumabOnNivolumab/blob/main/preNivolumabOnNivolumab.knit.md

###variance-stabilizing transformation----
sub = "Filt"
desc = "filtered data (variance stabilizing transformation)"

vsd <- vst(ddsFilt, blind=F)

#if you want to use only top 500 highest variance to make the PC
#rv <- rowVars(assay(vsd))
#pc <- prcomp(t(assay(vsd)[head(order(-rv),500),]))

pc <- prcomp(t(assay(vsd)))

#calculate percentages
pc.var<-pc$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
head(pc.per)

pca.res.df <- as_tibble(pc$x)

class = vsd$class

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
       path = "../analysis/dataExploration/",
       width = 8, height = 6, bg = "white")


ggplotly(myplot) #interactive plot

###find outliers ----

idx <- pc$x[,1] < -190 | pc$x[,2] > 100
sum(idx)

myplot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color = class, shape = idx, text = vsd$Run)+
  geom_point()+
  theme_bw() +
  #geom_label(nudge_y = 0, nudge_x = 0) +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")"))+
  ylab(paste0("PC2 (",pc.per[2],"%",")"))+ 
  labs(title=paste0("PCA plot of ", desc),
       caption=paste0("produced on ", Sys.time()))+
  coord_fixed()

ggsave(paste0("PCA", sub,"Data_Outliers.tiff"), 
       myplot, 
       device = "tiff",
       path = "../analysis/dataExploration/",
       width = 8, height = 6, bg = "white")


###remove outliers and redo vst ----

class <- vsd$class
class <- class[!idx]
ddsFiltF <- ddsFilt[,!idx]
colData(ddsFiltF) <- droplevels(colData(ddsFiltF))

sub = "vstF"
desc = "filtered data (vst, outliers removed)"
vsd <- vst(ddsFiltF, blind=F)
pc <- prcomp(t(assay(vsd)))


#calculate percentages
pc.var<-pc$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
head(pc.per)

pca.res.df <- as_tibble(pc$x)


mycolors <- c("#92C051", "#633372","#de597c","#FFC0CB","#fe9b00",
              "gray","#fbe183","#5192C0","#51C092")


avalues= c(.8,.8,.8,.8,.8,0.3,.8,.8,.8)
fillCOLS = sapply(1:9,function(i)alpha(mycolors[i],avalues[i]))

labs <- levels(disease)

desc = "Filtered Data (known variation)"
myplot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, fill = disease, color = disease, text = vsd$Run)+
  geom_point(shape = 21)+
  theme_bw() +
  #geom_label(nudge_y = 0, nudge_x = 0) +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")"))+
  ylab(paste0("PC2 (",pc.per[2],"%",")"))+ 
  labs(title=paste0("PCA plot of ", desc),
       caption=paste0("produced on ", Sys.time()))+
  coord_fixed()+
  scale_color_manual(values = mycolors, labels = labs)+
  scale_fill_manual(values = fillCOLS, labels = labs)+
  theme(legend.position = "bottom", legend.title = element_blank(),
        legend.justification.bottom = c(0.5,1),
        legend.key.size = unit(0.7, 'lines'),
        legend.text = element_text(size =11))+
  guides(colour = guide_legend(nrow = 3))
  #scale_color_met_d("Cross")

ggsave(paste0("PCA", sub,"DataDisease.tiff"), 
       myplot, 
       device = "tiff",
       path = "../analysis/dataExploration/",
       width = 6, height = 4, bg = "white")

ggplotly(myplot)


#pc plot

pca.res.df <- pc$x[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = vsd$sampleLabels,
             group = class, 
             hospital = vsd$Hospital,
             W5 = vsd$W5,
             disease = vsd$Patient_group)

pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of new variable (column) storing all the values (data)

ggplot(pca.pivot) +
  aes(x=sample, y=loadings, fill=disease) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_cowplot() +
  coord_flip()+
  theme(axis.text.y = element_blank())


ggsave(paste0("Top4PCs", sub,"Data_disease.tiff"), 
       device = "tiff",
       path = "../analysis/dataExploration/",
       width = 8, height = 6, bg = "white")




###estimate RUV factors + technical variation----
#https://github.com/mikelove/preNivolumabOnNivolumab/blob/main/preNivolumabOnNivolumab.knit.md

#First we see what happens when we run simple DESeq2 analysis comparing the two groups, 
#without attempting to control for the technical variation. Here we use glmGamPoi which 
#is an efficient method for estimating dispersion when we have many samples. This can 
#speed up the analysis by an order of magnitude. For details, see: https://doi.org/10.1093/bioinformatics/btaa1009

ddsIn <- ddsFiltF
ddsIn <- DESeq(ddsIn, test="LRT", reduced=~1, fitType="glmGamPoi")


res <- results(ddsIn)
table(res$padj < .05)
table(abs(res$log2FoldChange) < .05)

table(abs(res$log2FoldChange) < .05 & res$padj > .05) #empiricals

#FALSE  TRUE 
#7186  7359 

#There are thousands of DE genes. As we will see, these are the result of uncontrolled confounding.
#Ignoring the technical variation is not appropriate, if there is correlation 
#between the technical variation and the condition. 

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

#display.brewer.all(type="qual")

myplot <- ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color = vsd$class, text = vsd$Run)+
  geom_point(alpha = 0.9)+
  theme_bw() +
  #geom_label(nudge_y = 0, nudge_x = 0) +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc.per[1],"%",")"))+
  ylab(paste0("PC2 (",pc.per[2],"%",")"))+ 
  #labs(title=desc,
  #    caption=paste0("produced on ", Sys.time()))+
  guides(color=guide_legend(title="class"))+
  coord_fixed()
  

ggsave(paste0("PCA", sub,"Data.tiff"), 
       myplot, 
       device = "tiff",
       path = "../analysis/dataExploration/",
       width = 8, height = 6, bg = "white")


##visualize the unwanted variation

Archambault = c("#88a0dc", "#381a61", "#7c4b73",
                "#ed968c", "#ab3329", "#e78429", "#f9d14a")

Johnson = c("#a00e00", "#d04e00", "#f6c200", "#0086a8", "#132b69")

OKeeffe1 = c("#6b200c", "#973d21", "#da6c42", "#ee956a", "#fbc2a9",
             "#f6f2ee", "#bad6f9", "#7db0ea", "#447fdd", "#225bb2", "#133e7e")

VanGogh1 = c("#2c2d54", "#434475", "#6b6ca3", "#969bc7", "#87bcbd",
                  "#89ab7c", "#6f9954")

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
       path = "../analysis/dataExploration/",
       width = 10, height = 7, bg = "white")



