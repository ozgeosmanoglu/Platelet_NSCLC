# load libraries ----
suppressPackageStartupMessages({
  library(GEOquery)
  library(readr)
  library(tidyverse)
})

#GSE89843<- getGEO('GSE89843',GSEMatrix=T)
#phenodata_geo <- GSE89843[["GSE89843_series_matrix.txt.gz"]]@phenoData@data


#phenodata from geo and pub - combine ----
gse <- getGEO(filename="../data/meta_data/GSE89843_series_matrix.txt.gz", AnnotGPL=TRUE)
phenodata_geo <- gse@phenoData@data

phenodata_geo$'GEO ID (GSE89843)' <- gsub("Blood_Platelets_", "",phenodata_geo$title)
phenodata_geo$'GEO ID (GSE89843)' <- gsub("-Healthy Control-", "-HC-", phenodata_geo$'GEO ID (GSE89843)')

#phenodata from the publication Best et al., 2017
phenodata_pub <- readr::read_delim("../data/meta_data/GSE89843_Best_et_al_Cancer_Cell_2017_Table_S1_conversion_to_GEO.txt") 

#phenodata from ENA of bioproject PRJNA353588----
phenodata_ena <- read_csv(file = "../data/meta_data/SraRunTable.txt")
#find what features are needed

inputPhenos <- c("phenodata_ena", "phenodata_geo", "phenodata_pub")

importantFeatures <- list()
for (j in 1:length(inputPhenos)) {
  inputPheno <- get(inputPhenos[j])
  for (i in 1:ncol(inputPheno)) {
    #print(dim(table(inputPheno[i])))
    if (between(dim(table(inputPheno[i])), 2, 100)) {
      importantFeatures[[inputPhenos[j]]][i] <- colnames(inputPheno[i])
      importantFeatures <- lapply(importantFeatures, na.omit)
      
    }
  }
}

rm(i)
rm(j)

c(importantFeatures$phenodata_ena)
c(importantFeatures$phenodata_geo)
c(importantFeatures$phenodata_pub)

##

phenodataF <- merge(phenodata_geo, phenodata_pub, by = 'GEO ID (GSE89843)')
phenodataF <- phenodataF %>% dplyr::select(c(`GEO ID (GSE89843)`, geo_accession,
                                             `Sample name`, `Patient group`, submission_date,
                                             `Classification group`, `Storage time`, Hospital,
                                             Age, Gender, Smoking, Metastasis,
                                             Treatment, `Matched cohort`, `Full cohort`))


#phenodata_ena$Rungz <- paste0(phenodata_ena$Run, ".fastq.gz")

duplicates <- phenodata_ena[duplicated(phenodata_ena$`Sample Name`) | duplicated(phenodata_ena$`Sample Name`, fromLast=TRUE),]
duplist <- duplicates$Run

#RunN <- aggregate(Run ~ Sample.Name, phenodata_ena[phenodata_ena$Run %in% duplist,], paste, collapse="_")
#duplist_N <- merge(RunN, phenodata_ena, by = "Sample.Name")

#duplist <- paste0(duplicates$Run, ".fastq.gz")
# commandtocombine <- paste0("cat ", aggregate(Rungz ~ Sample.Name, phenodata_ena[phenodata_ena$Rungz %in% duplist,], paste, collapse=" ")$Rungz,
#        " > ", aggregate(Run ~ Sample.Name, phenodata_ena[phenodata_ena$Rungz %in% duplist,], paste, collapse="_")$Run, ".fastq.gz")
# 
# write.table(commandtocombine, "commandstocombine.txt", row.names = F, col.names = F)

agg.dup <- aggregate(. ~`Sample Name`, duplicates, paste, collapse = "_", na.action=NULL)

phenodata_ena <- phenodata_ena[!phenodata_ena$Run %in% duplist,]
phenodata_ena <-  phenodata_ena %>% mutate(across(everything(), as.character))
phenodata_ena <- data.frame(rbind(phenodata_ena, agg.dup))
names(phenodata_ena)[names(phenodata_ena) == 'GEO_Accession..exp.'] <- "geo_accession"
phenodata_ena<- phenodata_ena %>% dplyr::select(Run,BioSample,Experiment,geo_accession, AvgSpotLen,
                                                  DATASTORE.filetype, DATASTORE.provider, DATASTORE.region)
phenodata_ena <- phenodata_ena %>% separate(BioSample, c('Biosample', NA), sep = "_",
                                  remove = T)
phenodata_ena <- phenodata_ena %>% separate(Experiment, c('Experiment', NA), sep = "_",
                                                          remove = T)
phenodata_ena <- phenodata_ena %>% separate(geo_accession, c('geo_accession', NA), sep = "_",
                                                          remove = T)    
phenodata_ena <- phenodata_ena %>% separate(DATASTORE.filetype, c('filetype', NA), sep = "_",
                                              remove = T)    
phenodata_ena <- phenodata_ena %>% separate(DATASTORE.provider, c('provider', NA), sep = "_",
                                              remove = T)   
phenodata_ena <- phenodata_ena %>% separate(DATASTORE.region, c('region', NA), sep = "_",
                                              remove = T)   
phenodata_ena <- phenodata_ena %>% separate(AvgSpotLen, c('AvgSpotLen', NA), sep = "_",
                                              remove = T)   

#combine all info ----


phenodataF <- merge(phenodataF, phenodata_ena, by = "geo_accession")
phenodataF <- separate_rows(phenodataF, Run, sep = "_")

#phenodataMsub <- phenodataM[grep("_", phenodataM$Run),]
#phenodataF <- rbind(phenodataM, phenodataMsub)
#edit out special characters
colnames(phenodataF) <- gsub(" ","_",colnames(phenodataF))
colnames(phenodataF) <- gsub("\\(|\\)","",colnames(phenodataF))
phenodataF <- data.frame(sapply(phenodataF,function(x) {gsub(" ", "_", x)}))
phenodataF <- data.frame(sapply(phenodataF,function(x) {gsub(">", "longerT", x)}))
phenodataF <- data.frame(sapply(phenodataF,function(x) {gsub("<", "shorterT", x)}))
#change some colnames - create new name
colnames(phenodataF)[colnames(phenodataF) == "Classification_group"] <- "class"
phenodataF$sampleLabels <- c(paste0("LC",1:length(grep("NSCLC", phenodataF$class)))
                             ,paste0("HD", 1:length(grep("Non-cancer", phenodataF$class))))

phenodataF[sapply(phenodataF, is.character)] <- 
  lapply(phenodataF[sapply(phenodataF, is.character)],  as.factor)

phenodataF <- phenodataF %>% mutate_if(is.character,as.factor)
