

##Does the control centrality of critical, indispensable, and critical&indispensable nodes differ between non-ICU and ICU patients?

#### load packages ----
suppressPackageStartupMessages({ 
  library("readxl")
  library(dplyr)
  library("ggsci")
  library("ggpubr")
  library("scales")
  library(psych)
  library(cowplot)
  })



#function
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}   #http://www.sthda.com/english/wiki/ggplot2-barplots-quick-start-guide-r-software-and-data-visualization



#### import data ####
data1 <- indis %>% 
  mutate(node = rep("indispensable", nrow(indis)))

data2 <- critical  %>%
  mutate(node = rep("critical", nrow(critical)))

data <- rbind(data1,data2)
data$node.f2 <- factor(data$node, levels=c("critical", "indispensable"))

colnames(data)

features <- colnames(data %>% dplyr::select(where(is.numeric)))

#### check data, describe, summarize ####
for (i in seq_along(features)){
  var = features[i]
  df_sum<- data_summary(data, varname=var, 
                        groupnames=c("node.f2"))
  ##### compare means with Mann-Whitney U test ####
  stat.test <- compare_means(
    as.formula(paste(var, "~ node.f2")), data = data, 
    method = "wilcox.test", p.adjust.method = "BH"
  )
  
  #### make the bar plot ####
  p <- ggbarplot(data, x = "node.f2", y = var, 
                 ylab= var, 
                 xlab = "",
                 fill = "node.f2", 
                 width = 0.5,
                 add = c("mean_se","point"),
                 add.params = list(alpha = 0.2),
                 alpha = .5)
  
  fig2 <- p + stat_compare_means(aes(label = "p.format"),
                                 label.y.npc = "top",
                                 #label.y = 1.4,
                                 label.x.npc = "middle",
                                 hjust = 0.5,
                                 method = "wilcox.test")+
    scale_fill_grey(start=0, end=0.8)+
    theme_bw()+
    theme(legend.position = "none", axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 14))              
  
  fig2
  #### export to high resolution svg images ####
  ggsave(filename = paste0("../analysis/6_NetworkCons/PlateletNetwork/Barplots/", var, ".tiff"),
         device = "tiff",
         width = 5,
         height = 5)
  
}


#### check normality of groups to compare ####
qqnorm(data[data$node.f2=="critical",][[var]])
qqline(data[data$node.f2=="critical",][[var]])
hist(data[data$node.f2=="critical",][[var]])

qqnorm(data[data$node.f2=="indispensable",][[var]])
qqline(data[data$node.f2=="indispensable",][[var]])
hist(data[data$node.f2=="indispensable",][[var]])


