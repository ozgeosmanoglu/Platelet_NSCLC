# Load packages ----
suppressPackageStartupMessages({ 
  library(SANTA)
  library(igraph)
  library(dplyr)
  library(ggplot2)
  library(ggraph)
  library(gt)
})

# 2. DRUG REPURPOSING----

## load drugs data ----
drugbank <- read.csv("../analysis/7_DrugRepurposing/uniprot links.csv")
drug_targets_full <- drugbank %>% dplyr::select(Name, UniProt.ID)

subnetworkNC <- read.csv("../analysis/6_NetworkCons/PlateletNetwork/Subnetwork.txt default node.csv")

indisinsub <- subnetworkNC %>%
  filter(Classification == "2") %>%
  pull(shared.name)


## proteins of interest ----
drug_targets <- drug_targets_full %>%
  filter(UniProt.ID %in% indisinsub) 

drugCounts <- drug_targets %>% 
  dplyr::count(Name, sort = TRUE)

top5Drugs <- drug_targets %>% 
  dplyr::count(Name, sort = TRUE) %>% 
  head(n=5) %>%
  pull(Name)

DTT <- drug_targets %>%
  dplyr::left_join(PlNetworkNodeAnnots, by = c("UniProt.ID" = "uniprot")) %>%
  dplyr::group_by(Name) %>%
  dplyr::summarise(
    target_count = n(),
    targets = paste(symbol, collapse = ", ")
  ) %>%
  dplyr::arrange(desc(target_count))

DTT %>%
  gt() %>%
  cols_label(Name = md("**Drug Name**"),
             target_count = md("**# of Targets**"),
             targets = md("**Targets**")) %>%
  cols_align("center",
             columns = 2) %>%
  tab_options(data_row.padding = px(1)) %>%
  tab_footnote(md("FDA approved drugs from Drugbank v5.1.12."), 
               locations = cells_column_labels(columns = Name),
               placement = "right") %>%
  tab_footnote(md("Only reanalyzed indispensable nodes are targeted."), 
               locations = cells_column_labels(columns = targets),
               placement = "right") %>%
  tab_options(data_row.padding = px(0.3)) %>%
  tab_options(table.font.size = 14) %>%
  cols_width(Name ~ px(350))%>%
  gtsave("../analysis/7_DrugRepurposing/SubIndisDrugCounts.docx")


# tables of targets

TTT <- drug_targets %>% 
  dplyr::left_join(PlNetworkNodeAnnots, by = c("UniProt.ID" = "uniprot")) %>%
  dplyr::group_by(symbol) %>%
  dplyr::summarise(
    drug_count = n(),
    drugs = paste(Name, collapse = ", ")
  ) %>%
  dplyr::arrange(desc(drug_count))

TTT %>%
  gt()%>%
  cols_label(symbol = md("**Target Name**"),
             drug_count = md("**# of Drugs**"),
             drugs = md("**Drugs**")) %>%
  cols_align("center",
             columns = 2) %>%
  tab_options(data_row.padding = px(1)) %>%
  tab_footnote(md("FDA approved drugs from Drugbank v5.1.12."), 
               locations = cells_column_labels(columns = symbol),
               placement = "right") %>%
  tab_footnote(md("Only reanalyzed indispensable nodes are targeted."), 
               locations = cells_column_labels(columns = drugs),
               placement = "right") %>%
  tab_options(data_row.padding = px(0.3)) %>%
  tab_options(table.font.size = 12) %>%
  cols_width(drugs ~ px(500))%>%
  gtsave("../analysis/7_DrugRepurposing/SubIndisTargetCounts.docx")


# add drugs to the network

NodesCritIndis <- as_data_frame(critIndisNetwork, "v")
dfDrugs <- data.frame("name"= DTT$Name, "symbol" = DTT$Name, "Degree" = NA,
                      "logFC" = NA, "absLog2FC" = NA, "Weight" = NA, "node_type" = "drug")
NewNodesCritIndis <- rbind(NodesCritIndis,dfDrugs)

EdgesCritIndis <- as_data_frame(critIndisNetwork, "e")
df <-  drug_targets
df <- df %>%
  mutate(source_genesymbol = NA,
         target_genesymbol = NewNodesCritIndis$symbol[match(df$UniProt.ID, NewNodesCritIndis$name)])

dfDrugsInts <- cbind(df, as.data.frame(matrix(NA, nrow = 94, ncol = 17)))
colnames(dfDrugsInts) <- colnames(EdgesCritIndis)
NewEdgesCritIndis <- rbind(EdgesCritIndis,dfDrugsInts)

CritIndiswithDrugs <- graph_from_data_frame(d = NewEdgesCritIndis, 
                                            vertices = NewNodesCritIndis,
                                            directed = T)


V(CritIndiswithDrugs)$node_type = ifelse(V(CritIndiswithDrugs)$name %in% setdiff(indisinsub,unique(drug_targets$UniProt.ID)), "indispensable",
                                     ifelse(V(CritIndiswithDrugs)$name %in% criticalsInNet, "critical",
                                            ifelse(V(CritIndiswithDrugs)$name %in% unique(drug_targets$UniProt.ID),"direct drug target (indispensable)",
                                                   ifelse(V(CritIndiswithDrugs)$name %in% DTT$Name,"drug",
                                                   "intermediate node"))))
  
  
set.seed(2024)
ggraph(CritIndiswithDrugs, layout = "dh") +
  geom_edge_fan(
    aes(edge_colour = type,
        start_cap = label_rect(node1.name),
        end_cap = label_rect(node2.name)),
    arrow = arrow(length = unit(1, 'mm')),
    edge_width = .2,
    edge_alpha = .5
  ) +
  geom_node_point(aes(color = node_type, size = ifelse(node_type == "drug", 5, 8)), alpha = 0.5) +
  geom_node_text(
    aes(label = symbol),
    size = 3,
    color = "gray10",
    position = "identity",
    repel = FALSE
  ) +
  scale_color_manual(values = c("critical"="#A0CFEC",
                                "direct drug target (indispensable)"= "orange",
                                "indispensable" = "#F5E216", 
                                "intermediate node" = "gray30",
                                "drug" = "purple"),
                     guide = guide_legend(title = 'Node type')) +
  scale_edge_color_manual(values = c("gray70", "tomato")) +
  scale_size_continuous(range = c(5, 8), guide = "none") +  # Remove the size legend
  theme_graph() +
  xlab("") +
  ylab("")


ggsave(filename = "../analysis/7_DrugRepurposing/SubIndisDrugNetwork.tiff",
       device = "tiff",
       width = 20,
       height = 12)
     
  




