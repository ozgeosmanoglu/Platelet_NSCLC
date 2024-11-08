#Gene Prioritization following https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9732137/

# Load packages ----
suppressPackageStartupMessages({ 
  library(SANTA)
  library(igraph)
  library(dplyr)
  library(ggplot2)
  library(ggraph)
  library(gt)
  })


# 1. CALCULATE GENE SCORES ----
## Choose network ----

#networkIn <- read.delim("../analysis/6_NetworkCons/PlateletNetwork/PlPPINetwork.txt")
networkIn <- PlNetwork

networkdata <- networkIn %>%
  mutate(
    source_logFC = ifelse(is.na(source_logFC), 10^-6, source_logFC),
    target_logFC = ifelse(is.na(target_logFC), 10^-6, target_logFC),
    weight = 1/abs((source_logFC * target_logFC)) #edge weights
  )
#annotdata <- read.delim("../analysis/6_NetworkCons/PlateletNetwork/PlPPINetworkUniprotSymbol.txt")
annotdata <- PlNetworkNodeAnnots

## Convert to graph
net <- graph_from_data_frame(d = networkdata,
                             directed = T)

## remove nonconnected elements - get largest connected component
components <- igraph::components(net, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(net)[components$membership == biggest_cluster_id]
netSub <- igraph::induced_subgraph(net, vert_ids)

## Calculate node weights ----
degrees <- igraph::degree(netSub, v = V(netSub), mode = c("all"), loops = TRUE, normalized = FALSE)

annotdata <- annotdata %>%
  filter(uniprot %in% V(netSub)$name) %>%
  mutate(
    Degree = degrees[uniprot],
    logFC = resL_sig0$logFC[match(symbol, resL_sig0$gene)],
    absLog2FC = abs(resL_sig0$logFC[match(symbol, resL_sig0$gene)]),
    absLog2FC = replace_na(absLog2FC, 1e-6),
    Weight = Degree * absLog2FC
  )

## Create graph ----

net2 <- graph_from_data_frame(d = as_data_frame(netSub, what= c("edges")), 
                                        vertices = annotdata,
                              directed = T)

#plot(degree_distribution(net, cumulative = FALSE))
#as_data_frame(net2, what="edges")
#as_data_frame(net2, what="vertices")


## Gene Prioritization ----
#ranked by their strength of association with high-weight vertices 

GeneScores <- Knode(net2, dist.method="shortest.paths", 
                    vertex.attr="Weight", 
                    edge.attr = NULL,
                    correct.factor=1, nsteps=1000,verbose=TRUE)
summary(GeneScores)

gsDF <- as.data.frame(GeneScores)
gsDF$uniprot <- rownames(gsDF)
Q90 <- quantile(GeneScores, probs = seq(0, 1, 0.1))[10]

annotdata$GeneScore <- gsDF$GeneScores[match(annotdata$uniprot, gsDF$uniprot)]

write.table(annotdata, "../analysis/6_NetworkCons/PlateletNetwork/PlPPIGeneScores.txt", 
            row.names = F, sep = "\t")

PrioGenes <- annotdata[annotdata$GeneScore > Q90,]
PrioGenes <- PrioGenes %>% arrange(desc(GeneScore))


write.table(PrioGenes, "../analysis/6_NetworkCons/PlateletNetwork/PlPPIPrioGenes.txt", 
            row.names = F, sep = "\t")



QF <- quantile(annotdata$Weight, probs = seq(0, 1, 0.1))[10] #decide on threshold
TopGenes <- annotdata[annotdata$Weight > QF,]
TopGenes <- TopGenes %>% arrange(desc(Weight))

PrioGenes[c(1:4,6:7)] %>%
  gt() %>%
  fmt_number(
    columns = 4:6,
    decimals = 2) %>%
  sub_missing(
    columns = everything(),
    rows = everything(),
    missing_text = "---"
  ) %>%
  cols_label(uniprot = md("**Uniprot ID**"),
             symbol = md("**Name**"),
             Degree = md("**Degree**"),
             logFC = md("**log2FC**"),
             Weight = md("**Weight**"),
             GeneScore = md("**Gene Score**")) %>%
  cols_align("center",
             columns = 3:6)  %>%
  # tab_footnote(md("Weights are calculated by: W(i) = |log2FC(i)| x Degree(i)"), 
  #              locations = cells_column_labels(columns = Weight),
  #              placement = "right") %>%
  tab_options(data_row.padding = px(0.1)) %>%
  tab_options(table.font.size = 11)%>%
  gtsave("../analysis/6_NetworkCons/PlateletNetwork/AllPrioGenes.pdf")



Sources <- PrioGenes$uniprot
Targets <- TopGenes$uniprot

ggVennDiagram(list("gene s." = PrioGenes$uniprot, 
                   "weight" = TopGenes$uniprot, 
                   "critical" = criticalnodes,
                   "indisp." = indisnodes
                   ))+
  scale_fill_gradient(low = "seashell", high = "tomato")
  #coord_flip()

ggsave("../analysis/6_NetworkCons/PlateletNetwork/Venn_genes.tiff",
       device = "tiff",
       width = 6,
       height = 4)



commons <- intersect(TopGenes$uniprot, PrioGenes$uniprot)

collected_path_nodes = list()

for(i in 1:length(Sources)){
  
  paths <- shortest_paths(net2, from = Sources[[i]],
                          to = Targets,
                          weights = NA,
                          output = 'vpath')
  path_nodes <- lapply(paths$vpath,names) %>% unlist() %>% unique()
  collected_path_nodes[[i]] <- path_nodes
}

collected_path_nodes <- unlist(collected_path_nodes) %>% unique()


subnodes <- c(Sources,Targets, collected_path_nodes) %>%
  unique()
subnetwork <- induced_subgraph(graph = net2,vids = subnodes)

## remove nonconnected elements - get largest connected component
components <- igraph::components(subnetwork, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(subnetwork)[components$membership == biggest_cluster_id]
subnetwork <- igraph::induced_subgraph(subnetwork, vert_ids)

#visualize network
V(subnetwork)$node_type = ifelse(V(subnetwork)$name %in% Sources, "prio",
                                       ifelse(V(subnetwork)$name %in% Targets, "high-weight",
                                              "intermediate"))

set.seed(2024)
ggraph(
  subnetwork,
  layout = "dh"
) +
  geom_edge_fan(
    aes(
        start_cap = label_rect(node1.name),
        end_cap = label_rect(node2.name),
    ), color = "gray80",
    arrow = arrow(length = unit(0, 'mm')
    ),
    edge_width = .2,
    edge_alpha = .5
  ) +
  geom_node_point(aes(color = node_type), size = 8, alpha = 0.5) +
  geom_node_text(aes(label = symbol), color = "black",
                 position = "identity",
                 repel = F, size =3) +
  scale_color_manual(values = c("#F5E216","gray30","#A0CFEC"),
                     guide = guide_legend(title = 'Node type')) +
  theme_graph() +
  xlab("") +
  ylab("") 


ggsave(filename = "../analysis/6_NetworkCons/PlateletNetwork/PPIsubnetwork.tiff",
       device = "tiff",
       width = 20,
       height = 10)


set.seed(2024)
ggraph(
  subnetwork,
  layout = "dh"
) +
  geom_edge_fan(
    aes(
        start_cap = label_rect(node1.name),
        end_cap = label_rect(node2.name),
    ), color = "gray80",
    arrow = arrow(length = unit(0, 'mm')),
    edge_width = .2,
    edge_alpha = .5
  ) +
  geom_node_point(aes(color = logFC), size = 8, alpha = 0.7) +  # Color by logFC
  geom_node_text(aes(label = symbol), color = "black",
                 position = "identity",
                 repel = F, size = 3) +
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0, 
                        guide = guide_legend(title = 'logFC')) +  # Gradient scale for logFC
  #scale_edge_color_manual(values = c("gray70", "tomato")) +
  theme_graph() +
  xlab("") +
  ylab("")

ggsave(filename = "../analysis/6_NetworkCons/PlateletNetwork/PPISubnetwork_FCs.tiff",
       device = "tiff",
       width = 20,
       height = 10)

# 2. DRUG REPURPOSING----

## load drugs data ----
drugbank <- read.csv("../analysis/7_DrugRepurposing/uniprot links.csv")
drug_targets_full <- drugbank %>% dplyr::select(Name, UniProt.ID)


## proteins of interest ----
#Sources <- PrioGenes$uniprot
#Targets <- TopGenes$uniprot

target_nodes = Targets

drug_targets <- drug_targets_full %>%
  filter(UniProt.ID %in% Sources) 

drugCounts <- drug_targets %>% 
  dplyr::count(Name, sort = TRUE)

top5Drugs <- drug_targets %>% 
  dplyr::count(Name, sort = TRUE) %>% 
  head(n=5) %>%
  pull(Name)

drug_targetsIn <- drug_targets

DTT <- drug_targetsIn %>%
  dplyr::left_join(annotdata, by = c("UniProt.ID" = "uniprot")) %>%
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
  tab_footnote(md("Only nodes with high gene scores (prio genes) are targeted."), 
               locations = cells_column_labels(columns = targets),
               placement = "right") %>%
  tab_options(data_row.padding = px(0.3)) %>%
  tab_options(table.font.size = 14) %>%
  cols_width(Name ~ px(350)) %>%
  gtsave("../analysis/7_DrugRepurposing/PrioDrugCountsPPI.pdf")

# tables of targets

drug_targetsIn <- drug_targets

TTT <- drug_targetsIn %>% 
  dplyr::left_join(annotdata, by = c("UniProt.ID" = "uniprot")) %>%
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
  tab_footnote(md("Only nodes with high gene scores (prio genes) are targeted."), 
               locations = cells_column_labels(columns = drugs),
               placement = "right") %>%
  tab_options(data_row.padding = px(0.3)) %>%
  tab_options(table.font.size = 12) %>%
  cols_width(drugs ~ px(500))%>%
  gtsave("../analysis/7_DrugRepurposing/PrioTargetCounts.docx")


#get drug induced networks

layoutOpts <- c(
  "layout_with_dh", "layout_with_drl", "layout_with_fr",
  "layout_with_gem", "layout_with_graphopt", "layout_with_kk",
  "layout_with_lgl", "layout_with_mds", "layout_with_sugiyama",
  "layout_as_bipartite", "layout_as_star", "layout_as_tree"
)



for (i in 1:dim(drugCounts)[1]){
  drugOI <- drugCounts[i, "Name"]
  
  source_nodes <- drug_targets %>%
    filter(Name == drugOI) %>%
    pull(UniProt.ID)
  ## get shortest paths
  collected_path_nodes = list()
  
  for(i_source in 1:length(source_nodes)){
    paths <- shortest_paths(subnetwork, from = source_nodes[[i_source]],
                            to = target_nodes,
                            weights = NA,
                            output = 'vpath')
    
    # Check if paths$vpath is NULL or empty
    if (!is.null(paths$vpath)) {
      path_nodes <- lapply(paths$vpath, names) %>% unlist() %>% unique()
      collected_path_nodes[[i_source]] <- path_nodes
    }
  }
 
   # Combine all collected path nodes and check if it's empty
  collected_path_nodes <- unlist(collected_path_nodes) %>% unique()
  
  if (length(collected_path_nodes) == 0) {
    next  # Skip the current iteration of the outer loop if no paths were collected
  }
  
  ## get drug network ----
  drug_nodes <- c(source_nodes,target_nodes, collected_path_nodes) %>%
    unique()
  drug_network <- induced_subgraph(graph = subnetwork,vids = drug_nodes)
  drug_network <- igraph::induced_subgraph(drug_network, 
                                           which(igraph::components(drug_network, mode = "weak")$membership == 
                                                   which.max(igraph::components(drug_network, mode = "weak")$csize)))
  
  
  V(drug_network)$node_type = ifelse(V(drug_network)$name %in% setdiff(Sources,source_nodes), "prio node",
                                     ifelse(V(drug_network)$name %in% source_nodes, "direct drug target",
                                            ifelse(
                                              V(drug_network)$name %in% target_nodes,"high-weight","intermediate node")))
  
  set.seed(123)
  ggraph(
    drug_network,
    layout = "dh"
  ) +
    geom_edge_fan(
      aes(
          start_cap = label_rect(node1.name),
          end_cap = label_rect(node2.name),
      ), color = "gray80",
      arrow = arrow(length = unit(0, 'mm')
      ),
      edge_width = .2,
      edge_alpha = .5
    ) +
    geom_node_point(aes(color = node_type), size = 8, alpha = 0.5) +
    geom_node_text(aes(label = symbol), color = "black",
                   position = "identity",
                   repel = F, size =3) +
    scale_color_manual(values = c("high-weight"="#F5E216","direct drug target"= "blue", "prio node" ="#A0CFEC" , 
                                  "intermediate node" = "gray30"),
                       guide = guide_legend(title = 'Node type')) +
    #scale_edge_color_manual(values = c("gray70","tomato"))+
    theme_graph() +
    xlab("") +
    ylab("") +
    ggtitle(paste0(drugOI, " induced network"))
  
  ggsave(filename = paste0("../analysis/7_DrugRepurposing/drug_networks_PPI/", drugOI, "DrugNetwork.tiff"),
         device = "tiff",
         width = 14,
         height = 8)
  
  
}  

