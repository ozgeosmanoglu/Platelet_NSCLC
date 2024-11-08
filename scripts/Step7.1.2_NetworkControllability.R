# Load packages ----
suppressPackageStartupMessages({ 
  library(dplyr)
  library(ggplot2)
  library(igraph)
  library(ggraph)
  library(gt)
})

#1.  NETWORK CONTROLLABILITY ---- 

## load the network controllability results ----
NCresults <- readxl::read_excel("../analysis/6_NetworkCons/PlateletNetwork/NetworkAnalysis.xlsx")
NCresults <- left_join(NCresults, annotdata, by = "symbol")

indis <- NCresults %>%
  filter(Classification == "indispensable")

indisnodes <- indis %>% pull(`shared name`)

critical <- NCresults %>%
  filter(Group == "critical") 

criticalnodes <- critical%>%
  pull(`shared name`)


EdgeAtt <- as_data_frame(net2, "edges")

# get shortest paths from critical to indispensable with edge weights
collected_path_nodes = list()

for(i in 1:length(criticalnodes)){
  
  paths <- shortest_paths(net2, from = criticalnodes[i],
                          to = indisnodes,
                          weights = E(net2)$weight,
                          mode = "out",
                          algorithm = "dijkstra", #positive weights
                          output = 'vpath')
  path_nodes <- lapply(paths$vpath,names) %>% unlist() %>% unique()
  collected_path_nodes[[i]] <- path_nodes
}

collected_path_nodes <- unlist(collected_path_nodes) %>% unique()
#collected_path_nodes <- path_nodes

# create network
subnodes <- c(criticalnodes, indisnodes, collected_path_nodes) %>%
  unique()
critIndisNetwork <- induced_subgraph(graph = net2,vids = subnodes)

## remove nonconnected elements - get largest connected component
components <- igraph::components(critIndisNetwork, mode="weak")
biggest_cluster_id <- which.max(components$csize)
vert_ids <- V(critIndisNetwork)[components$membership == biggest_cluster_id]
critIndisNetwork <- igraph::induced_subgraph(critIndisNetwork, vert_ids)

#visualize network
V(critIndisNetwork)$node_type = ifelse(V(critIndisNetwork)$name %in% indisnodes, "indispensable",
                                   ifelse(V(critIndisNetwork)$name %in% criticalnodes, "critical",
                                         "intermediate"))

criticalsInNet <- V(critIndisNetwork)[V(critIndisNetwork)$node_type == "critical"]$name
indisInNet <- V(critIndisNetwork)[V(critIndisNetwork)$node_type == "indispensable"]$name

write.table(as_data_frame(critIndisNetwork, "e"), "../analysis/6_NetworkCons/PlateletNetwork/Subnetwork.txt",
            row.names = F, sep = "\t")

write.table(as_data_frame(critIndisNetwork, "v"), "../analysis/6_NetworkCons/PlateletNetwork/SubnetworkNodes.txt",
            row.names = F, sep = "\t")

set.seed(297)
ggraph(
  critIndisNetwork,
  layout = "dh"
) +
  geom_edge_fan(
    aes(edge_colour = type,
        start_cap = label_rect(node1.name),
        end_cap = label_rect(node2.name),
    ),
    arrow = arrow(length = unit(1, 'mm')
    ),
    edge_width = .2,
    edge_alpha = .5
  ) +
  geom_node_point(aes(color = node_type), size = 8, alpha = 0.5) +
  geom_node_text(aes(label = symbol), color = "black",
                 position = "identity",
                 repel = F, size =3) +
  scale_color_manual(values = c("#A0CFEC","#F5E216","gray30"),
                     guide = guide_legend(title = 'Node type')) +
  scale_edge_color_manual(values = c("gray70","tomato"))+
  theme_graph() +
  xlab("") +
  ylab("") 

ggsave(filename = "../analysis/6_NetworkCons/PlateletNetwork/SubnetworkStyle1.tiff",
       device = "tiff",
       width = 20,
       height = 10)

set.seed(297)
ggraph(
  critIndisNetwork,
  layout = "dh"
) +
  geom_edge_fan(
    aes(edge_colour = type,
        start_cap = label_rect(node1.name),
        end_cap = label_rect(node2.name),
    ),
    arrow = arrow(length = unit(1, 'mm')),
    edge_width = .2,
    edge_alpha = .5
  ) +
  geom_node_point(aes(color = logFC), size = 8, alpha = 0.7) +  # Color by logFC
  geom_node_text(aes(label = symbol), color = "black",
                 position = "identity",
                 repel = F, size = 3) +
  scale_color_gradient2(low = "blue", high = "red", midpoint = 0, 
                        guide = guide_legend(title = 'logFC')) +  # Gradient scale for logFC
  scale_edge_color_manual(values = c("gray70", "tomato")) +
  theme_graph() +
  xlab("") +
  ylab("")

ggsave(filename = "../analysis/6_NetworkCons/PlateletNetwork/SubnetworkStyle1_FCs.tiff",
       device = "tiff",
       width = 20,
       height = 10)



set.seed(297)
ggraph(
  critIndisNetwork,
  layout = "dh"
) +
  geom_edge_fan(
    aes(edge_colour = type,
        start_cap = label_rect(node1.name),
        end_cap = label_rect(node2.name),
    ),
    arrow = arrow(length = unit(1, 'mm')
    ),
    edge_width = .3,
    edge_alpha = .5
  ) +
  geom_node_point(size = 0.1, color = "white", alpha = 0) +
  geom_node_label(aes(label = symbol, color = node_type),
                  position = "identity",
                  repel = F, size =4, label.size = 1) +
  geom_node_label(aes(label = symbol),
                  position = "identity",
                  repel = F, size =4, label.size = NA) +
  scale_color_manual(values = c("#A0CFEC","#F5E216","gray30"),
                     guide = guide_legend(title = 'Node type')) +
  scale_edge_color_manual(values = c("gray80","tomato"))+
  theme_graph() +
  xlab("") +
  ylab("")

ggsave(filename = "../analysis/6_NetworkCons/PlateletNetwork/SubnetworkStyle2.tiff",
       device = "tiff",
       width = 24,
       height = 12)

# 2. DRUG REPURPOSING----

## load drugs data ----
drugbank <- read.csv("../analysis/7_DrugRepurposing/uniprot links.csv")
drug_targets_full <- drugbank %>% dplyr::select(Name, UniProt.ID)

## proteins of interest ----
drug_targets <- drug_targets_full %>%
  filter(UniProt.ID %in% criticalsInNet) 

drugCounts <- drug_targets %>% 
  dplyr::count(Name, sort = TRUE)

top5Drugs <- drug_targets %>% 
  dplyr::count(Name, sort = TRUE) %>% 
  head(n=5) %>%
  pull(Name)


indisdrugTargets <- drug_targets_full %>%
  filter(UniProt.ID %in% indisInNet) 

#drug_targetsIn <- drug_targets
drug_targetsIn <- indisdrugTargets


DTT <- drug_targetsIn %>%
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
  tab_footnote(md("Only indispensable nodes are targeted."), 
               locations = cells_column_labels(columns = targets),
               placement = "right") %>%
  tab_options(data_row.padding = px(0.3)) %>%
  tab_options(table.font.size = 14) %>%
  cols_width(Name ~ px(350))%>%
  gtsave("../analysis/7_DrugRepurposing/IndisDrugCounts.pdf")

 
# tables of targets
TTT <- drug_targetsIn %>% 
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
  tab_footnote(md("Only indispensable nodes are targeted."), 
               locations = cells_column_labels(columns = drugs),
               placement = "right") %>%
  tab_options(data_row.padding = px(0.3)) %>%
  tab_options(table.font.size = 12) %>%
  cols_width(drugs ~ px(500))%>%
  gtsave("../analysis/7_DrugRepurposing/IndisTargetCounts.pdf")


# get subnetworks of drugs and visualize

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
  
  
  ## get shortest paths ----
  collected_path_nodes = list()

  for(i_source in 1:length(source_nodes)){
    paths <- shortest_paths(critIndisNetwork, from = source_nodes[[i_source]],
                            to = indisInNet,
                            weights = NA,
                            mode = "out",
                            algorithm = "dijkstra",
                            output = 'vpath')
    path_nodes <- lapply(paths$vpath,names) %>% unlist() %>% unique()
    collected_path_nodes[[i_source]] <- path_nodes
  }
  
  collected_path_nodes <- unlist(collected_path_nodes) %>% unique()

  ## get drug network ----
  drug_nodes <- c(source_nodes,indisInNet, collected_path_nodes) %>%
    unique()
  drug_network <- induced_subgraph(graph = critIndisNetwork,vids = drug_nodes)
  drug_network <- igraph::induced_subgraph(drug_network, 
                                           which(igraph::components(drug_network, mode = "weak")$membership == 
                                                   which.max(igraph::components(drug_network, mode = "weak")$csize)))
  
  
  V(drug_network)$node_type = ifelse(V(drug_network)$name %in% setdiff(criticalsInNet,source_nodes), "critical",
                                     ifelse(V(drug_network)$name %in% source_nodes, "direct drug target (critical)",
                                            ifelse(V(drug_network)$name %in% indisInNet,"indispensable","intermediate node")))

  
  set.seed(123)
  ggraph(
    drug_network,
    layout = "dh"
  ) +
    geom_edge_fan(
      aes(edge_colour = type,
          start_cap = label_rect(node1.name),
          end_cap = label_rect(node2.name),
      ),
      arrow = arrow(length = unit(1, 'mm')
      ),
      edge_width = .2,
      edge_alpha = .5
    ) +
    geom_node_point(aes(color = node_type), size = 8, alpha = 0.5) +
    geom_node_text(aes(label = symbol), color = "black",
                    position = "identity",
                    repel = F, size =3) +
    scale_color_manual(values = c("critical"="#A0CFEC","direct drug target (critical)"= "blue", "indispensable" = "#F5E216", 
                                  "intermediate node" = "gray30"),
                       guide = guide_legend(title = 'Node type')) +
    scale_edge_color_manual(values = c("gray70","tomato"))+
    theme_graph() +
    xlab("") +
    ylab("") +
    ggtitle(paste0(drugOI, " induced network"))
  
  ggsave(filename = paste0("../analysis/7_DrugRepurposing/drug_networks/", drugOI, "DrugNetwork.tiff"),
         device = "tiff",
         width = 14,
         height = 8)
  
  
}  



