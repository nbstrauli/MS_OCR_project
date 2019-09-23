#!/usr/bin/ Rscript
library(stringr)
library(igraph)
library(rjson)
library(dplyr)

ptm <- proc.time()
source("init_shared.R")


el <-
  read.table(
    overall_cluster_edge_list_file,
    sep = "\t",
    header = T,
    stringsAsFactors = FALSE
  )
data <-
  read.table(
    data_for_clustering,
    sep = "\t",
    header = T,
    stringsAsFactors = FALSE
  )

subsets <- unique(data$subset)

# do quality control here
#rearrange the column, so first column is node_id
vertex <-
  data[, c(
    "Node_ID",
    "cluster_criteria",
    "compartment",
    "patient",
    "subset",
    "vgene",
    "jgene",
    "CDR3",
    "PATIENT",
    "count"
  )]

graph_all <- graph.data.frame(el, directed = F, vertices = vertex)

# Calculate the maximal (weakly or strongly) connected components of a graph
vertex$cluster_membership <- as.factor(clusters(graph_all)$membership)

graph_all <- graph.data.frame(el, directed = F, vertices = vertex)


# find cluster that appears in at least 2 subset
# i.e. membership count > 1
write_clusters_with_multiple_timepoint <-
  function (timepoint_count, vertex, graph_all) {
    # select count(subset), membership from vertex groupby membership
    # select membership where count_subset>= count
    membership_with_more_than_one_subset <- character()
    for (data_cluster in split(vertex, vertex$cluster_membership)) {
      timepoints <-sapply(str_split(data_cluster$subset, '_'), '[',2)
      number_of_timepoints <- length(unique(timepoints))
      if (number_of_timepoints >= timepoint_count) {
        print(unique(timepoints))
        membership_with_more_than_one_subset <-
          c(
            membership_with_more_than_one_subset,
            first(data_cluster$cluster_membership)
          )
      }
    }
    
    clusters_with_multiple_timepoint <- induced.subgraph(
      graph_all,
      which(
        V(graph_all)$cluster_membership %in% membership_with_more_than_one_subset
      )
    )
    
    gml_file_name <- file.path(
      result_dir,
      paste0(
        patient,
        "_clusters_at_least_",
        timepoint_count,
        "_timepoint_",
        config_suffix,
        ".gml"
      )
    )
    
    write.graph(clusters_with_multiple_timepoint, gml_file_name, format = "gml")
  }


write_clusters_with_multiple_timepoint(2, vertex, graph_all)
write_clusters_with_multiple_timepoint(3, vertex, graph_all)
write_clusters_with_multiple_timepoint(4, vertex, graph_all)
write_clusters_with_multiple_timepoint(5, vertex, graph_all)

paste(commandArgs(), collapse = " ")
print(proc.time() - ptm)
