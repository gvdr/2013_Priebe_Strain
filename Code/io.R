#' input/output functions

#' `adjacency_from_paj()` return a weighted adjacency matrix from
#' a Pajek project file
#' it gets `name` of the file (with the explicit extension)
#' and the `base` path (ending with a /) as inputs
adjacency_from_paj <- function(Pajek_project){
  library(network)
  library(igraph)
  library(intergraph)
  
#' We invoke explicitly all the functions environemnt,
#' so to avoid any conflict between the three graph libraries
  Net_paj <- network::read.paj(Pajek_project)
  Net_net <- Net_paj$networks[[1]]
  Web <- intergraph::asIgraph(Net_net)
  igraph::V(Web)$name <- igraph::V(Web)$vertex.names
  weight_name <- igraph::list.edge.attributes(Web)[
    igraph::list.edge.attributes(Web) != "na"
  ]
  weights <- as.vector(igraph::get.edge.attribute(Web,weight_name))
  Web <- igraph::set_edge_attr(Web,"weight", value = weights)
  if("na" %in% igraph::list.edge.attributes(Web)){
    Web <- igraph::remove.edge.attribute(Web,"na")
  }
  Web <- igraph::remove.edge.attribute(Web,weight_name)
  if("na" %in% igraph::list.vertex.attributes(Web)){
    Web <- igraph::remove.vertex.attribute(Web,"na")
  }
  Web <- igraph::remove.vertex.attribute(Web,"vertex.names")
  
  Adjacency_test <- igraph::get.adjacency(Web,attr = "weight",sparse = F)

  return(Adjacency_test)
}