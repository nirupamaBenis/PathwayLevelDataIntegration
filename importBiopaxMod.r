"
@author Nirupama Benis
@version 0.2 Modified from import_biopax.r
@since 02-11-2015
<p>
This script can convert a biopax object (from rBiopaxParser) to a pathway.catalogue object
<p>
"

library(rBiopaxParser)
library(hash)
library(graph)
library(igraph)

ImportBiopaxMod = function (data, removePathways = NULL) {
  
  # == param
  # -biopax           The biopax data obtained from readBiopax from rBiopaxParser.
  # -rempvePathways   Any pathway names that you would like to remove (eg. Reactome categories)
  #
  # == details
  # The output from this functoin can be used in the cepa.all function to calculate pathway significance. 
  # The rBiopaxParser function pathway2RegulatoryGraph can be used instead of the modified function Pathway2GraphAllTypes. 
  # The difference is that, while building the pathway graph, the latter uses all the subclasses of biopax level 3 class Interaction, 
  # while the former starts from only the Control interactions.
  #
  # == return
  # A ``pathway.catalogue`` class object
  
  if(!inherits(data, "biopax")) {
    message("importing biopax")
    suppressMessages(biopax <- readBiopax(data, verbose = FALSE))
  } else {
    biopax = data
  }
  removePathways <- removePathways
#   if(biopax$biopaxlevel != 2) {
#     stop("Only Biopax level 2 is supported.")
#   }
  
  pathway_df = listPathways(biopax)
  pathway_df <- pathway_df[!(pathway_df[,2] %in% removePathways),]
  pathway_list = sapply(pathway_df[[1]], function(pid) {
   #changed it to accomodate changes in recent biopax level 3 files
  suppressWarnings(graph <- Pathway2GraphAll(biopax, pid, expandSubpathways = TRUE, splitComplexMolecules = FALSE, useIDasNodenames = TRUE, verbose = FALSE, withSubGraphs = T))
    
    if(!is.null(graph)) {
      if(length(graph::edges(graph)) == 0) {
        graph = NULL
      } else {
        edge = graph::edges(graph)
        input = rep(names(edge), sapply(edge, length))
        output = unlist(edge)
        interaction_id = paste(pid, seq_along(output), sep = "_")
        graph = data.frame(interaction.id = interaction_id, input = input, output = output, stringsAsFactors = FALSE)
      }
    }
    return(graph)
  })
  names(pathway_list) <- gsub(" ", ".", pathway_df[,2]) # added to make it look like the PID.db data
  pathway_list = pathway_list[!sapply(pathway_list, is.null)]
  pathList = lapply(pathway_list, function(pathway) pathway[[1]])
  interactionList = do.call("rbind", pathway_list)

  # nodes in pathway_list are complex ids
  all_nodes = c(interactionList[[2]], interactionList[[3]])
  all_nodes = unique(all_nodes)
  
  # mapping from nodes to molecules
  bp2 = selectInstances(biopax, id = all_nodes)
  l = isOfClass(bp2, "Complex")
  complex_nodes = unique(bp2[l]$id)
  non_complex_nodes = all_nodes[! all_nodes %in% complex_nodes]
  
  nl = c(lapply(complex_nodes, function(nid) {
    splitComplex(biopax, nid, returnIDonly = TRUE)
  }),
  lapply(non_complex_nodes, function(nid) {
    nid
  }))
  names(nl) = c(complex_nodes, non_complex_nodes)
  node_list = data.frame(node.id = rep(names(nl), sapply(nl, length)), molecule.id = unlist(nl), stringsAsFactors = FALSE)
  #changed NAME to displayName
  mapping = data.frame(node.id = node_list$node.id,
                       name = sapply(node_list$molecule.id, function(nid) getInstanceProperty(biopax, nid, property = "displayName", includeAllNames = F)),
                       class = sapply(node_list$molecule.id, function(nid) getInstanceClass(biopax, nid)),
                       stringsAsFactors = FALSE)
  mapping = mapping[mapping$class == "Protein", 1:2] #change this to only the ones you are able to map using ncbi2r
  
  node.type = sapply(all_nodes, function(nid) getInstanceClass(biopax, nid))
  node.name = sapply(all_nodes, function(nid) getInstanceProperty(biopax, nid, property = "displayName"))
  
  res = list(pathList = pathList,
             interactionList = interactionList,
             mapping = mapping,
             node.type = node.type,
             node.name = node.name)
  class(res) = "pathway.catalogue"
  return(res)
}