centrality = function(graph, method="equal.weight") {
  
  if(length(method) > 1) {
    stop("Length of method must be equal to 1.\n")
  }
  
  if(is.function(method)) {
    return(method(graph))
  } else if(mode(method) == "name") {
    method = eval(method)
    return(method(graph))
  } else if(method == "equal.weight") {
    return(rep(1, vcount(graph)))
  } else if(method == "in.degree") {
    return(igraph::degree(graph, mode="in"))
  } else if(method == "out.degree") {
    return(igraph::degree(graph, mode="out"))
  } else if(method == "degree") {
    return(igraph::degree(graph))
  } else if(method == "betweenness") {
    return(igraph::betweenness(graph))
  } else if(method == "in.reach") {
    return(reach(graph, mode="in"))
  } else if(method == "out.reach") {
    return(reach(graph, mode="out"))
  } else if(method == "reach") {
    return(reach(graph))
  } else if(method == "in.spread") {
    return(spread(graph, mode="in"))
  } else if(method == "out.spread") {
    return(spread(graph, mode="out"))
  } else if(method == "spread") {
    return(spread(graph))
  } else {
    stop("Wrong centrality method.\n")
  }
}

reach = function(graph, weights=E(graph)$weight, mode=c("all", "in", "out")) {
  mode = mode[1]
  sp = shortest.paths(graph, weights=weights, mode=mode)
  s = apply(sp, 1, function(x) {
    if(all(x == Inf)) {
      return(0)
    }
    else {
      return(max(x[x != Inf]))
    }
  })
  return(s)
}

find.nlevelFun = function(method) {
  if(is.character(method)) {
    f = switch(method,
               tvalue     = nodeLevelFun.tvalue,
               tvalue_sq  = nodeLevelFun.tvalue_sq,
               tvalue_abs = nodeLevelFun.tvalue_abs,
               stop("Default node-level functions in string format are only tvalue, tvalue_sq and tvalue_abs."))
  } else if ( is.function(method) ) {
    if(length(as.list(args(method))) != 3) {
      stop("Self-defined node-level function only can have two arguments.")
    }
    f = method
  }
  return(f)
}

find.plevelFun = function(method) {
  if(is.character(method)) {
    f = switch(method,
               max    = pathwayLevelFun.max,
               min    = pathwayLevelFun.min,
               median = pathwayLevelFun.median,
               sum    = pathwayLevelFun.sum,
               mean   = pathwayLevelFun.mean,
               rank   = pathwayLevelFun.rank,
               stop("Wrong pathway level method."))
  } else if ( is.function(method) ) {
    if(length(as.list(args(method))) != 2) {
      stop("Self-defined pathway-level function only can have one argument.")
    }
    f = method
  }
  return(f)
}

node.score = function(gene.level = NULL, gene.in.node = NULL, nlevelFun = NULL) {
  
  node.level.from.expr = numeric(length(gene.in.node))
  for(i in 1:length(gene.in.node)) {
    node.level.from.expr[i] = ifelse(length(gene.in.node[[i]]), nlevelFun(gene.level[names(gene.level) %in% gene.in.node[[i]]]), 0)
  }
  return(node.level.from.expr)
}

nodeLevelFun.tvalue = function(x1, x2, ...) {
  n1 = length(x1)
  n2 = length(x2)
  v1 = var(x1)
  v2 = var(x2)
  ifelse(v1 + v2 == 0, 0, (mean(x1) - mean(x2)) / sqrt(v1/n1 + v2/n2))
}
nodeLevelFun.tvalue_sq = function(x1, x2, ...) {
  nodeLevelFun.tvalue(x1, x2)^2
}
nodeLevelFun.tvalue_abs = function(x1, x2, ...) {
  abs(nodeLevelFun.tvalue(x1, x2))
}


pathwayLevelFun.max = function(x) {
  max(x, na.rm = TRUE)
}
pathwayLevelFun.min = function(x) {
  min(x, na.rm = TRUE)
}
pathwayLevelFun.median = function(x) {
  median(x, na.rm = TRUE)
}
pathwayLevelFun.sum = function(x) {
  sum(x, na.rm = TRUE)
}
pathwayLevelFun.mean = function(x) {
  mean(x, na.rm = TRUE)
}
pathwayLevelFun.rank = function(x) {
  wilcox.test(x, exact = FALSE)$statistic
}

sampleLabel = function (label, treatment, control) {
  if (sum(label == treatment) == 0 | sum(label == control) == 0) {
    stop("Can not find treatment label or control label.")
  }
  res = list(label = label, treatment = treatment, control = control)
  class(res) = "sampleLabel"
  return(res)
}

.treatment = function(sl) {
  return(which(sl$label == sl$treatment))
}

.control = function(sl) {
  return(which(sl$label == sl$control))
}

.permutate = function(sl) {
  sl$label = sample(sl$label, length(sl$label), replace = FALSE)
  return(sl)
}

.factor = function(sl) {
  return(factor(sl$label))
}

generate.pathway <- function (el) 
{
  if (dim(el)[2] != 2) {
    stop("Second dimension of edgelist should be 2.\n")
  }
  el.vector = apply(el, 1, paste, collapse = "--")
  el.vector = unique(el.vector)
  el = t(as.matrix(as.data.frame(strsplit(el.vector, "--"))))
  g = graph.edgelist(el, directed = TRUE)
  return(g)
}

pathway.nodes <- function (pathway) 
{
  if (class(pathway) != "igraph") {
    stop("Wrong class for pathway.\n")
  }
  n = vcount(pathway)
  name = get.vertex.attribute(pathway, "name")
  if (length(name)) {
    return(name)
  }
  else {
    return(1:n)
  }
}

p.table <- function (x, adj.method = NA, cutoff = ifelse(adj.method == "none", 
                                                         0.01, 0.05)) 
{
  if (class(x) != "cepa.all") {
    stop("x should be cepa.all object.\n")
  }
  n.pathway = length(x)
  p.value = matrix(0, nrow = length(x), ncol = length(x[[1]]))
  for (i in 1:length(x)) {
    p.value[i, ] = sapply(x[[i]], function(x) x$p.value)
  }
  rownames(p.value) = names(x)
  colnames(p.value) = names(x[[1]])
  if (!is.na(adj.method)) {
    p.value = apply(p.value, 2, p.adjust, adj.method)
    l = apply(p.value, 1, function(x) {
      sum(x <= cutoff) > 0
    })
    p.value = p.value[l, , drop = FALSE]
  }
  return(p.value)
}

set.pathway.catalogue <- function (pathList, interactionList, mapping, min.node = 5, max.node = 5000, 
                                   min.gene = min.node, max.gene = max.node, ...) 
{
  if (!is.list(pathList)) {
    stop("pathList should be a list.\n")
  }
  if (is.null(names(pathList))) {
    stop("pathList should have names.\n")
  }
  if (length(dim(interactionList)) != 2) {
    stop("interactionList should be two dimension matrix.\n")
  }
  if (dim(interactionList)[2] != 3) {
    stop("interactinList should contain 3 columns.\n")
  }
  l = sapply(pathList, function(x) {
    it = interactionList[interactionList[, 1] %in% x, 2:3]
    node = unique(c(it[, 1], it[, 2]))
    l.node = length(node)
    gene = unique(mapping[mapping[, 1] %in% node, 2])
    l.gene = length(gene)
    return(c(l.node, l.gene))
  })
  pathList = pathList[l[1, ] >= min.node & l[1, ] <= max.node & 
                        l[2, ] >= min.gene & l[2, ] <= max.gene]
  if (length(pathList) == 0) {
    warning("Your pathway catalogue is empty!")
  }
  pc = list(pathList = pathList, interactionList = interactionList, 
            mapping = mapping, ...)
  class(pc) = "pathway.catalogue"
  return(pc)
}

sampleLabel <- function (label, treatment, control) 
{
  if (sum(label == treatment) == 0 | sum(label == control) == 
        0) {
    stop("Can not find treatment label or control label.")
  }
  res = list(label = label, treatment = treatment, control = control)
  class(res) = "sampleLabel"
  return(res)
}