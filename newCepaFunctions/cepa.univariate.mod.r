cepa.univariate.mod <- function (mat, allScaledData, label, pc, pathway = NULL, id = NULL, cen = "equal.weight", cen.name = if (is.function(cen)) deparse(substitute(cen)) else if (mode(cen) == "name") deparse(cen) else cen, iter = 1000, nlevel = "tvalue_abs", plevel = "mean", node.level.from.expr = NULL, node.level.t.value = NULL, r.node.level.from.expr = NULL) {
  
  nlevelFun = find.nlevelFun(nlevel)
  plevelFun = find.plevelFun(plevel)
#   node.level.from.expr = NULL
#   r.node.level.from.expr = NULL
  if (!is.null(pathway)) {
    if (is.matrix(pathway) || is.data.frame(pathway)) {
      if (length(dim(pathway)) != 2 || dim(pathway)[2] != 
            2) {
        stop("if pathway is a matrix or data frame, it should be 2 dimension and the number of columns is 2.\n")
      }
      pathway = generate.pathway(pathway)
    } else if (class(pathway) != "igraph") {
      stop("Since pathway is not formatted as edge list, it should be an igraph object.")
    }
  } else if (!is.null(id)) {
    path = pc$pathList[[id]]
    inter = pc$interactionList[pc$interactionList[, 1] %in% path, 2:3]
    pathway = generate.pathway(inter)
  } else {
    stop("You should specify pathway argument or id argument.")
  }
  if (iter < 100) {
    stop("Iterations should not be smaller than 100.\n")
  }
  if (length(cen) != 1) {
    stop("Length of cen must be equal to 1.\n")
  }
  weight = centrality(pathway, cen) 
  add = 0
  if (sum(weight == 0) != length(weight)) {
    add = min(weight[weight > 0])/100
  }
  weight = weight + ifelse(sum(weight == weight[1]) == length(weight), 0, add)
  node = pathway.nodes(pathway)
  mapping = pc$mapping[pc$mapping[, 1] %in% node, ]
  node.name = node
  member = character(0)
  for (i in 1:length(node)) {
    l = mapping[, 1] == node[i]
    if (sum(l)) {
      member = sort(unique(mapping[l, 2]))
      node.name[i] = paste(member, collapse = "\n")
    }
  }
  if (is.null(node.level.from.expr)) {
    l = rownames(mat) %in% mapping[, 2]
    if (sum(l) == 0) {
      node.level.from.expr = rep(0, length(node))
      node.level.t.value = rep(0, length(node))
    }     else {
      mat.gene = mat[l, , drop = FALSE]
      mat.gene = t(apply(mat.gene, 1, scale))
      mat.node = matrix(0, nrow = length(node), ncol = dim(mat)[2])
      rownames(mat.node) = node
      for (i in 1:length(node)) {
        l = mapping[, 1] == node[i]
        gene.in.node = unique(mapping[l, 2])
        gene.in.node.in.mat = gene.in.node[gene.in.node %in% rownames(mat.gene)]
        if (length(gene.in.node.in.mat) == 1) {
          mat.node[i, ] = mat.gene[gene.in.node.in.mat, ]
        }   else if (length(gene.in.node.in.mat) > 1) {
          mm = t(mat.gene[gene.in.node.in.mat, ])
          pcar = prcomp(mm)
          mat.node[i, ] = predict(pcar, mm)[, 1]
        }
      }
      if (any(is.na(mat.node))) {
        tmpPos <- which(is.na(rowSums(mat.node)))
        mat.node <- mat.node[-tmpPos, ]
        weight <- weight[-tmpPos]
      }
      node.level.from.expr = apply(mat.node, 1, function(x) nlevelFun(x[.treatment(label)], x[.control(label)]))
      node.level.t.value = apply(mat.node, 1, function(x) nodeLevelFun.tvalue(x[.treatment(label)], x[.control(label)]))
    }
  }
  node.level.from.expr[is.na(node.level.from.expr)] = 0
  node.level.t.value[is.na(node.level.t.value)] = 0
  node.level = weight * node.level.from.expr
  ds = quantile(node.level, c(1, 0.75, 0.5, 0))
  names(ds) = c("max", "q75", "median", "min")
  s = plevelFun(node.level)
  s.random = numeric(iter)
  ds.random = matrix(0, iter, 4)
  colnames(ds.random) = c("max", "q75", "median", "min")
  if (sum(rownames(mat) %in% mapping[, 2]) == 0) { ## change to more than 1
    s.random = rep(0, iter)
  } else {
    for (i in 1:iter) {
      if (is.null(r.node.level.from.expr)) {
        #         r.label = .permutate(label) ## no need because we do not have too many samples
        r.mat.node = mat.node
        tmp <- allScaledData[sample(1:dim(allScaledData)[1], sum(rowSums(mat.node) != 0)), ]
        r.mat.node[rowSums(mat.node) != 0, ] = tmp
        #         print(rownames(tmp))
        r.node.level.from.expr.current = apply(r.mat.node, 1, function(x) nlevelFun(x[.treatment(label)], x[.control(label)]))
      } else {
        r.node.level.from.expr.current = r.node.level.from.expr[, i]
      }
      r.node.level.from.expr.current[is.na(r.node.level.from.expr.current)] = 0
      r.node.level = weight * r.node.level.from.expr.current
      s.random[i] = plevelFun(r.node.level)
      ds.random[i, ] = quantile(r.node.level, c(1, 0.75, 0.5, 0))
    }
  }
  p.value = (sum(s.random >= s))/(iter)
  res = list(score = s, score.distribution = ds, score.random = s.random, 
             score.distribution.random = ds.random, p.value = p.value, 
             centrality = cen.name, weight = weight, node.level.t.value = node.level.t.value, 
             node.level = node.level, node.name = node.name, pathway = pathway, 
             framework = "gsa.univariate")
  class(res) = "cepa"
  return(invisible(res))
}