cepa.univariate.all.mod <- function (mat, label, pc, cen = c("equal.weight", "in.degree", "out.degree", "betweenness", "in.reach", "out.reach"), cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), nlevel = "tvalue_abs", plevel = "mean", iter = 1000) 
{
  nlevelFun = find.nlevelFun(nlevel)
  plevelFun = find.plevelFun(plevel)
  for (ce in cen) {
    if (is.function(ce)) {
      stop("Functions cannot be used directly because we need the function name, use quote or substitute.\n")
    }
  }
  if (class(pc) != "pathway.catalogue") {
    stop("pc argument should be a pathway.catalogue object.")
  }
  cat("  Calculate gene level values.\n")
  allScaledData = t(apply(mat, 1, scale))
  if (any(is.na(allScaledData))) {
    allScaledData <- allScaledData[-which(is.na(rowSums(allScaledData))),]
  }
#   write("Inside New Function", append = T, file = "testParallel.txt", sep = "\t")
  n.pathway = length(pc$pathList)
  pathway.name = names(pc$pathList)
  pathway.result = list()
  length(pathway.result) = n.pathway
  pathway.result = lapply(pathway.result, function(x) {
    y = list()
    length(y) = length(cen.name)
    names(y) = cen.name
    return(y)
  })
  names(pathway.result) = pathway.name
  cat("  Calculate pathway score...\n")
  # i = grep(id, names(pc$pathList))
  for (i in 1:length(pc$pathList)) {
    cat("    ", i, "/", length(pc$pathList), ", ", pathway.name[i], "...\n", sep = "")
    path = pc$pathList[[i]]
    inter = pc$interactionList[pc$interactionList[, 1] %in% path, 2:3]
    pathway = generate.pathway(as.matrix(inter))
#     cat("      Calculate node level value and permutate sample labels...\n")
    cat("      Calculate node level value and permutate gene values...\n")
    node = pathway.nodes(pathway)
    mapping = pc$mapping[pc$mapping[, 1] %in% node, ]
    l = rownames(mat) %in% mapping[, 2]
    cat("      ", sum(l), " genes measured in the pathway...\n", sep = "")
    if (sum(l) == 0) {
      node.level.from.expr = rep(0, length(node))
      node.level.t.value = rep(0, length(node))
      r.node.level.from.expr = matrix(0, nrow = length(node), ncol = iter)
    } else {
      mat.gene = mat[l, , drop = FALSE]
      mat.gene = t(apply(mat.gene, 1, scale))
      if (any(is.na(mat.gene))) {
        mat.gene <- mat.gene[-which(is.na(rowSums(mat.gene))),]
      }
      mat.node = matrix(0, nrow = length(node), ncol = dim(mat)[2])
      rownames(mat.node) = node
      for (k in 1:length(node)) {
        l = mapping[, 1] == node[k]
        gene.in.node = unique(mapping[l, 2])
        gene.in.node.in.mat = gene.in.node[gene.in.node %in% rownames(mat.gene)]
        if (length(gene.in.node.in.mat) == 1) {
          mat.node[k, ] = mat.gene[gene.in.node.in.mat, ]
        } else if (length(gene.in.node.in.mat) > 1) {
          mm = t(mat.gene[gene.in.node.in.mat, ])
          pcar = prcomp(mm)
          mat.node[k, ] = predict(pcar, mm)[, 1]
        }
      }
      node.level.from.expr = apply(mat.node, 1, function(x) nlevelFun(x[.treatment(label)], x[.control(label)]))
      node.level.t.value = apply(mat.node, 1, function(x) nodeLevelFun.tvalue(x[.treatment(label)], x[.control(label)]))
      node.level.from.expr[is.na(node.level.from.expr)] = 0
      node.level.t.value[is.na(node.level.t.value)] = 0
      r.node.level.from.expr = matrix(0, nrow = length(node), ncol = iter)
      for (k in 1:iter) {
        r.mat.node = mat.node
        tmp <- allScaledData[sample(1:dim(allScaledData)[1], sum(rowSums(mat.node) != 0)), ]
        r.mat.node[rowSums(mat.node) != 0, ] = tmp
        r.node.level.from.expr[, k] = apply(r.mat.node, 1, function(x) nlevelFun(x[.treatment(label)], x[.control(label)]))
        r.node.level.from.expr[is.na(r.node.level.from.expr[, k]), k] = 0
      }
    }
    j = 0
    for (ce in cen) {
      j = j + 1
      pathway.result[[i]][[j]] = cepa.univariate.mod(mat = mat, allScaledData = allScaledData,
                                                 label = label, pc = pc, pathway = pathway, cen = ce, 
                                                 iter = iter, nlevel = nlevel, plevel = plevel, 
                                                 node.level.from.expr = node.level.from.expr, 
                                                 node.level.t.value = node.level.t.value, r.node.level.from.expr = r.node.level.from.expr)
      cat("      - ", ce, ": ", round(pathway.result[[i]][[j]]$p.value, 3), "\n", sep = "")
    }
  }
  class(pathway.result) = "cepa.all"
  return(pathway.result)
}