cepa.all.parallel <- function (mat = NULL, label = NULL, pc, 
          cen = default.centralities, cen.name = sapply(cen, function(x) ifelse(mode(x) == 
                                                                                  "name", deparse(x), x)), nlevel = "tvalue_abs", plevel = "mean", 
          iter = 1000, ncores = 2) 
{
  if (length(pc$pathList) < ncores) {
    stop("Number of cores should not be larger than the number of pathways.")
  }
  cat("Use snow package to do parallel computing (", ncores, 
      "cores ) ...\n")
  cepa.all.parallel.by.snow(mat = mat, 
                            label = label, pc = pc, cen = cen, cen.name = cen.name, 
                            nlevel = nlevel, plevel = plevel, iter = iter, ncores = ncores)
}

cepa.all.parallel.by.snow = function(mat = NULL, label = NULL, pc, cen = default.centralities,
                                                cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), 
                                                nlevel = "tvalue_abs", plevel = "mean", iter = 1000, ncores = 2) {
  # to avoid lazy evaluation
  mode(mat)
  mode(label)
  mode(pc)
  mode(cen)
  mode(cen.name)
  mode(nlevel)
  mode(plevel)
  mode(iter)
  
  d = divide(1:length(pc$pathList), ncores)
  if(dim(d)[1] < ncores) {
    ncores = dim(d)[1]
    cat("Since your task can only be divided into", ncores, "parts, modify the number of cores to", ncores, "\n")
  }
  
  cl = makeCluster(ncores, type="SOCK")
#   ignore = clusterCall(cl, function() {source("allNeededFunctions.r")})
ignore = clusterCall(cl, function() {sd <- "/home/user/allInfoDir/workspace/scripts/newCepaFunctions/"; file.sources <- list.files(sd, pattern = ".r", full.names = T);
sapply(file.sources, source, .GlobalEnv)
})
#   clusterExport(cl, c("pc", "d")) ## for testing
  ignore = clusterCall(cl, function() {library("igraph")})
  
  res = clusterApply(cl, 1:ncores, function(i) {
    pc <- set.pathway.catalogue(pathList = pc$pathList[d[i, 1]:d[i, 2]], interactionList = pc$interactionList, mapping = pc$mapping)
 
#     write(length(pc$pathList), file = paste0("testParallel.txt"), append = T, sep = "\t")
    cepa.all(mat = mat, label = label, pc = pc, cen = cen, nlevel = nlevel, plevel = plevel, iter = iter) #, cen.name = cen.name
  })
  stopCluster(cl)
  res.p = combine.cepa.all(res)
  return(res.p) 
}

divide = function(x, k) {
  if(length(x) ==1 && is.numeric(x)) {
    x = 1:x
  }
  if(length(x) < k) {
    k = length(x)
  }
  w = floor(length(x)/k)
  q = length(x) - k*w
  d = matrix(0, nrow=k, ncol=2)
  n = 1
  for(i in 1:k) {
    d[i, 1] = n
    d[i, 2] = n+w-1+ifelse(q>0, 1, 0)
    n = d[i,2]+1
    q = ifelse(q > 0, q-1, 0)
  }
  d[k,2] = length(x)
  return(d)
}

combine.cepa.all = function(res) {        
  obj = list()
  for(i in 1:length(res)) {
    obj = c(obj, res[[i]])
  }
  class(obj) = "cepa.all"
  return(obj)
}
