cepa.all <- function (mat = NULL, label = NULL, pc, cen = default.centralities, cen.name = sapply(cen, function(x) ifelse(mode(x) == "name", deparse(x), x)), nlevel = "tvalue_abs", plevel = "mean", iter = 1000) 
{
  res = cepa.univariate.all.mod(mat = mat, label = label, pc = pc, 
                              cen = cen, cen.name = cen.name, nlevel = nlevel, 
                              plevel = plevel, iter = iter)
  return(res)
}