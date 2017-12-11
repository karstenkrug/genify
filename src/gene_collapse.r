
gene.collapse <- function(exprs, anno, meth='median', id.col){
  
  #cat('----- :',id.col, '\n', mode(id.col), '\n')  
  
  ## expression data
  exprs.str <- apply(exprs, 1, paste, collapse='\t')
  id <- anno[, id.col]
  #save(exprs.str, exprs, file = 'tmp.RData')
  
  exprs.gc <- tapply(exprs.str, id, function(x){
    strsplit(x,'\t') %>% 
      unlist() %>% 
      as.numeric() %>% 
      matrix(nrow=length(x), byrow=T) %>% 
      apply(MARGIN = 2, FUN = match.fun(meth), na.rm=T)
    })
  exprs.gc <- matrix(unlist(exprs.gc), nrow=length(exprs.gc), byrow = T, dimnames=list(names(exprs.gc), colnames(exprs))) %>% data.matrix()
  
  ## annotation data
  anno.gc <- anno[match(row.names(exprs.gc), anno[, id.col]),]
  ##View(anno.gc)
  anno.gc <- data.frame( geneSymbol=anno.gc[, id.col],anno.gc[, -which(colnames(anno.gc) == id.col)], stringsAsFactors = F)
  
  ## remove not quantified rows
  rm.idx <- which(apply(exprs.gc, 1, function(x) sum(is.na(x))/length(x) ) == 1)
  if(length(rm.idx) > 0 ){
    exprs.gc <- exprs.gc[-rm.idx, ]
    anno.gc <- anno.gc[-rm.idx, ]
  }
  ## combine and return
  return(data.frame(anno.gc, exprs.gc, stringsAsFactors = F))
}
