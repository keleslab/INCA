mapSWToVariants = function(variants, SW) {
  variants = as.data.table(variants)
  SW = as.data.table(SW)[pos!='pos',]

  SW = SW[, pos := as.numeric(pos)]
  SW = SW[variants, on=.(chrom=Chr, pos=Pos, ref=Ref, alt=Alt), mult='first'][, 1:ncol(SW)]
  return(SW)
}

mapDEResultToGenes = function(variants, DEresult, CellName=NULL) {
  data = as.data.table(variants)[,'Gene']
  suppressMessages(
    suppressWarnings(
      data[, Gene := HGNChelper::checkGeneSymbols(Gene,unmapped.as.na=FALSE)[,'Suggested.Symbol']]))

  ## Separate rows to have one gene in each row
  DEresult = as.data.table(DEresult)
  DEresult = DEresult[, .(gene = unlist(strsplit(gene, '[^[:alnum:]\\-\\_\\.]', perl=TRUE))),
                      by=eval(colnames(DEresult)[-1])] # gene moves to the last column
  suppressMessages(
    suppressWarnings(
      DEresult[, gene := HGNChelper::checkGeneSymbols(gene,unmapped.as.na=FALSE)[,'Suggested.Symbol']]))

  DEresult = DEresult[, .SD[which.min(1)], by = gene] ## 1st column is p-value/q-value
  DEresult[, 'significant'] = (DEresult[[2]] <= 0.05) * 1 ## 2nd column is p-value/q-value
  data = DEresult[data, on=.(gene=Gene)][,-1] ## match to data

  if (!is.null(CellName)) {
    colnames(data) = paste(CellName, colnames(data), sep='_')
  }
  return(data)
}
