mapCellGenoToVar = function(variants, CellLine1, CellLine2, CellNames) {
  variants = as.data.table(variants)

  if (!'Zyg' %in% colnames(variants)) {
    variants[,'Zyg'] = 1
    message('Zygosity not provided. Assume all variants are heterozygous.')
  }

  ## Cell line 1
  colnames(CellLine1) = c('Chr','Pos','Ref', paste('C1',c('Alt','Alt2','Zyg'),sep='_'))
  CellLine1 = CellLine1[variants[,.(Chr,Pos,Ref)], on=.(Chr=Chr, Pos=Pos, Ref=Ref)][,-c(1:3)]
  ## Cell line 2
  colnames(CellLine2) = c('Chr','Pos','Ref', paste('C2',c('Alt','Alt2','Zyg'),sep='_'))
  CellLine2 = CellLine2[variants[,.(Chr,Pos,Ref)], on=.(Chr=Chr, Pos=Pos, Ref=Ref)][,-c(1:3)]

  type = paste(apply(CellLine1, 1, function(x) any(!is.na(x))*1),
               apply(CellLine2, 1, function(x) any(!is.na(x))*1), sep=':')

  cdata = cbind(variants[,.(Alt,Zyg)], CellLine1, CellLine2, AltCellLine='Neither')
  cdata[type=='1:0' & Alt==C1_Alt & C1_Alt2=='', 'AltCellLine'] = CellNames[1]
  cdata[type=='0:1' & Alt==C2_Alt & C2_Alt2=='', 'AltCellLine'] = CellNames[2]

  cdata[type=='1:1' & Alt==C1_Alt & Alt==C2_Alt, 'AltCellLine'] = 'Both'
  cdata[AltCellLine=='Both' & Zyg==C1_Zyg & Zyg!=C2_Zyg, 'AltCellLine'] = CellNames[1]
  cdata[AltCellLine=='Both' & Zyg==C2_Zyg & Zyg!=C1_Zyg, 'AltCellLine'] = CellNames[2]

  cdata[C1_Alt2!='', C1_Alt := paste(C1_Alt,C1_Alt2,sep=',')]
  cdata[C2_Alt2!='', C2_Alt := paste(C2_Alt,C2_Alt2,sep=',')]

  cdata = cdata[,.(C1_Alt,C1_Zyg,C2_Alt,C2_Zyg,AltCellLine)]
  colnames(cdata) = gsub('C1',CellNames[1],colnames(cdata))
  colnames(cdata) = gsub('C2',CellNames[2],colnames(cdata))

  variants = cbind(variants, cdata)
  return(variants)
}


mapDEResultToGenes = function(variants, DEresult, CellName=NULL) {
  data = as.data.table(variants)[,'Gene']
  suppressMessages(
    suppressWarnings(
      data[, Gene := HGNChelper::checkGeneSymbols(Gene,unmapped.as.na=FALSE)[,'Suggested.Symbol']]))

  DEresult = as.data.table(DEresult)
  DEresult = DEresult[, .(gene = unlist(strsplit(gene, '[^[:alnum:]\\-\\_\\.]', perl=TRUE))),
                      by=eval(colnames(DEresult)[-1])]
  suppressMessages(
    suppressWarnings(
      DEresult[, gene := HGNChelper::checkGeneSymbols(gene,unmapped.as.na=FALSE)[,'Suggested.Symbol']]))

  DEresult = DEresult[, .SD[which.min(1)], by = gene] ## 1: column of p-value/q-value
  DEresult[, 'significant'] = (DEresult[[2]] <= 0.05) * 1 ## 2: column of p-value/q-value
  data = DEresult[data, on=.(gene=Gene)][,-1] ## match to data

  if (!is.null(CellName)) {
    colnames(data) = paste(CellName, colnames(data), sep='_')
  }
  return(data)
}


mapSWToVariants = function(variants, SW) {
  variants = as.data.table(variants)
  SW = as.data.table(SW)

  SW[, pos := as.numeric(pos)]
  SW = SW[variants, on=.(chrom=Chr, pos=Pos, ref=Ref, alt=Alt), mult='first'][, 1:ncol(SW)]
  return(SW)
}
