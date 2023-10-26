scoreClinVarQSW = function(variants, SW, target, pathogenic=TRUE) {
  variants = as.data.table(variants)
  SW = as.data.table(SW)

  SW = mapSWToVariants(variants,SW)
  SW = SW[,grep(paste0(target,'[|]'), colnames(SW)),with=FALSE]
  colnames(SW) = gsub('(.+?)[|](.+?)[|].*rep(.*).hg19','\\2_\\1_SW_Rep\\3', colnames(SW)) # reformat colnames
  variants = cbind(variants, SW)

  ## Load quantiles of SeqWeaver from ClinVar
  directory = 'https://raw.github.com/jduan607/INCA/master/ClinVar_SeqWeaver_scores'
  if (pathogenic) {
    clinvar = fread(file.path(directory,'clinvar_pathogenic_snv_sw_quantile.txt.gz'))
  } else {
    clinvar = fread(file.path(directory,'clinvar_snv_sw_quantile.txt.gz'))
  }
  clinvar = clinvar[order(Quantile),] # smallest to largest

  ## ClinVar-quantiled SeqWeaver scores
  q = sapply(colnames(SW), function(x) {
    dt = SW[,x,with=FALSE]; setnames(dt, x, 'V')
    qt = clinvar[, c(x,'Quantile'), with=FALSE]; setnames(qt, x, 'V')
    qt = qt[dt, on=.(V<V), mult='last'][, Quantile := ifelse(is.na(Quantile), 0, Quantile)]
    qt[,Quantile]
  })

  if (!'INCAscore' %in% colnames(variants)) {
    variants[,'INCAscore'] = 0
  }
  variants[, INCAscore := INCAscore + apply(q, 1, max)]

  return(variants)
}


scoreAllelicEffect = function(variants, CellLine1, CellLine2, CellNames, window=15, parallel=FALSE) {
  variants = as.data.table(variants)
  enriched = matrix(0, ncol=2, nrow=nrow(variants))
  ## Cell line 1
  enriched[varInEnrichedRegion(variants, CellLine1, window, parallel), 1] = 1
  ## Cell line 2
  enriched[varInEnrichedRegion(variants, CellLine2, window, parallel), 2] = 1

  variants[, 'EnrichedRegion'] = as.character(factor(paste(enriched[,1], enriched[,2], sep=':'),
                                                     levels = c('0:0','0:1','1:0','1:1'),
                                                     labels = c('Neither',CellNames,'Both')))

  if ('AltCellLine' %in% colnames(variants)) {
    ALT = variants[,AltCellLine] %in% CellNames
  } else {
    ALT = TRUE
  }
  if (!'INCAscore' %in% colnames(variants)) {
    variants[,'INCAscore'] = 0
  }
  variants[ALT & EnrichedRegion %in% CellNames, INCAscore := INCAscore + 1]

  return(variants)
}


scoreVarImpactOnGE = function(variants, CellLine1, CellLine2, CellNames, window=15, parallel=FALSE) {
  variants = as.data.table(variants)

  ## One gene in each row
  variants[, Gene := gsub('\\s+','',Gene)]
  variants[Gene=='', 'Gene'] = '-'
  variants = variants[, .(Gene = unlist(strsplit(Gene, '[^[:alnum:]\\-\\_\\.]', perl=TRUE))),
                      by = eval(grep('^Gene$', colnames(variants), invert=TRUE, value=TRUE))]

  if (!'list' %in% class(CellLine1)) {
    CellLine1 = list(CellLine1)
  }
  if (!'list' %in% class(CellLine2)) {
    CellLine1 = list(CellLine2)
  }

  ## (i) The gene displays significant changes in expression upon knockdown in one of the cell lines
  de = data.table(mapDEResultToGenes(variants, CellLine1[[1]], CellNames[1]), ## Cell line 1
                  mapDEResultToGenes(variants, CellLine2[[1]], CellNames[2])) ## Cell line 2

  variants[, 'DESignificant'] =
    as.character(factor(apply(de[,grep('significant',colnames(de)),with=FALSE], 1, paste, collapse=':'),
                        levels = c('NA:NA','0:0','0:1','1:0','1:1'),
                        labels = c(NA,'Neither',CellNames,'Both')))

  if ('AltCellLine' %in% colnames(variants)) {
    ALT = variants[,AltCellLine] %in% CellNames
  } else {
    ALT = TRUE
  }
  if (!'INCAscore' %in% colnames(variants)) {
    variants[,'INCAscore'] = 0
  }
  variants[ALT & DESignificant %in% CellNames, INCAscore := INCAscore + 1]

  ## (ii) Whether the cell line with the expression change has a peak covering the location of the variant
  if (length(CellLine1) > 1 & length(CellLine2) > 1) {
    enriched = matrix(0, ncol=2, nrow=nrow(variants))
    ## Cell line 1
    enriched[varInEnrichedRegion(variants, CellLine1[-1], window, parallel), 1] = 1
    ## Cell line 2
    enriched[varInEnrichedRegion(variants, CellLine2[-1], window, parallel), 2] = 1

    variants[, 'DEBinding'] = as.character(factor(paste(enriched[,1], enriched[,2], sep=':'),
                                                  levels = c('0:0','0:1','1:0','1:1'),
                                                  labels = c('Neither',CellNames,'Both')))

    variants[ALT & DESignificant %in% CellNames & DESignificant == DEBinding, INCAscore := INCAscore + 1]
  } else if ('EnrichedRegion' %in% colnames(variants)) {
    variants[ALT & DESignificant %in% CellNames & DESignificant == EnrichedRegion, INCAscore := INCAscore + 1]
  } else {
    message('No `peak` or `counts` files provided for assessing binding at the variant location.
        Max score for RBP-SNV impact on gene expression is 1.')
  }
  return(variants)
}
