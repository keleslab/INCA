alignCellGenoToVar = function(variants, CellLine1, CellLine2, CellNames) {
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


scoreClinVarQSW = function(variants, SW, target, pathogenic=TRUE) {
  variants = as.data.table(variants)
  SW = as.data.table(SW)

  SW = mapSWToVariants(variants,SW)
  SW = SW[,grep('HepG2|K562',colnames(SW)), with=FALSE]
  SW = SW[,grep(paste0(target,'[|]'), colnames(SW)), with=FALSE]
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
  q = apply(q, 1, max)

  if (!'INCAscore' %in% colnames(variants)) {
    variants[,'INCAscore'] = 0
  }
  variants[,'CVQSWscore'] = q
  variants[,'INCAscore'] = variants[,'INCAscore'] + q

  return(variants)
}


scoreAllelicEffect = function(variants, CellLine1, CellLine2, CellNames, window=15, parallel=FALSE) {
  variants = as.data.table(variants)

  if (!'INCAscore' %in% colnames(variants)) {
    variants[,'INCAscore'] = 0
  }
  if ('AltCellLine' %in% colnames(variants)) {
    ALT = variants[,AltCellLine] %in% CellNames
  } else {
    ALT = TRUE
  }

  data = copy(variants)[ALT,]
  scores = list(varInEnrichedRegion(data, CellLine1, window, parallel), ## Cell line 1
                varInEnrichedRegion(data, CellLine2, window, parallel)) ## Cell line 2
  scores = sapply(scores, function(x) apply(x,1,max))

  enriched = apply(scores, 2, function(y) 1*(y>0))
  enriched = as.character(factor(paste(enriched[,1], enriched[,2], sep=':'),
                                 levels = c('0:0','1:0','0:1','1:1'),
                                 labels = c('Neither',CellNames,'Both')))
  variants[ALT, 'EnrichedRegion'] = enriched ## add to data

  AEscore = rep(0,nrow(variants))
  AEscore[ALT][enriched %in% CellNames] = apply(scores, 1, max)[enriched %in% CellNames]
  variants[,'AEscore'] = AEscore
  variants[, 'INCAscore'] = variants[,'INCAscore'] + AEscore

  return(variants)
}


scoreVarImpactOnGE = function(variants, CellLine1, CellLine2, CellNames, window=15, parallel=FALSE) {
  variants = as.data.table(variants)

  ## One gene in each row
  variants[, Gene := ifelse(gsub('\\s+','',Gene)=='','-',Gene)]
  variants[Gene=='', 'Gene'] = '-'
  variants = variants[, .(Gene = unlist(strsplit(Gene, '[^[:alnum:]\\-\\_\\.]', perl=TRUE))),
                      by = eval(grep('^Gene$', colnames(variants), invert=TRUE, value=TRUE))]

  if (!'INCAscore' %in% colnames(variants)) {
    variants[,'INCAscore'] = 0
  }
  if ('AltCellLine' %in% colnames(variants)) {
    ALT = variants[,AltCellLine] %in% CellNames
  } else {
    ALT = TRUE
  }
  data = copy(variants)[ALT,]

  if (!'list' %in% class(CellLine1)) {
    CellLine1 = list(CellLine1)
  }
  if (!'list' %in% class(CellLine2)) {
    CellLine2 = list(CellLine2)
  }

  DEscore = rep(0, nrow(variants))
  ## (i) The gene displays significant changes in expression upon knockdown in one of the cell lines
  de = data.table(mapDEResultToGenes(data, CellLine1[[1]], CellNames[1]), ## Cell line 1
                  mapDEResultToGenes(data, CellLine2[[1]], CellNames[2])) ## Cell line 2
  de = de[,grep('significant',colnames(de)),with=FALSE]

  DEgene = as.character(factor(paste(de[[1]], de[[2]], sep=':'),
                               levels = c('NA:NA','0:0','1:0','0:1','1:1'),
                               labels = c(NA,'Neither',CellNames,'Both')))
  variants[ALT, 'DEgene'] = DEgene ## add to data

  ## (ii) Whether the cell line with the expression change has a peak covering the location of the variant
  if (length(CellLine1) > 1 & length(CellLine2) > 1) {
    scores = list(varInEnrichedRegion(data, CellLine1[-1], window, parallel), ## Cell line 1
                  varInEnrichedRegion(data, CellLine2[-1], window, parallel)) ## Cell line 2
    scores = sapply(scores, function(x) apply(x,1,max))

    enriched = apply(scores, 2, function(y) 1*(y>0))
    enriched = as.character(factor(paste(enriched[,1], enriched[,2], sep=':'),
                                   levels = c('0:0','1:0','0:1','1:1'),
                                   labels = c('Neither',CellNames,'Both')))
    variants[ALT, 'DEbinding'] = enriched ## add to data
  } else if ('EnrichedRegion' %in% colnames(variants)) {
    enriched = variants[ALT,EnrichedRegion]
  } else {
    message('No `peak` files provided for assessing binding at the variant location.
            Max score for RBP-SNV impact on gene expression is 1.')
    enriched = NULL
  }

  DEscore[ALT][DEgene %in% CellNames] = 1
  DEscore[ALT][DEgene %in% CellNames & DEgene==enriched] = 2
  variants[, 'DEscore'] = DEscore
  variants[, 'INCAscore'] = variants[,'INCAscore'] + DEscore

  return(variants)
}
