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

## ClinVar-quantiled SeqWeaver scores
scoreClinVarQSW = function(variants, SW, target, empirical = 0) {
    variants = as.data.table(variants)
    SW = as.data.table(SW)
    SW = mapSWToVariants(variants,SW)

    if (!'INCAscore' %in% colnames(variants)) {
        variants[,'INCAscore'] = 0
    }

    if (is.numeric(empirical)) {
        ## Load quantiles of SeqWeaver computed using ClinVar SNVs
        directory = 'https://raw.github.com/jduan607/INCA/master/ClinVar_SeqWeaver_scores'
        if (empirical == 0) {
            quantiles = fread(file.path(directory,'clinvar_snv_sw_quantile.txt.gz'))
        } else if (empirical == 1) {
            quantiles = fread(file.path(directory,'clinvar_pathogenic_snv_sw_quantile.txt.gz'))
        } else {
            message('Invalid `empirical` option. Use empirical quantiles of SeqWeaver scores from all ClinVar SNVs.')
            quantiles = fread(file.path(directory,'clinvar_snv_sw_quantile.txt.gz'))
        }
        ## Reformat column names to match column names in quantiles
        colnames(SW) = gsub('(.+?)[|](.+?)[|].*rep(.*).hg19','\\2_\\1_SW_Rep\\3', colnames(SW))
    } else {
        ## User-provided quantiles of SeqWeaver
        quantiles = as.data.table(empirical)
    }

    SW = SW[, grep(paste0(target,'[^A-Za-z0-9]'), colnames(SW)), with=FALSE] # extract target
    quantiles = quantiles[order(Quantile),] # smallest to largest

    ## ClinVar-quantiled SeqWeaver scores
    QSWscore = lapply(colnames(SW), function(x) {
      if (!x %in% colnames(SW)) {
          message(paste(paste('Empirical/estimated quantiles of', x, 'are not provided.'),
                        paste(x, 'is excluded for computation of ClinVar-quantiled SeqWeaver scores.'), sep='\n'))
          rep(0, nrow(SW))
      } else {
          dt = SW[,x,with=FALSE]; setnames(dt, x, 'V')
          qt = quantiles[, c(x,'Quantile'), with=FALSE]; setnames(qt, x, 'V')
          qt = qt[dt, on=.(V<V), mult='last'][, Quantile := ifelse(is.na(Quantile), 0, Quantile)]
          qt[,Quantile]
      } })
    QSWscore = apply(do.call(cbind, QSWscore), 1, max)

    variants = cbind(variants, SW)
    variants[,'QSWscore'] = QSWscore
    variants[,'INCAscore'] = variants[,'INCAscore'] + QSWscore

    return(variants)
}


scoreAllelicEffect = function(variants, CellLine1, CellLine2, CellNames=paste('Cell type',1:2), window=15, parallel=FALSE) {
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
  scores = sapply(scores, function(x)
    apply(x, 1, function(i) ifelse(any(i<0), min(i), max(i))))

  ## IDR peak
  enriched = apply(scores, 2, function(j) 1*(j==1))
  enriched = as.character(factor(paste(enriched[,1], enriched[,2], sep=':'),
                                  levels = c('0:0','1:0','0:1','1:1'),
                                  labels = c('Neither',CellNames,'Both')))

  ## Read counts/non-IDR peak
  enriched2 = apply(scores, 2, function(j) ifelse(j<0, -1, ifelse(j>0,1,0)))
  enriched2 = factor(paste(enriched2[,1], enriched2[,2], sep=':'))
  if (any(scores<0)) {
    enriched[enriched2=='1:-1' & enriched=='Neither'] = CellNames[1]
    enriched[enriched2=='-1:1' & enriched=='Neither'] = CellNames[2]
  } else {
    enriched[enriched2=='1:0' & enriched=='Neither'] = CellNames[1]
    enriched[enriched2=='0:1' & enriched=='Neither'] = CellNames[2]
  }
  enriched[enriched2=='1:1' & enriched=='Neither'] = 'Both'

  AEscore = rep(0,nrow(variants))
  AEscore[ALT][enriched %in% CellNames] = apply(scores, 1, max)[enriched %in% CellNames]
  variants[ALT, 'EnrichedRegion'] = enriched
  variants[,'AEscore'] = AEscore
  variants[, 'INCAscore'] = variants[,'INCAscore'] + AEscore

  return(variants)
}

## Scores for the RBP-variant impact on the gene expression
scoreVarImpactOnGE = function(variants, CellLine1, CellLine2, CellNames=paste('Cell type',1:2)) {
  variants = as.data.table(variants)

  ## Separate rows to have one gene in each row
  variants[, Gene := ifelse(gsub('\\s+','',Gene)=='','-', gsub('\\s+','',Gene))]
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

  ## (ii) Whether the cell line with the expression change has a peak covering the location of the variant
  if (length(CellLine1) > 1 & length(CellLine2) > 1) {
    scores = list(varInEnrichedRegion(data, CellLine1[-1]), ## Cell line 1
                  varInEnrichedRegion(data, CellLine2[-1])) ## Cell line 2
    scores = sapply(scores, function(x) apply(x,1,max))

    ## IDR peak
    enriched = apply(scores, 2, function(j) 1*(j==1))
    enriched = as.character(factor(paste(enriched[,1], enriched[,2], sep=':'),
                                   levels = c('0:0','1:0','0:1','1:1'),
                                   labels = c('Neither',CellNames,'Both')))

    ## non-IDR peak
    enriched2 = apply(scores, 2, function(j) 1*(j>0))
    enriched2 = factor(paste(enriched2[,1], enriched2[,2], sep=':'))
    enriched[enriched2=='1:0' & enriched=='Neither'] = CellNames[1]
    enriched[enriched2=='0:1' & enriched=='Neither'] = CellNames[2]
    enriched[enriched2=='1:1' & enriched=='Neither'] = 'Both'
  } else if ('EnrichedRegion' %in% colnames(variants)) {
    enriched = variants[ALT,EnrichedRegion]
  } else {
    message('No `peak` files provided for assessing binding at the variant location.
            Max score for RBP-SNV impact on gene expression is 1.')
    enriched = NULL
  }

  DEscore[ALT][DEgene %in% CellNames] = 1
  DEscore[ALT][DEgene %in% CellNames & DEgene==enriched] = 2
  variants[ALT, 'DEgene'] = DEgene
  variants[ALT, 'DEbinding'] = enriched
  variants[, 'DEscore'] = DEscore
  variants[, 'INCAscore'] = variants[,'INCAscore'] + DEscore

  return(variants)
}
