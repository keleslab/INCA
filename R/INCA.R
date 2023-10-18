standardizeChr = function(chrom) {
  chrom = paste0('chr', gsub('^chr|chrom|chromosome', '', chrom, ignore.case=TRUE))
  return(chrom)
}


# Convert data.table or data.frame to RLE class.
convertToRLE = function(input) {
  if (is.data.frame(input)) {
    stopifnot('Invalid: `input` must contain three columns in the following order: Chr, Length, and Value.'
              = ncol(input)==3)

    input = as.data.table(input)
    colnames(input) = c('Chr','Length','Value')

    rle_object = sapply(unique(input[,Chr]), function(chr) with(input[Chr==chr,], S4Vectors::Rle(Value, Length)))
  } else if ('list' %in% class(input)) {
    ## List elements that are not RLE objects
    for (i in which(!sapply(input, is, 'Rle'))) {
      input[i] = convertToRLE(input[[i]])
    }
    rle_object = input
  } else {
    message('Invalid: `input` need to be either a data.frame or a list of RLE or data.frame.')
  }

  return(rle_object)
}


extractPeaksBySignals = function(peaks, threshold=0.5) {
  idr_col = grep('IDR', colnames(peaks), ignore.case=TRUE) ## column of IDR peaks
  if (length(idr_col) == 0) {
    message('Column of IDR peaks not found. Return the original file.')
    return(peaks)
  }
  if (length(idr_col) > 1) {
    message('Multiple columns of IDR peaks. Use the first column.')
    idr_col = idr_col[1]
  }

  signal_col = grep('signal(?!.*IDR)', colnames(peaks), ignore.case=TRUE, perl=TRUE) ## columns of non-IDR peaks

  if (threshold < 0 | threshold > 1) {
    message('Invalid quantile for `threshold`. Use default = 0.5.')
    threshold = 0.5
  }

  new_peaks= lapply(signal_col, function(x) { ## new non-IDR peaks based on threshold
    cutoff = quantile(peaks[[x]][peaks[[idr_col]]>0], threshold)
    peaks[peaks[[x]]>cutoff, ]
  })
  new_peaks = c(new_peaks, list(peaks[peaks[[idr_col]]>0,])) ## add IDR peaks
  new_peaks = unique(rbindlist(new_peaks))

  return(new_peaks)
}


readCountsWRTControl = function(experiment, control, method=1, pseudocount=0.1) {
  rle_exp = convertToRLE(experiment)
  rle_ctrl = convertToRLE(control)

  if (!method %in% c(1,2)) {
    message('Invalid `method` option. Use default = 1 (subtraction).')
    method = 1
  }
  if (method == 2 & pseudocount <= 0) {
    message('Invalid `pseudocount` for fold change. Use default = 0.1.')
    pseducount = 0.1

  }

  if (method==1) {
    counts = sapply(names(rle_exp), function(x)
      suppressWarnings( rle_exp[[x]] - rle_ctrl[[x]] ))
  } else {
    counts = sapply(names(rle_exp), function(x)
      suppressWarnings( log2((rle_exp + pseudocount) / (rle_ctrl + pseudocount)) ))
  }
  return(counts)
}


findReadCountsCutoff = function(counts, peaks, threshold=0.5, parallel=FALSE) {
  counts = convertToRLE(counts)

  if (parallel) {
    peaks_RC = foreach(i=1:nrow(peaks), .combine=c) %dopar% {
      x = unlist(peaks[i,])
      counts[[x[1]]][x[2]:x[3]]
    }
  } else {
    peaks_RC = apply(peaks, 1, function(x) counts[[x[1]]][x[2]:x[3]])
    peaks_RC = do.call(c, peaks_RC)
  }

  if (threshold < 0 | threshold > 1) {
    message('Invalid quantile for `threshold`. Use default = 0.5.')
    threshold = 0.5
  }
  cutoff = quantile(as.vector(peaks_RC), threshold)

  return(cutoff)
}


varInEnRegionByPeaks = function(vdata, peaks) {
  vdata[,'Index'] = 1:nrow(vdata)

  colnames(peaks)[1:3] = c('Chr','Start','End')
  indices = c(peaks[vdata, on=.(Chr=Chr, Start<=Start, End>=Start), nomatch=NULL][,Index],
              peaks[vdata, on=.(Chr=Chr, Start<=End, End>=End), nomatch=NULL][,Index],
              peaks[vdata, on=.(Chr=Chr, Start>=Start, End<=End), nomatch=NULL][,Index])

  return(unique(indices))
}

varInEnRegionByCounts = function(vdata, counts, cutoff, parallel=FALSE) {
  counts = convertToRLE(counts)

  if (parallel) {
    indices = foreach(i = 1:nrow(vdata), .combine=c) %dopar% {
      x = unlist(vdata[i,])
      max(counts[[x[1]]][x[2]:x[3]]@values) > cutoff
    }
  } else {
    indices = apply(vdata, 1, function(x) max(counts[[x[1]]][x[2]:x[3]]@values) > cutoff)
  }
  return(which(indices))
}

varInEnrichedRegion = function(variants, input, window=15, parallel=FALSE) {
  if (window < 0) {
    message('Invalid `window` size. Use default = 15.')
    window = 15
  }
  data = as.data.table(variants)[, .(Chr,Start,End)]
  data[, ':='(Chr = standardizeChr(Chr), Start = Start - window, End = End + window)]

  if (!all(sapply(input, function(x) 'list' %in% class(x)))) {
    input = list(input)
  }

  indices = foreach(ls = input, .combine=c) %dopar% {
    cases = c('counts','peaks','cutoff','threshold') %in% names(ls)
    if (cases[1]) { ## varInEnRegionByCounts()
      if (!cases[3]) { ### cutoff not provided
        if (cases[2]) { #### peaks provided, findReadCountsCutoff()
          ls$cutoff = findReadCountsCutoff(ls$counts, ls$peaks, ifelse(cases[4],ls$threshold,-1))
        } else {
          message('`cutoff` not provided. Use default cutoff = 30.')
          ls$cutoff = 30
        }
      }
      varInEnRegionByCounts(data, ls$counts, ls$cutoff, parallel)
    } else if (cases[2]) { ## varInEnRegionByPeaks()
      if (cases[4]) { ### threshold provided, extractPeaksBySignals()
        ls$peaks = extractPeaksBySignals(ls$peaks, ls$threshold)
      }
      varInEnRegionByPeaks(data, ls$peaks)
    } else {
      message('No files for `peaks` or `counts` provided.')
      return(NULL)
    }
  }

  return(unique(indices))
}


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
