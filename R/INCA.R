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
    print('Invalid: `input` need to be either a data.frame or a list of RLE or data.frame.')
  }

  return(rle_object)
}


extractPeaksBySignals = function(peaks, threshold=0.5) {
  idr_col = grep('IDR', colnames(peaks), ignore.case=TRUE) ## column of IDR peaks
  if (length(idr_col) == 0) {
    print('Column of IDR peaks not found. Return the original file.')
    return(peaks)
  }
  if (length(idr_col) > 1) {
    print('Multiple columns of IDR peaks. Use the first column.')
    idr_col = idr_col[1]
  }

  signal_col = grep('signal(?!.*IDR)', colnames(peaks), ignore.case=TRUE, perl=TRUE) ## columns of non-IDR peaks

  if (threshold < 0 | threshold > 1) {
    print('Invalid quantile for `threshold`. Use default = 0.5.')
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
    print('Invalid `method` option. Use default = 1 (subtraction).')
    method = 1
  }
  if (method == 2 & pseudocount <= 0) {
    print('Invalid `pseudocount` for fold change. Use default = 0.1.')
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
    print('Invalid quantile for `threshold`. Use default = 0.5.')
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
    print('Invalid `window` size. Use default = 15.')
    window = 15
  }
  data = variants[, .(Chr,Start,End)]
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
          print('`cutoff` not provided. Use default cutoff = 30.')
          ls$cutoff = 30
        }
      }
      varInEnRegionByCounts(data, ls$counts, ls$cutoff, parallel)
    } else if (cases[2]) { ## varInEnRegionByPeaks()
      if (cases[4]) { #### threshold provided, extractPeaksBySignals()
        ls$peaks = extractPeaksBySignals(ls$peaks, ls$threshold)
      }
      varInEnRegionByPeaks(data, ls$peaks)
    } else {
      print('No files for `peaks` or `counts` provided.')
      return(NULL)
    }
  }

  return(unique(indices))
}


computeAllelicEffect = function(variants, CellLine1, CellLine2, CellLineNames, window=15, parallel=FALSE) {
  enriched = matrix(0, ncol=2, nrow=nrow(variants))

  ## Cell line 1
  indices = varInEnrichedRegion(variants, CellLine1, window, parallel)
  enriched[indices,1] = 1

  ## Cell line 2
  indices = varInEnrichedRegion(variants, CellLine1, window, parallel)
  enriched[indices,2] = 1

  variants[, 'EnrichedRegion'] = as.character(factor(paste(enriched[,1], enriched[,2], sep=':'),
                                                     levels = c('0:0','0:1','1:0','1:1'),
                                                     labels = c('Neither',CellLineNames,'Both')))
  return(variants)
}
