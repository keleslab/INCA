standardizeChr = function(chrom) {
  chrom = paste0('chr', gsub('^chr|chrom|chromosome', '', chrom, ignore.case=TRUE))
  return(chrom)
}


## Convert data frame(s) to a list of Rle objects.
convertToRle = function(input) {
  if (is.data.frame(input)) {
    stopifnot('Invalid: `input` must contain three columns in the order of chromosomes, lengths, and values.'
              = ncol(input)==3)

    input = as.data.table(input)
    colnames(input) = c('Chr','Length','Value')
    input[, Chr := standardizeChr(Chr)]

    rle_object = sapply(unique(input[,Chr]), function(chr) with(input[Chr==chr,], S4Vectors::Rle(Value, Length)))
  } else if ('list' %in% class(input)) {
    for (i in which(!sapply(input, is, 'Rle'))) {
      input[i] = convertToRle(input[[i]]) # list elements that are not Rle objects
    }
    names(input) = standardizeChr(names(input))
    rle_object = input
  } else {
    message('Invalid: `input` must be a data frame or a list of data frames.')
  }

  return(rle_object)
}


## Calculate aligned read counts with respect to control
compareRCToControl = function(experiment, control, method=1, pseudocount=0.1) {
  rle_exp = convertToRle(experiment)
  rle_ctrl = convertToRle(control)

  if (!method %in% c(1,2)) {
    message('Invalid method option. Use default `method` = 1 (subtraction).')
    method = 1
  }
  if (method == 2 & pseudocount <= 0) {
    message('Invalid pseudocount for fold change. Use default `pseudocount` = 0.1.')
    pseducount = 0.1
  }

  if (method==1) {
    counts = sapply(names(rle_exp), function(x)
      suppressWarnings( rle_exp[[x]] - rle_ctrl[[x]] ))
  } else {
    counts = sapply(names(rle_exp), function(x)
      suppressWarnings( (rle_exp + pseudocount) / (rle_ctrl + pseudocount) ))
  }
  return(counts)
}


## Select IDR peaks and non-IDR peaks with signal values exceeding the threshold
selectPeaksBySignals = function(peaks, threshold=0.5) {
  idr_col = grep('IDR', colnames(peaks), ignore.case=TRUE) # column of IDR peaks
  if (length(idr_col) == 0) {
    message('Column of IDR peaks not found. Return the original file.')
    return(peaks)
  }
  if (length(idr_col) > 1) {
    message('Multiple columns of IDR peaks. Use the first column.')
    idr_col = idr_col[1]
  }

  signal_col = grep('signal(?!.*IDR)', colnames(peaks), ignore.case=TRUE, perl=TRUE) # columns of non-IDR peaks

  if (threshold < 0 | threshold > 1) {
    message('Invalid threshold. Use default `threshold` = 0.5.')
    threshold = 0.5
  }

  peaks = as.data.table(peaks)
  new_peaks= lapply(signal_col, function(x) { ## new non-IDR peaks based on threshold
    cutoff = quantile(peaks[[x]][peaks[[idr_col]]>0], threshold)
    peaks[peaks[[x]]>cutoff, ]
  })
  new_peaks = c(new_peaks, list(peaks[peaks[[idr_col]]>0,])) ## add IDR peaks
  new_peaks = unique(rbindlist(new_peaks))

  return(new_peaks)
}


## Calculate the quantile of aligned read counts
computeRCQuantile = function(counts, peaks, threshold=0.5, parallel=FALSE) {
  counts = convertToRle(counts)

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
    message('Invalid threshold. Use default `threshold` = 0.5.')
    threshold = 0.5
  }
  cutoff = quantile(as.vector(peaks_RC), threshold)

  return(cutoff)
}
