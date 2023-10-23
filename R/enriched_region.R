## Variants that overlap with the specified peaks
varInEnRegionByPeaks = function(vdata, peaks) {
  vdata[,'Index'] = 1:nrow(vdata)

  peaks = as.data.table(peaks)
  colnames(peaks)[1:3] = c('Chr','Start','End')
  indices = c(peaks[vdata, on=.(Chr=Chr, Start<=Start, End>=Start), nomatch=NULL][,Index],
              peaks[vdata, on=.(Chr=Chr, Start<=End, End>=End), nomatch=NULL][,Index],
              peaks[vdata, on=.(Chr=Chr, Start>=Start, End<=End), nomatch=NULL][,Index])

  return(unique(indices))
}

## Variants that have aligned read counts exceeding the cutoff
varInEnRegionByRC = function(vdata, counts, cutoff, parallel=FALSE) {
  counts = convertToRle(counts)

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
    if (cases[1]) { ## varInEnRegionByRC()
      if (!cases[3]) { ### cutoff not provided
        if (cases[2]) { #### peaks provided, computeRCQuantile()
          ls$cutoff = computeRCQuantile(ls$counts, ls$peaks, ifelse(cases[4],ls$threshold,-1))
        } else {
          message('`cutoff` not provided. Use default cutoff = 30.')
          ls$cutoff = 30
        }
      }
      varInEnRegionByRC(data, ls$counts, ls$cutoff, parallel)
    } else if (cases[2]) { ## varInEnRegionByPeaks()
      if (cases[4]) { ### threshold provided, selectPeaksBySignals()
        ls$peaks = selectPeaksBySignals(ls$peaks, ls$threshold)
      }
      varInEnRegionByPeaks(data, ls$peaks)
    } else {
      message('No files for `peaks` or `counts` provided.')
      return(NULL)
    }
  }

  return(unique(indices))
}
