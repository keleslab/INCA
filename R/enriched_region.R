## Variants that overlap with the specified peaks
varInEnRegionByPeaks = function(vdata, peaks) {
  vdata[,'Index'] = 1:nrow(vdata)

  peaks = as.data.table(peaks)
  colnames(peaks)[1:3] = c('Chr','Start','End')
  indices = c(peaks[vdata, on=.(Chr=Chr, Start<=Start, End>=Start), nomatch=NULL][,Index],
              peaks[vdata, on=.(Chr=Chr, Start<=End, End>=End), nomatch=NULL][,Index],
              peaks[vdata, on=.(Chr=Chr, Start>=Start, End<=End), nomatch=NULL][,Index])

  v = rep(0,nrow(vdata)); v[indices] = 1
  return(v)
}

## Variants that have aligned read counts below the lower cutoff and over the upper cutoff
varInEnRegionByRC = function(vdata, counts, cutoff, parallel=FALSE) {
  counts = convertToRle(counts)

  if (parallel) {
    v = foreach(i = 1:nrow(vdata), .combine=c) %dopar% {
      x = unlist(vdata[i,])
      max(counts[[x[1]]][x[2]:x[3]]@values)
    }
  } else {
    v = apply(vdata, 1, function(x) max(counts[[x[1]]][x[2]:x[3]]@values))
  }
  v = cut(v, breaks=c(-Inf,sort(cutoff),Inf), labels=c(('-1')[length(cutoff)-1],'0','1'))
  return(as.numeric(as.character(v)))
}

## Variants that reside in an enriched region.
varInEnrichedRegion = function(vdata, input, window=15, parallel=FALSE) {
  if (window < 0) {
    message('Invalid half-window size. Use default `window` = 15.')
    window = 15
  }
  data = as.data.table(vdata)[, .(Chr,Start,End)]
  data[, ':='(Chr = standardizeChr(Chr), Start = Start - window, End = End + window)]

  if (!all(sapply(input, function(x) 'list' %in% class(x)))) {
    input = list(input)
  }

  scores = foreach(ls = input) %dopar% {
    cases = c('counts','peaks','cutoff','threshold') %in% names(ls)
    scale = 1
    if (cases[1]) { ## varInEnRegionByRC()
      if (all(!cases[2:3])) { ### default cutoff
        message('cutoff for varInEnRegionByRC() not provided. Use default `cutoff` = c(5,15).')
        ls$cutoff = c(5,15); scale = 0.5
      } else if (!cases[3]) { ### cutoff not provided, peaks provided, computeRCQuantile()
        if (!cases[4]) { ### threshold not provided
          message('threshold for computeRCQuantile() not provided. Use default `threshold` = c(0.1,0.8).')
          ls$threshold = c(0.1,0.8)
        }
        ls$cutoff = computeRCQuantile(ls$counts, ls$peaks, ls$threshold, parallel)
        scale = max(ls$threshold)
      }
      varInEnRegionByRC(data, ls$counts, ls$cutoff, parallel) * scale
    } else if (cases[2]) { ## varInEnRegionByPeaks()
      if (cases[4]) { ### threshold provided, selectPeaksBySignals()
        ls$peaks = selectPeaksBySignals(ls$peaks, ls$threshold)
        scale = ls$threshold
      }
      varInEnRegionByPeaks(data, ls$peaks) * scale
    } else {
      message('No files for `peaks` or `counts` provided.')
      return(NA)
    }
  }

  return(do.call(cbind, scores))
}
