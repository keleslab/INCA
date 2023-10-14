# Convert data.table or data.frame to RLE class.
convertToRLE = function(input) {
  if (is.data.frame(input)) {
    stopifnot('Invalid: `input` must contain three columns in the following order: Chr, Length, and Value.'
              = ncol(input)==3)
    
    input = as.data.table(input)
    colnames(input) = c('Chr','Length','Value')
    
    rle_object = sapply(unique(input[,Chr]), function(chr) with(input[Chr==chr,], S4Vectors::Rle(Value, Length)))
  } else if (is.list(input)) {
    ## List elements that are not RLE objects 
    for (i in which(!sapply(input, is, 'Rle'))) {
      input[i] = convertToRLE(input[[i]])
    }
    rle_object = input
  } else {
    print('Invalid: `input` need to be either a data.table or a list of RLE or data.table.')
  }
  
  return(rle_object)
}

standardizeChr = function(chrom) {
  chrom = paste0('chr', gsub('^chr|chrom|chromosome', '', chrom, ignore.case=TRUE))
  return(chrom)
}  

#peak = fread('/z/Proj/pg_collaboration/RBP-analysis/ENCODE_eCLIP/PeakSignals/HNRNPK_K562_PeakSignals.txt.gz')

peaksBySignals = function(peak, threshold=0.5) {
  idr_col = grep('IDR', colnames(peak), ignore.case=TRUE) # column of IDR peaks
  stopifnot('No IDR peaks found in the peak file.' = length(idr_col)!=0)
  
  if (length(idr_col) > 1) {
    print('Multiple columns of IDR peaks found in the peak file. Use the first one.')
    idr_col = idr_col[1]
  }
  
  signal_col = grep('signal(?!.*IDR)', colnames(peak), ignore.case=TRUE, perl=TRUE) ## columns of non-IDR peaks
  
  if (threshold < 0 | threshold > 1) {
    print('Invalid quantile for threshold. Use default = 0.5.')
    threshold = 0.5
  }
  
  new_peak= lapply(signal_col, function(x) { ## new non-IDR peaks based on threshold
    cutoff = quantile(peak[[x]][peak[[idr_col]]>0], threshold)
    peak[peak[[x]]>cutoff, ] 
  }) 
  new_peak = c(new_peak, list(peak[peak[[idr_col]]>0,])) ## add IDR peaks
  new_peak = unique(rbindlist(new_peak))
  
  return(new_peak)
}

#input1 = fread('/z/Proj/pg_collaboration/RBP-analysis/ENCODE_eCLIP/NormRC/HNRNPK_K562_NormRC_Rep1.txt.gz') 
#input2 = fread('/z/Proj/pg_collaboration/RBP-analysis/ENCODE_eCLIP/NormRC/HNRNPK_K562_NormRC_Ctrl.txt.gz') 

#readCountsWRTControl(input1, input2, 1)

readCountsWRTControl = function(experiment, control, method=1, pseudocount=0.1) {
  rle_exp = convertToRLE(experiment)
  rle_ctrl = convertToRLE(control)
  
  if (!method %in% c(1,2)) {
    print('Invalid method option. Use default = 1 (subtraction).')
    method = 1
  }
  
  if (method == 2 & pseudocount <= 0) {
    print('Invalid pseudocount for fold change. Use default = 0.1.')
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

#readCountsCutoff(rc, peak[signalValue_IDR>0,], 0.5)

readCountsCutoff = function(counts, peak, threshold=0.5, parallel=FALSE) {
  counts = convertToRLE(counts)
  
  if (parallel) {
    peak_RC = foreach(i=1:nrow(peak), .combine=c) %dopar% {
      x = unlist(peak[i,])
      counts[[x[1]]][x[2]:x[3]]
    }
  } else {
    peak_RC = apply(peak, 1, function(x) counts[[x[1]]][x[2]:x[3]])
    peak_RC = do.call(c, peak_RC)
  }
  
  if (threshold < 0 | threshold > 1) {
    print('Invalid quantile for threshold. Use default = 0.5.')
    threshold = 0.5
  }   
  
  return(quantile(as.vector(peak_RC), threshold))                                                           
}


variantsInEnrichedRegion = function(variants, input, window=15, cutoff=30, parallel=FALSE) {
  variants[, Chr := standardizeChr(Chr)]
  
  if (window < 0) {
    print('Invalid window size. Use default = 15.')
    window = 15
  }
  data = variants[, .(Chr,Start,End)]
  data[, ':='(Start = Start - window, End = End + window, Index = 1:nrow(data))]
  
  if (is.data.frame(input)) {
    if (all(c('Chr','Start','End') %in% colnames(input))) { ## select by peak
      indices = c(input[data, on=.(Chr=Chr, Start<=Start, End>=Start), nomatch=NULL][,Index],
                  input[data, on=.(Chr=Chr, Start<=End, End>=End), nomatch=NULL][,Index],
                  input[data, on=.(Chr=Chr, Start>=Start, End<=End), nomatch=NULL][,Index])
      return(unique(indices))
    } else {
      input = convertToRLE(input)
    }
  } else if (is.list(input)) {
    input = convertToRLE(input)
  } else {
    print('Invalid input')
  }
  
  if (parallel) { ## select by aligned read counts
    indices = foreach(i = 1:nrow(data), .combine=c) %dopar% {
      x = unlist(data[i,])
      max(input[[x[1]]][x[2]:x[3]]@values) > cutoff
    }
  } else {
    indices = apply(data, 1, function(x) max(input[[x[1]]][x[2]:x[3]]@values) > cutoff)
  }
  indices = data[,Index][indices]   
  
  return(indices)
}


computeAllelicEffect = function(variants, CellLine1, CellLine2, CellLineNames, window=15, parallel=FALSE) {
  # Functions to find variants in enriched region given different input files
  cases_fn = function(ls) {
    if ('counts' %in% names(ls)) { ## determine whether in enriched region by read counts
      if (!'cutoff' %in% names(ls)) { ## cutoff not provided for read counts
        if ('peak' %in% names(ls)) { ## peak provided, then find cutoff by readCountsCutoff()
          if ('threshold' %in% names(ls)) {
            ls$cutoff = readCountsCutoff(ls$counts, ls$peak, ls$threshold)
          } else {
            print('Threshold not provided for readCountsCutoff(). Use default threshold quantile = 0.5 to compute cutoff.')
            ls$cutoff = readCountsCutoff(ls$counts, ls$peak)
          }
        } else {
          print('Cutoff not provided. Peak file not provided to compute cutoff by readCountsCutoff(). Use default cutoff = 30.')
          ls$cutoff = 30
        }
      }
      indices = variantsInEnrichedRegion(variants, ls$counts, window=window, cutoff=ls$cutoff, parallel=parallel)
    } else if ('peak' %in% names(ls)) { ## determine whether in enriched region by peak signals
      if ('threshold' %in% names(ls)) { ## threshold provide, then subset peaks by peaksBySignals()
        ls$peak = peaksBySignals(ls$peak, ls$threshold)
      }
      indices = variantsInEnrichedRegion(variants, ls$peak, window=window, parallel=parallel)
    } else {
      print('No peak or counts files provided .')
    }
    return(indices)
  }
  
  
  enriched = matrix(0, ncol=2, nrow=nrow(variants))
  
  ## Cell line 1
  if (parallel) {
    indices = foreach(x = CellLine1, .combine=c) %dopar% {
      cases_fn(x)
    }
  } else {
    indices = do.call(c, lapply(CellLine1, function(x) cases_fn(x)))
  }
  enriched[indices,1] = 1                             
  
  ## Cell line 2
  if (parallel) {
    indices = foreach(x = CellLine2, .combine=c) %dopar% {
      cases_fn(x)
    }
  } else {
    indices = do.call(c, lapply(CellLine2, function(x) cases_fn(x)))
  }
  enriched[indices,2] = 1 
  
  variants[, 'EnrichedRegion'] = factor(apply(enriched,1,paste,collapse=':'), levels=c('0:0','0:1','1:0','1:1'), 
                                        labels=c('Neither',CellLineNames,'Both'))
  return(variants)
}

