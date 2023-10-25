#install.packages('/ua/jduan27/GVI/INCAData', repos=NULL, type='source')
library(INCAData)
library(doMC)
cores=detectCores(); registerDoMC(cores=cores)

directory = '/z/Proj/pg_collaboration/RBP-analysis/ENCODE_eCLIP'

dir.create(file.path(directory, 'NormRC'))
dir.create(file.path(directory, 'PeakSignals'))

file = fread(file.path(directory,'files.txt'), header=FALSE,nrow=1)[[1]]
metadata = fread(file)
metadata = metadata[`File assembly`=='hg19',] # hg19 genome assembly only
metadata = metadata[`File Status`=='released',] # Omit archived files
metadata[,`Experiment target` := gsub('-human','', `Experiment target`)] # Drop '-human' suffix

target = unique(metadata[,.(`Experiment target`,`Biosample term name`)]) # Experiment targets

## Peak signal values from multiple replicates
time = foreach(i = 1:nrow(target), .combine=c) %dopar% {
    files = suppressMessages( dplyr::left_join(target[i,], metadata[`File type`=='bed',], multiple='all') )

    experiment = files[,`File download URL`]
    suffix = sapply(files[,`Biological replicate(s)`], function(x) ifelse(x=='1, 2','IDR',paste0('Rep',x)))
    output = paste(target[[1]][i], target[[2]][i], 'PeakSignals.txt.gz', sep='_')                                   

    ## Computation time
    system.time({
        dt = summarizePeakSignals(experiment, suffix, output=file.path(directory, 'PeakSignals', output))
    })[3]

}
fwrite(cbind(target, ComputationTime = time), file.path(directory,'ComputationTime_summarizePeakSignals'), sep='\t')
                    
                    
## Normalized aligned read counts from eCLIP-seq 
time1 = foreach(i = 1:nrow(target), .combine=c) %dopar% {
    files = suppressMessages( dplyr::left_join(target[i,], metadata[`File type`=='bam',], multiple='all') )
    
    foreach(j = 1:nrow(files), .combine=c) %dopar% {
        experiment = files[j,`File download URL`]
        output = paste(target[[1]][i], target[[2]][i], 'NormRC', paste0('Rep',files[j,`Biological replicate(s)`],'.txt.gz'), sep='_')
        
        ## Computation time
        system.time({
            dt = summarizeReadCounts(experiment, output=file.path(directory, 'NormRC', output))
        })[3]
    } 
}
time1 = cbind(target[rep(1:nrow(target),each=2),], ComputationTime = time1)         
                    
                    
## Normalized aligned read counts from control eCLIP-seq 
time2 = foreach(i = 1:nrow(target), .combine=c) %dopar% {
    files = suppressMessages( dplyr::left_join(target[i,], metadata[`File type`=='bed',], multiple='all') )
    accession = gsub('.*/(E.*)/', '\\1', Reduce(intersect, strsplit(files[,`Derived from`], split=', '))) # eCLIP control 

    experiment = paste0('https://www.encodeproject.org/files/',accession,'/@@download/',accession,'.bam')
    output = paste(target[[1]][i], target[[2]][i], 'NormRC_Ctrl.txt.gz', sep='_')
    
    ## Computation time
    system.time({ 
        dt = summarizeReadCounts(experiment, output=file.path(directory, 'NormRC', output))
    })[3]
}
time2 = cbind(target, ComputationTime = time2)

fwrite(rbind(time1,time2), file.path(directory,'ComputationTime_summarizeReadCounts'), sep='\t')                    
