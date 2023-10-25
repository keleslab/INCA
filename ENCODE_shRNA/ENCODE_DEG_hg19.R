library(data.table)
library(doMC)
cores=detectCores(); registerDoMC(cores=cores-5)

directory = '/z/Proj/pg_collaboration/RBP-analysis/ENCODE_shRNA'

file = fread(file.path(directory,'files.txt'), header=FALSE,nrow=1)[[1]]
metadata = fread(file)
metadata = metadata[`File assembly`=='hg19',] # hg19 genome assembly only
metadata = metadata[`File Status`=='released',] # Omit archived files
metadata = metadata[`Output type`=='differential expression quantifications',] # DEG analysis
metadata[,`Experiment target` := gsub('-human','', `Experiment target`)] # Drop '-human' suffix

target = unique(metadata[,.(`Experiment target`,`Biosample term name`)]) # Experiment targets

## Differential expression of genes from shRNA-seq
foreach(i = 1:nrow(target), .combine=c) %dopar% {
    files = suppressMessages( dplyr::left_join(target[i,], metadata, multiple='all') )

    deg = sapply(files[,`File download URL`], function(x) 'gene' %in% fread(x, nrow=1, header=FALSE))
    data = fread(files[deg, `File download URL`])                 
    
    output = paste(target[[1]][i], target[[2]][i], 'DGE.txt.gz', sep='_')
    fwrite(data, file.path(directory,'DEG',output), sep='\t')
}
