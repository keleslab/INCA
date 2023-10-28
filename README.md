# INCA

## Whole-genome Sequencing (WGS) Data

We extracted  the required information (Chr, Pos, Ref, Alt, Zyg) from the VCF files of WGS in HepG2 and K562 cell lines from the ENCODE project:
  + HepG2: ENCFF713BPG
  + K562: ENCFF752OAX
    
The corresponding code can be found in 'ExtractVCFInfo.py'.

Note: Zyg = 0 - Homozygous reference; 1 - Heterozygous; 2 - Unknown; 3 - Homozygous alternative

## Main Example in the Paper

Variants in the example are in high linkage-disequilibrium with the SNP _rs1057868_. INCAscore is computed using the data for the RBP _HNRNPK_.

### Load required data

```{r}
directory = 'https://raw.github.com/jduan607/INCA/master'

## GWAS
variants = fread(file.path(directory, 'GWAS', 'final_data.txt'))
SW = fread(file.path(directory,'GWAS','seqweaver_results.tsv'))
SW = SW[,grep('eCLIP',colnames(SW)),with=FALSE]

## DGE
dge1 = fread('https://www.encodeproject.org/files/ENCFF382MHL/@@download/ENCFF382MHL.tsv')[status=='OK',]
dge2 = fread('https://www.encodeproject.org/files/ENCFF527KUQ/@@download/ENCFF527KUQ.tsv')[status=='OK',]

## WGS
wgs1 = fread(file.path(directory,'ENCODE_WGS','HepG2_WGS.txt.gz'))
wgs2 = fread(file.path(directory,'ENCODE_WGS','K562_WGS.txt.gz'))

## ENCODE - HepG2
peak1 = fread(file.path(directory,'ENCODE_eCLIP/PeakSignals','HNRNPK_HepG2_PeakSignals.txt.gz'))

# Load if needed
#exp1.1 = fread(file.path(directory,'ENCODE_eCLIP/NormRC','HNRNPK_HepG2_NormRC_Rep1.txt.gz')) 
#exp1.2 = fread(file.path(directory,'ENCODE_eCLIP/NormRC','HNRNPK_HepG2_NormRC_Rep2.txt.gz')) 
#ctrl1 = fread(file.path(directory,'ENCODE_eCLIP/NormRC','HNRNPK_HepG2_NormRC_Ctrl.txt.gz'))

#rc1.1 = compareRCToControl(exp1.1, ctrl1)
#rc1.2 = compareRCToControl(exp1.2, ctrl1)

## ENCODE - K562
peak2 = fread(file.path(directory,'ENCODE_eCLIP/PeakSignals','HNRNPK_K562_PeakSignals.txt.gz'))

# Load if needed
#exp2.1 = fread(file.path(directory,'ENCODE_eCLIP/NormRC','HNRNPK_K562_NormRC_Rep1.txt.gz')) 
#exp2.2 = fread(file.path(directory,'ENCODE_eCLIP/NormRC','HNRNPK_K562_NormRC_Rep2.txt.gz')) 
#ctrl2 = fread(file.path(directory,'ENCODE_eCLIP/NormRC','HNRNPK_K562_NormRC_Ctrl.txt.gz'))

#rc2.1 = compareRCToControl(exp2.1, ctrl2)
#rc2.2 = compareRCToControl(exp2.2, ctrl2)
```

### Parallel backend

Parallel backend is not required to register, but it will reduce computing time. If not registered, set `parallel`=FALSE (as the default).

```{r}
library(doMC)
registerDoMC()
```

### (A) ClinVar-quantiled SeqWeaver scores

```{r}
variants = scoreClinVarQSW(variants, SW, 'HNRNPK', pathogenic=TRUE)
```

### (B) eCLIP-seq allelic effects

```{r}
epg1 = list(list(peaks=peak1[signalValue_IDR>0,]),
                list(peaks=peak1, threshold=0.5)) 

epg2 = list(list(peaks=peak2[signalValue_IDR>0,]),
            list(peaks=peak2, threshold=0.5))

# parallel = FALSE if no parallel backend registered
variants = scoreAllelicEffect(variants, epg1, epg2, c('HepG2','K562'), parallel=TRUE)
```

### (C) RBP-SNV impact on gene expression

```{r}
de1 = list(dge1[,.(gene,q_value)], # first element must be DE
           list(peaks=peak1, threshold=0.5)) 
de2 = list(dge2[,.(gene,q_value)], 
           list(peaks=peak2, threshold=0.5))
# parallel = FALSE if no parallel backend registered
variants = scoreVarImpactOnGE(variants, de1, de2, c('HepG2','K562'), parallel=TRUE)
```
