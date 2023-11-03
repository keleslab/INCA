# INCA

## Installation

```{r}
devtools::install_github('keleslab/INCA')
```

### Parallel backend

Parallel backend is not required to register, but it will reduce computing time. If not registered, set `parallel`=FALSE (as the default).

```{r}
library(doMC)
registerDoMC()
```

## Input

INCA takes as input a variant file with the following fields: _Chr_, _Pos_, _Start_, _End_, _Ref_, _Alt_, and _Gene_. It then derives three scores for the effect of each variant on each RBP activity using available cell lines (K562 and HepG2 for this implementation): (i) ClinVar-quantiled SeqWeaver scores; (ii) a score of allelic effect on RBP binding derived from a pre-computed library of RBP eCLIP-seq experiments (Van Nostrand et al., 2020); and (iii) a gene impact score for the gene that the SNV is proximal to based on a pre-computed library of accompanying RNA-seq experiments (wild type and RBP knockdown by shRNA).

Chr|Pos|Start|End|Ref|Alt|Gene
---|---|---|---|---|---|---
chr7|75607155|75607155|75607155|A|G|POR
chr7|75611756|75611756|75611756|C|T|POR
...
chr7|75636240|75636240|75636240|T|C|STYXL1


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
library(INCA)

directory = 'https://raw.github.com/jduan607/INCA/master'

## GWAS
variants = fread(file.path(directory, 'GWAS', 'final_data.txt'))
SW = fread(file.path(directory,'GWAS','seqweaver_results.tsv'))
SW = SW[,c(1:8, grep('HepG2|K562',colnames(SW))),with=FALSE] # the example focuses on HepG2 and K562

## DEG
deg1 = fread(file.path(directory,'ENCODE_shRNA/DEG','HNRNPK_HepG2_DEG.txt.gz'))
deg2 = fread(file.path(directory,'ENCODE_shRNA/DEG','HNRNPK_K562_DEG.txt.gz'))

## WGS
wgs1 = fread(file.path(directory,'ENCODE_WGS','HepG2_WGS.txt.gz'))
wgs2 = fread(file.path(directory,'ENCODE_WGS','K562_WGS.txt.gz'))

## ENCODE - peak
peak1 = fread(file.path(directory,'ENCODE_eCLIP/PeakSignals','HNRNPK_HepG2_PeakSignals.txt.gz'))
peak2 = fread(file.path(directory,'ENCODE_eCLIP/PeakSignals','HNRNPK_K562_PeakSignals.txt.gz'))

## ENCODE - read counts (Load if needed)
exp1.1 = fread(file.path(directory,'ENCODE_eCLIP/NormRC','HNRNPK_HepG2_NormRC_Rep1.txt.gz')) 
exp1.2 = fread(file.path(directory,'ENCODE_eCLIP/NormRC','HNRNPK_HepG2_NormRC_Rep2.txt.gz')) 
ctrl1 = fread(file.path(directory,'ENCODE_eCLIP/NormRC','HNRNPK_HepG2_NormRC_Ctrl.txt.gz'))

rc1.1 = compareRCToControl(exp1.1, ctrl1)
rc1.2 = compareRCToControl(exp1.2, ctrl1)

exp2.1 = fread(file.path(directory,'ENCODE_eCLIP/NormRC','HNRNPK_K562_NormRC_Rep1.txt.gz')) 
exp2.2 = fread(file.path(directory,'ENCODE_eCLIP/NormRC','HNRNPK_K562_NormRC_Rep2.txt.gz')) 
ctrl2 = fread(file.path(directory,'ENCODE_eCLIP/NormRC','HNRNPK_K562_NormRC_Ctrl.txt.gz'))

rc2.1 = compareRCToControl(exp2.1, ctrl2)
rc2.2 = compareRCToControl(exp2.2, ctrl2)
```

### Align genotypes of K562 and HepG2 cell lines to variants

```{r}
variants = alignCellGenoToVar(variants, wgs1, wgs2, c('HepG2','K562'))
```

### (A) ClinVar-quantiled SeqWeaver scores

```{r}
variants = scoreClinVarQSW(variants, SW, 'HNRNPK', empirical=1)
```

### (B) eCLIP-seq allelic effects

```{r}
epg1 = list(list(peaks=peak1[signalValue_IDR>0,]),
            list(peaks=peak1, threshold=0.5), # optional
            list(counts=rc1.1, peaks=peak1[signalValue_IDR>0,], threshold=c(0.1,0.8)), # optional
            list(counts=rc1.2, peaks=peak1[signalValue_IDR>0,], threshold=c(0.1,0.8))) # optional
            
epg2 = list(list(peaks=peak2[signalValue_IDR>0,]),
            list(peaks=peak2, threshold=0.5), # optional
            list(counts=rc2.1, peaks=peak2[signalValue_IDR>0,], threshold=c(0.1,0.8)), # optional
            list(counts=rc2.2, peaks=peak2[signalValue_IDR>0,], threshold=c(0.1,0.8))) # optional

# parallel = FALSE if no parallel backend registered
variants = scoreAllelicEffect(variants, epg1, epg2, c('HepG2','K562'), window=15, parallel=TRUE)
```

### (C) RBP-SNV impact on gene expression

```{r}
epg1 = list(deg1[,.(gene,q_value)], # first element must be DE
           list(peaks=peak1, threshold=1)) 
epg2 = list(deg2[,.(gene,q_value)], 
           list(peaks=peak2, threshold=1))
# parallel = FALSE if no parallel backend registered
variants = scoreVarImpactOnGE(variants, epg1, epg2, c('HepG2','K562'))
```
