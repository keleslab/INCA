# INCA

## Whole-genome Sequencing (WGS) Data

We extracted  the required information (Chr, Pos, Ref, Alt, Zyg) from the VCF files of WGS in HepG2 and K562 cell lines from the ENCODE project:
  + HepG2: ENCFF713BPG
  + K562: ENCFF752OAX
    
The corresponding code can be found in 'ExtractVCFInfo.py'.

Note: Zyg = 0 - Homozygous reference; 1 - Heterozygous; 2 - Unknown; 3 - Homozygous alternative

## Main example in the paper

```{r}
directory = 'https://raw.github.com/jduan607/INCA/master/ENCODE_eCLIP'

## HepG2
peak1 = fread(file.path(directory,'PeakSignals','HNRNPK_HepG2_PeakSignals.txt.gz'))

exp1.1 = fread(file.path(directory,'NormRC','HNRNPK_HepG2_NormRC_Rep1.txt.gz')) 
exp1.2 = fread(file.path(directory,'NormRC','HNRNPK_HepG2_NormRC_Rep2.txt.gz')) 
ctrl1 = fread(file.path(directory,'NormRC','HNRNPK_HepG2_NormRC_Ctrl.txt.gz'))

rc1.1 = readCountsWRTControl(exp1.1, ctrl1)
rc1.2 = readCountsWRTControl(exp1.2, ctrl1)

## K562
peak2 = fread(file.path(directory,'PeakSignals','HNRNPK_K562_PeakSignals.txt.gz'))

exp2.1 = fread(file.path(directory,'NormRC','HNRNPK_K562_NormRC_Rep1.txt.gz')) 
exp2.2 = fread(file.path(directory,'NormRC','HNRNPK_K562_NormRC_Rep2.txt.gz')) 
ctrl2 = fread(file.path(directory,'NormRC','HNRNPK_K562_NormRC_Ctrl.txt.gz'))

rc2.1 = readCountsWRTControl(exp2.1, ctrl2)
rc2.2 = readCountsWRTControl(exp2.2, ctrl2)

dge1 = fread('https://www.encodeproject.org/files/ENCFF382MHL/@@download/ENCFF382MHL.tsv')[status=='OK',]
dge2 = fread('https://www.encodeproject.org/files/ENCFF527KUQ/@@download/ENCFF527KUQ.tsv')[status=='OK',]

## WGS
directory = 'https://raw.github.com/jduan607/INCA/master'

wgs1 = fread(file.path(directory,'ENCODE_WGS','HepG2_WGS.txt.gz'))
wgs2 = fread(file.path(directory,'ENCODE_WGS','K562_WGS.txt.gz'))
```
