# INCA

## Whole-genome Sequencing (WGS) Data

We extracted  the required information (Chr, Pos, Ref, Alt, Zyg) from the VCF files of WGS in HepG2 and K562 cell lines from the ENCODE project:
  + HepG2: ENCFF713BPG
  + K562: ENCFF752OAX
    
The corresponding code can be found in 'ExtractVCFInfo.py'.

Note: Zyg = 1 - Homozygous reference; 2 - Heterozygous; 3 - Unknown; 4 - Homozygous alternative
