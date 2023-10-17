from cyvcf2 import VCF
import numpy as np
import gzip

## gt_types is array of 0,1,2,3==HOM_REF, HET, UNKNOWN, HOM_ALT

def extract_vcf_information(vcf_path, output): 
    vcf_file = VCF(vcf_path)

    data = []

    for variant in vcf_file:
        if len(variant.gt_types) > 1:
            print('Issues with genotypes.')
            break
            
        if len(variant.ALT) == 1:
            v = [variant.CHROM, variant.POS, variant.REF, *variant.ALT, '', variant.gt_types[0]]
        elif len(variant.ALT) == 2:
            v = [variant.CHROM, variant.POS, variant.REF, *variant.ALT, variant.gt_types[0]]
        else:
            print('Issues with ALT.')
            break
            
        data.append(v)
        
    with gzip.open(output, 'wt') as f:
        for row in data:
            f.write('\t'.join(map(str, row)) + '\n')
            
## HepG2
extract_vcf_information('https://www.encodeproject.org/files/ENCFF713BPG/@@download/ENCFF713BPG.vcf.gz',
                        '/z/Proj/pg_collaboration/RBP-analysis/WGS/HepG2_WGS.txt.gz')  
            
## K562
extract_vcf_information('https://www.encodeproject.org/files/ENCFF752OAX/@@download/ENCFF752OAX.vcf.gz',
                        '/z/Proj/pg_collaboration/RBP-analysis/WGS/K562_WGS.txt.gz')        