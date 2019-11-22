#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: natalia castejon fernandez
"""

## --  Parametros globales -- ##

fields=['Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','genotype','maf','gene',
        'causal','SIFT_score','ljb23_sift','sift','Polyphen2_HDIV_score',
        'Polyphen2_HVAR_score','polyphen',
        'Func.refGene','Func.ensGene','ExonicFunc.refGene','ExonicFunc.ensGene','MutationTaster_score',
        'PROVEAN_score','CADD_phred','phyloP20way_mammalian',
        'SiPhy_29way_logOdds','FATHMM_score','CLNSIG','allele_coverage','function',
        'grantham','5000Exomes','FATHMM','clinvar','phylop','gnomAD_genome_ALL','1000g2015aug_all','1000G_ALL',
        'ExAC_ALL','AF','PopFreqMax','alelos']


## --  Modulos -- ##
import pandas as pd
import glob


## --  Funciones -- ##

all_files = glob.glob("csv/*.txt")
li = []
for filename in all_files:
    print('Loading file',filename)
    df = pd.read_csv(filename,  header=0, sep='\t',dtype={'causal': 'bool',
                                    'SIFT_score':'float',
                                    'Polyphen2_HDIV_score':'float','Polyphen2_HVAR_score':'float',
                                    'MutationTaster_score':'float','PROVEAN_score':'float','CADD_phred':'float',
                                    'phyloP20way_mammalian':'float',
                                    'SiPhy_29way_logOdds':'float',
                                    'FATHMM_score':'float',
                                    'gnomAD_genome_ALL':'float','1000g2015aug_all':'float','1000G_ALL':'float',
                                    'ExAC_ALL':'float','AF':'float','PopFreqMax':'float','alelos':'float',
                                    'allele_coverage':'str','Alt_Annovar':'str','Alt_IR':'str','avsnp147':'str','Chr':'str','clinvar':'str','CLNSIG':'str','ExonicFunc.ensGene':'str','ExonicFunc.refGene':'str'
                                    ,'FATHMM':'str','Func.ensGene':'str','Func.refGene':'str','function':'str','gene':'str','genotype':'str','grantham':'str','ljb23_sift':'str'
                                    ,'maf':'str','phylop':'str','polyphen':'str','Ref':'str','sift':'str','Start':'str'},
                                     decimal=',', usecols=fields)
    
    li.append(df)

frame = pd.concat(li, axis=0, ignore_index=True)
print('\n\n')
print("Loaded Non causal variants: % d"% frame.causal.value_counts()[0])
print("Loaded causal variants: % d"% frame.causal.value_counts()[1])
frame = frame.drop_duplicates()

print('\nNon causal variants before drop dups: %d \nCausal variants before drop dups: %d'% (frame.causal.value_counts()[0],frame.causal.value_counts()[1]))

frame.to_csv('AllCausalVariants.txt', sep='\t', index = False)