#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 14:47:58 2019
@author: natalia castejon fernandez
"""

## --   Modulos -- ##
import pandas as pd
import numpy as np
from scipy.stats import spearmanr

## --  Funciones -- ##

df=pd.read_csv('VARIANTES', sep='\t' ,header =0, 
               dtype={'Chr':'object','Start':'int64','Ref':'object','Alt_Annovar':'object',
                      'Alt_IR':'object','avsnp147':'object','genotype':'object','maf':'float64',
                      'gene':'object','causal':'bool','SIFT_score':'float64','sift':'float64',
                      'Polyphen2_HDIV_score':'float64','Polyphen2_HVAR_score':'float64','polyphen':'float64',
                      'Func.refGene':'float64','ExonicFunc.refGene':'float64','MutationTaster_score':'float64',
                      'PROVEAN_score':'float64','CADD_phred':'float64','phyloP20way_mammalian':'float64',
                      'SiPhy_29way_logOdds':'float64','FATHMM_score':'float64','allele_coverage':'object',
                      'function':'float64','grantham':'float64','FATHMM':'float64','clinvar':'float64',
                      'phylop':'float64','gnomAD_genome_ALL':'float64','1000g2015aug_all':'float64',
                      '1000G_ALL':'float64','ExAC_ALL':'float64','AF':'float64','PopFreqMax':'float64',
                      'alelos':'float64','ljb23_sift':'float64','5000Exomes':'float64'},
                      low_memory=False)

## Primero estudiamos relacion entre variables supuestamente iguales: Eliminamos valores nan y calculamos correlacion
    # maf, gnomAD_genome_ALL, 1000g2015aug_all, 1000G_ALL, ExAC_ALL, AF, PopFreqMax, 5000Exomes
variables=['maf', 'gnomAD_genome_ALL', '1000g2015aug_all', '1000G_ALL', 'ExAC_ALL', 'AF', 'PopFreqMax', '5000Exomes']
dic_pearson={}
dic_covar={}
dic_spearman={}

for i,x in enumerate (variables):
    for y in variables[i+1:]:
        df2=df.copy()
        df2=df2[np.isfinite(df2[str(x)])] 
        df2=df2[np.isfinite(df2[str(y)])]  
        X=df2[str(x)].values
        Y=df2[str(y)].values
        dic_pearson[str(x)+','+str(y)]=np.corrcoef(X, Y)[0,1]
        dic_covar[str(x)+','+str(y)]=np.cov(X, Y)[0,1]      
        corr, _ = spearmanr(X, Y)
        dic_spearman[str(x)+','+str(y)]=corr
    #   print (str(x)+' | '+str(y),
    #           "\n    Pearson: %.3f"% np.corrcoef(X, Y)[0,1],
    #           "\n    Covariate: %.3f"% np.cov(X, Y)[0,1],
    #           "\n    Spearman: %.3f"% corr)




