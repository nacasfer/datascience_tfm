#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 14:47:58 2019
@author: natalia castejon fernandez
"""

## --   Modulos -- ##
import pandas as pd
import dataformat as dfmt
import estatistics as est
from imblearn.over_sampling import SMOTE





## --  Parametros globales -- ##

path='../../DATOS/TOTALVariants.txt'
header=0
fields=['Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','genotype','maf','gene',
        'causal','SIFT_score','ljb23_sift','sift','Polyphen2_HDIV_score',
        'Polyphen2_HVAR_score','polyphen',
        'Func.refGene','Func.ensGene','ExonicFunc.refGene','ExonicFunc.ensGene','MutationTaster_score',
        'PROVEAN_score','CADD_phred','phyloP20way_mammalian',
        'SiPhy_29way_logOdds','FATHMM_score','CLNSIG','allele_coverage','function',
        'grantham','5000Exomes','FATHMM','clinvar','phylop','gnomAD_genome_ALL','1000g2015aug_all','1000G_ALL',
        'ExAC_ALL','AF','PopFreqMax','alelos']
      


###
### -- Lectura del archivo de variantes
variantes_DF=pd.read_csv(path , sep='\t', 
                             header=header,
                             dtype={'causal': 'bool',
                                    'SIFT_score':'float',
                                    'Polyphen2_HDIV_score':'float','Polyphen2_HVAR_score':'float',
                                    'MutationTaster_score':'float','PROVEAN_score':'float','CADD_phred':'float',
                                    'phyloP20way_mammalian':'float',
                                    'SiPhy_29way_logOdds':'float',
                                    'FATHMM_score':'float',
                                    'gnomAD_genome_ALL':'float','1000g2015aug_all':'float','1000G_ALL':'float',
                                    'ExAC_ALL':'float','AF':'float','PopFreqMax':'float','alelos':'float'},
                                     decimal=',', usecols=fields,low_memory=False)


###
### -- Formateo de los datos
variantes_DF=dfmt.formato_datos_origen(variantes_DF)


###
### -- Aplicacion estadisticas para completar valores nulos y reducir dimensionalidad
variantes_DF_reducido=est.relleno_y_reduccion(variantes_DF)


###
### -- Aumento de la proporcion de eventos exito (variantes causales = 1)
#!#! print('Numero de variantes pertenecientes a cada categoría:\n',variantes_DF_reducido.causal.value_counts())
#!#! print('Porcentaje de variantes causales:\n',variantes_DF_reducido.causal.value_counts()[1]/((variantes_DF_reducido.causal.value_counts()[1]+variantes_DF_reducido.causal.value_counts()[0])*100))


X = pd.DataFrame(variantes_DF_reducido, columns=['Chr' 'Start' 'Ref' 'Alt_Annovar' 'Alt_IR' 'avsnp147' 'genotype' 'gene' 'polyphen' 'Func.refGene' 'MutationTaster_score' 'allele_coverage' 'grantham' 'clinvar' 'phylop' 'alelos' 'FRECUENCIA_t' 'POBLACION_t' 'SIFT_t' 'POLYPHEN_t' 'PHYLOP_t' 'FUNCTION_t'])

y=variantes_DF_reducido.causal


smote = SMOTE(random_state=0)
X_resampled, y_resampled = smote.fit_sample(X, y)

X_resampled = pd.DataFrame(X_resampled, columns=X.columns)
X_resampled['causal']=y_resampled







print('Numero de variantes pertenecientes a cada categoría:\n',
      variantes_DF_reducido.causal.value_counts())





#confusion_matrix(y_test, y_pred)










