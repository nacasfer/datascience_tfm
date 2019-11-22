#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: natalia castejon fernandez
"""

## --   Modulos -- ##
import pandas as pd
import dataformat as dfmt
import estatistics as est
from imblearn.over_sampling import SMOTE
import numpy as np
import random




## --  Parametros globales -- ##

path='../../DATOS/AllCausalVariants.txt'
header=0
fields=['Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','genotype','maf','gene',
        'causal','SIFT_score','ljb23_sift','sift','Polyphen2_HDIV_score',
        'Polyphen2_HVAR_score','polyphen',
        'Func.refGene','Func.ensGene','ExonicFunc.refGene','ExonicFunc.ensGene','MutationTaster_score',
        'PROVEAN_score','CADD_phred','phyloP20way_mammalian',
        'SiPhy_29way_logOdds','FATHMM_score','CLNSIG','allele_coverage','function',
        'grantham','5000Exomes','FATHMM','clinvar','phylop','gnomAD_genome_ALL','1000g2015aug_all','1000G_ALL',
        'ExAC_ALL','AF','PopFreqMax','alelos']
      
random.seed(7)

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
print ('Data format processing')
variantes_DF=dfmt.formato_datos_origen(variantes_DF)


###
### -- Aplicacion estadisticas para completar valores nulos y reducir dimensionalidad
print('Model fiting for NaN')
variantes_DF_reducido=est.relleno_y_reduccion(variantes_DF)


###
### -- Aumento de la proporcion de eventos exito (variantes causales = 1)
#!#! print('Numero de variantes pertenecientes a cada categoría:\n',variantes_DF_reducido.causal.value_counts())
#!#! print('Porcentaje de variantes causales:\n',variantes_DF_reducido.causal.value_counts()[1]/((variantes_DF_reducido.causal.value_counts()[1]+variantes_DF_reducido.causal.value_counts()[0])*100))




# Generamos set datos donde NAn = -1
#variantes_Nan_transformed=variantes_DF_reducido.copy()
#variantes_Nan_transformed=variantes_Nan_transformed.drop(columns=['Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','gene','allele_coverage','genotype'])
#variantes_Nan_transformed = variantes_Nan_transformed.replace(np.nan, -1)


# Generamos set datos sin NAn 
#variantes_SinNan=pd.DataFrame(variantes_DF_reducido, columns=[ 'genotype' 'polyphen' 'Func.refGene' 'MutationTaster_score' 'grantham' 'clinvar' 'phylop' 'alelos' 'FRECUENCIA_t' 'POBLACION_t' 'SIFT_t' 'POLYPHEN_t' 'PHYLOP_t' 'FUNCTION_t' 'causal'])
#variantes_SinNan.dropna()




#APlicar SMOTE
#def aplica_SMOTE(df):
 #   X = df.copy()
  #  X=X.drop(columns=['causal'])
   # y=df['causal'].copy()

   # smote = SMOTE(random_state=0)
   # X_resampled, y_resampled = smote.fit_sample(X, y)

#    X_resampled = pd.DataFrame(X_resampled, columns=X.columns)
 #   X_resampled['causal']=y_resampled

  #  return X_resampled



#df=aplica_SMOTE(variantes_Nan_transformed)


#print('Numero de variantes pertenecientes a cada categoría:\n',
 #     df.causal.value_counts())





#confusion_matrix(y_test, y_pred)










