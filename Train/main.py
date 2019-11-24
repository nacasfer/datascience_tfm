#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: natalia castejon fernandez
"""

## --   Modulos -- ##
import pandas as pd
import dataformat as dfmt
import estatistics as est
from imblearn.over_sampling import SMOTENC
import numpy as np
import random
import filtradoVariantes as flt
import classification as mdl
from sklearn.metrics import confusion_matrix


## --  Parametros globales -- ##

path='../../DATOS/AllCausalVariants.txt'
header=0
fields=['causal','Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','genotype','maf','gene',
        'SIFT_score','ljb23_sift','sift','Polyphen2_HDIV_score',
        'Polyphen2_HVAR_score','polyphen','genotype',
        'Func.refGene','Func.ensGene','ExonicFunc.refGene','ExonicFunc.ensGene','MutationTaster_score',
        'PROVEAN_score','CADD_phred','phyloP20way_mammalian',
        'SiPhy_29way_logOdds','FATHMM_score','CLNSIG','allele_coverage','function',
        'grantham','5000Exomes','FATHMM','clinvar','phylop','gnomAD_genome_ALL','1000g2015aug_all','1000G_ALL',
        'ExAC_ALL','AF','PopFreqMax','alelos']
      
random.seed(7)


## --  Funciones -- ##

def aplica_SMOTENC(df):
    X = df.copy()
    X=X.drop(columns=['causal'])
    y=df['causal'].copy()

    smote_nc = SMOTENC(categorical_features=[0, 2, 5, 12], random_state=0)
    X_resampled, y_resampled = smote_nc.fit_resample(X, y)
    X_resampled = pd.DataFrame(X_resampled, columns=X.columns)
    X_resampled['causal']=y_resampled
    return X_resampled


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
                                    'FATHMM_score':'float','genotype':'category',
                                    'gnomAD_genome_ALL':'float','1000g2015aug_all':'float','1000G_ALL':'float',
                                    'ExAC_ALL':'float','AF':'float','PopFreqMax':'float','alelos':'float'},
                                     decimal=',', usecols=fields,low_memory=False)


###
### -- Formateo de los datos
print ('Data format processing')
##!@@ variantes_DF=dfmt.formato_datos_origen(variantes_DF)


###
### -- Aplicacion estadisticas para completar valores nulos y reducir dimensionalidad
print('Saving model for filling NaN')
##!@@ variantes_DF=est.relleno_y_reduccion(variantes_DF)


###
### --Disminucion de la proporcion de eventos no causales (variantes == 0). (Annadimos de nuevo variantes causales, por si se descartaron con el filtro y quitamos duplicados)
#!#! print('Variants belonging to each category:\n',variantes_DF.causal.value_counts())
#!#! print('Percentage:\n',variantes_DF.causal.value_counts()[1]/((variantes_DF.causal.value_counts()[1]+variantes_DF.causal.value_counts()[0])*100))

print('Filtering Non-causal variants')
##!@@ variantes_causales=variantes_DF.copy()
##!@@ variantes_causales=variantes_causales[variantes_causales['causal']== 1]

##!@@ variantes_DF=variantes_DF[variantes_DF['causal'] == False]
##!@@ variantes_DF=flt.filtradoVariantes(variantes_DF)
##!@@ variantes_DF=pd.concat([variantes_DF, variantes_causales], axis=0)
#!#! print('Non causal variants before filter: %d \nCausal variants before filter: %d'% (variantes_DF.causal.value_counts()[0],variantes_DF.causal.value_counts()[1]))

##!@@ 


##!@@ variantes_DF.to_csv('variantesDF.txt', sep='\t', index = False)

variantes_DF=pd.read_csv('variantesDF.txt' , sep='\t', 
                             header=0,
                             dtype={'Chr':'object','Start':'int64','Ref':'object','Alt_Annovar':'object',
                                    'Alt_IR':'object','avsnp147':'object','genotype':'category','gene':'object',
                                    'causal':'bool','polyphen':'float64','Func.refGene':'category',
                                    'MutationTaster_score':'float64','allele_coverage':'object','grantham':'float64',
                                    'clinvar':'category','phylop':'float64','alelos':'float64','FRECUENCIA_t':'float64',
                                    'POBLACION_t':'float64','SIFT_t':'float64','POLYPHEN_t':'float64',
                                    'PHYLOP_t':'float64','FUNCTION_t':'category'})

###
### -- Aumento de la proporcion de eventos exito (variantes causales = 1): Sustituyendo por NaN -1
print('Upsampling')
# Eliminamos columnas tipo string
variantes_DF=variantes_DF.drop(columns=['Chr','Start','Ref','Alt_Annovar','Alt_IR','alelos','avsnp147','gene','allele_coverage'])

# Generamos set datos donde NAn = -1
#!#! print('\nNon causal variants before dealing NaN: %d \nCausal variants: %d'% (variantes_DF.causal.value_counts()[0],variantes_DF.causal.value_counts()[1]))
variantes_Nan_transformed=variantes_DF.copy()
variantes_Nan_transformed['polyphen'] = variantes_Nan_transformed['polyphen'].replace(np.nan, -1)
variantes_Nan_transformed['genotype'] = variantes_Nan_transformed['genotype'].replace(np.nan, -1)
variantes_Nan_transformed['genotype']=variantes_Nan_transformed['genotype'].astype('float64')
variantes_Nan_transformed['Func.refGene'] = variantes_Nan_transformed['Func.refGene'].replace(np.nan, -1)
variantes_Nan_transformed['Func.refGene']=variantes_Nan_transformed['Func.refGene'].astype('float64')
variantes_Nan_transformed['MutationTaster_score'] = variantes_Nan_transformed['MutationTaster_score'].replace(np.nan, -1)
variantes_Nan_transformed['grantham'] = variantes_Nan_transformed['grantham'].replace(np.nan, -1)
variantes_Nan_transformed['clinvar'] = variantes_Nan_transformed['clinvar'].replace(np.nan, -1)
variantes_Nan_transformed['clinvar']=variantes_Nan_transformed['clinvar'].astype('float64')
variantes_Nan_transformed['phylop'] = variantes_Nan_transformed['phylop'].replace(np.nan, -1)
variantes_Nan_transformed['FRECUENCIA_t'] = variantes_Nan_transformed['FRECUENCIA_t'].replace(np.nan, -1)
variantes_Nan_transformed['POBLACION_t'] = variantes_Nan_transformed['POBLACION_t'].replace(np.nan, -1)
variantes_Nan_transformed['SIFT_t'] = variantes_Nan_transformed['SIFT_t'].replace(np.nan, -1)
variantes_Nan_transformed['POLYPHEN_t'] = variantes_Nan_transformed['POLYPHEN_t'].replace(np.nan, -1)
variantes_Nan_transformed['PHYLOP_t'] = variantes_Nan_transformed['PHYLOP_t'].replace(np.nan, -1)
variantes_Nan_transformed['FUNCTION_t'] = variantes_Nan_transformed['FUNCTION_t'].replace(np.nan, -1)
variantes_Nan_transformed['FUNCTION_t']=variantes_Nan_transformed['FUNCTION_t'].astype('float64')
#!#! print('\nNon causal variants after NaN transformation: %d \nCausal variants: %d'% (variantes_Nan_transformed.causal.value_counts()[0],variantes_Nan_transformed.causal.value_counts()[1]))

# Aplicamos SMOTE
variantes_Nan_transformed_SMOTEnc=aplica_SMOTENC(variantes_Nan_transformed)
#!#! print('Numero de variantes pertenecientes a cada categoría:\n',variantes_Nan_transformed.causal.value_counts())


###
### -- Entrenamos los modelos de clasificacion

# Naive con set original (nan==-1)
print('\nTraining Naive Model with imbalanced class (original set)')
percent=variantes_DF.causal.value_counts()[1]/((variantes_DF.causal.value_counts()[1]+variantes_DF.causal.value_counts()[0])*100)
modelo, result=mdl.modelo_bernouilli(percent ,variantes_Nan_transformed,'imbalanced')

y_predict=modelo.predict(variantes_Nan_transformed.drop(columns=['causal']))
tn, fp, fn, tp = confusion_matrix(variantes_Nan_transformed['causal'], y_predict).ravel()
#!#! print ('Confusion matrix:\n TP: %d \t FN: %d \n FP: %d \t TN: %d'% (tp,fn,fp,tn) )
#!#! print('Recall:',result['test_recall'])
#!#! print('Accuracy:',result['test_accuracy'])
#!#! print('Precision:',result['test_precision'])


# Naive con set_SMOTE 
print('\nTraining Naive Model with balanced class set')                                          
modelo, result=mdl.modelo_bernouilli(0.5,variantes_Nan_transformed_SMOTEnc,'balanced')

y_predict=modelo.predict(variantes_Nan_transformed_SMOTEnc.drop(columns=['causal']))
tn, fp, fn, tp = confusion_matrix(variantes_Nan_transformed_SMOTEnc['causal'], y_predict).ravel()
#!#! print ('Confusion matrix:\n TP: %d \t FN: %d \n FP: %d \t TN: %d'% (tp,fn,fp,tn) )
#!#! print('Recall:',result['test_recall'])
#!#! print('Accuracy:',result['test_accuracy'])
#!#! print('Precision:',result['test_precision'])


# Decission Tree 
print('\nTraining Decission Tree Model with imbalanced class (original set)')  
modelo, result=mdl.modelo_DecTree(variantes_Nan_transformed,'imbalanced')

y_predict=modelo.predict(variantes_Nan_transformed.drop(columns=['causal']))
tn, fp, fn, tp = confusion_matrix(variantes_Nan_transformed['causal'], y_predict).ravel()
print ('Confusion matrix:\n TP: %d \t FN: %d \n FP: %d \t TN: %d'% (tp,fn,fp,tn) )
print('Recall:',result['test_recall'])
print('Accuracy:',result['test_accuracy'])
print('Precision:',result['test_precision'])


# Decission Tree
print('\nTraining Decission Tree Model with balanced class set')                                                        
modelo, result=mdl.modelo_DecTree(variantes_Nan_transformed_SMOTEnc,'imbalanced')

y_predict=modelo.predict(variantes_Nan_transformed_SMOTEnc.drop(columns=['causal']))
tn, fp, fn, tp = confusion_matrix(variantes_Nan_transformed_SMOTEnc['causal'], y_predict).ravel()
print ('Confusion matrix:\n TP: %d \t FN: %d \n FP: %d \t TN: %d'% (tp,fn,fp,tn) )
print('Recall:',result['test_recall'])
print('Accuracy:',result['test_accuracy'])
print('Precision:',result['test_precision'])










