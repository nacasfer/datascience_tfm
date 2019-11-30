#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: natalia castejon fernandez
"""

## --  Modulos -- ##

import pandas as pd
import pickle
import Train.dataformat as dfmt
import glob
import numpy as np
import Train.estatistics as est
import Train.filtradoVariantes as flt
import os

## --  Parametros globales -- ##

modelsforpredict= ['imbalancednaive','balancednaive','imbalancedDT','balancedDT',
                       'imbalancedGB','balancedGB']#,'imbalancedCNN','balancedCNN']

path='../DATOS/parapredecir'
outfile='../DATOS/predictedVars.csv'
header=0
fields=['Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','genotype','maf','gene',
        'SIFT_score','ljb23_sift','sift','Polyphen2_HDIV_score',
        'Polyphen2_HVAR_score','polyphen','genotype',
        'Func.refGene','Func.ensGene','ExonicFunc.refGene','ExonicFunc.ensGene','MutationTaster_score',
        'PROVEAN_score','CADD_phred','phyloP20way_mammalian',
        'SiPhy_29way_logOdds','FATHMM_score','CLNSIG','allele_coverage','function',
        'grantham','5000Exomes','FATHMM','clinvar','phylop','gnomAD_genome_ALL','1000g2015aug_all','1000G_ALL',
        'ExAC_ALL','AF','PopFreqMax','alelos']



## --  Funciones -- ##

def aplica_modelos_fillNaN(variables,df):
    principal=str(variables[0])
    for var1 in variables[1:]:
        
        filename = 'MLmodel/'+str(principal)+'-'+str(var1)+'-model.sav'
        loaded_model = pickle.load(open(filename, 'rb'))
        df[principal]=[loaded_model.predict([[df.loc[i,str(var1)]]])[0] if (str(row) in ['nan','']) & (np.isfinite(df.loc[i,str(var1)])) else est.vale_valor(row) for  i,row in enumerate(df[principal])] 
        return df[principal]
    
def print_causals(filename,df, y_predict,typemod):
    if not os.path.exists('../DATOS/PredictedCausals'):
        os.makedirs('../DATOS/PredictedCausals')
    print('Printing causal variants')
    y_predict=y_predict.astype('bool')
    df=df[y_predict]
    df.to_csv('../DATOS/PredictedCausals/'+file+'_'+typemod+'_Causalvariants.txt', sep='\t', index = False)
    

###
### -- Lectura del archivo de variantes
all_files = glob.glob(path+"/*.txt")

for filename in all_files:

    file=filename.replace(path+'/','')
    print ('\n\n#########################\nFile: ', filename)
    variantes_DF=pd.read_csv(filename , sep='\t', 
                             header=header,
                             dtype={'SIFT_score':'float','genotype':'str','Polyphen2_HDIV_score':'float','Polyphen2_HVAR_score':'float', 'MutationTaster_score':'float','PROVEAN_score':'float','CADD_phred':'float', 'phyloP20way_mammalian':'float', 'SiPhy_29way_logOdds':'float', 'FATHMM_score':'float', 'gnomAD_genome_ALL':'float','1000g2015aug_all':'float','1000G_ALL':'float', 'ExAC_ALL':'float','AF':'float','PopFreqMax':'float','alelos':'float'},
                                     decimal=',', usecols=fields,low_memory=False)
    
    ###
    ### -- Formateo de los datos
    print ('Data format processing')
    variantes_DF=dfmt.formato_datos_origen(variantes_DF)
    
    ###
    ### -- Filling NaN: Aplicacion modelos (ojo, la lista de variables debe tener el indice en el mismo orden en que fueron entrendas!!)
    print('Filling NaN')
    variantes_DF['FRECUENCIA_t']=aplica_modelos_fillNaN(['maf','5000Exomes'],variantes_DF)     
    variantes_DF['POBLACION_t']=aplica_modelos_fillNaN(['1000G_ALL','1000g2015aug_all','gnomAD_genome_ALL','AF','ExAC_ALL','PopFreqMax'],variantes_DF)    
    variantes_DF['SIFT_previo']=aplica_modelos_fillNaN(['SIFT_score','sift','ljb23_sift'],variantes_DF)    
    variantes_DF['SIFT_t']=aplica_modelos_fillNaN(['SIFT_previo','PROVEAN_score'],variantes_DF)    
    variantes_DF['POLYPHEN_t']=aplica_modelos_fillNaN(['Polyphen2_HDIV_score','Polyphen2_HVAR_score','CADD_phred'],variantes_DF)    
    variantes_DF['PHYLOP_t']=aplica_modelos_fillNaN(['phyloP20way_mammalian','FATHMM','FATHMM_score','SiPhy_29way_logOdds'],variantes_DF)    
    variantes_DF['FUNCTION_t']=aplica_modelos_fillNaN(['ExonicFunc.refGene', 'function'],variantes_DF)    

    #!#! print('Nan after and before')    
    #!#! print('FRECUENCIA_t:',pd.isna(variantes_DF['1000G_ALL']).sum(),pd.isna(variantes_DF['FRECUENCIA_t']).sum())
    #!#! print('Nan en POBLACION_t:',pd.isna(df['1000G_ALL']).sum(),pd.isna(df['POBLACION_t']).sum())
    #!#! print('Nan en SIFT_t:',pd.isna(df['SIFT_score']).sum(),pd.isna(df['SIFT_t']).sum())
    #!#! print('Nan en POLYPHEN_t:',pd.isna(df['Polyphen2_HDIV_score']).sum(),pd.isna(df['POLYPHEN_t']).sum())
    #!#! print('Nan en PHYLOP_t:',pd.isna(df['phyloP20way_mammalian']).sum(),pd.isna(df['PHYLOP_t']).sum())
    #!#! print('Nan en FUNCTION_t:',pd.isna(df['ExonicFunc.refGene']).sum(),pd.isna(df['FUNCTION_t']).sum())


    ### Eliminamos columnas sobrantes
    variantes_DF=variantes_DF.drop(columns=['maf','5000Exomes','alelos','1000G_ALL','1000g2015aug_all',
                                            'gnomAD_genome_ALL','AF','ExAC_ALL','PopFreqMax','SIFT_score',
                                            'sift','ljb23_sift','SIFT_previo', 'Polyphen2_HDIV_score',
                                            'Polyphen2_HVAR_score', 'phyloP20way_mammalian','FATHMM',
                                            'FATHMM_score','SiPhy_29way_logOdds','ExonicFunc.refGene', 
                                            'function'])


    variantes_DF['Func.refGene']=variantes_DF['Func.refGene'].astype('category')
    variantes_DF['FUNCTION_t']=variantes_DF['FUNCTION_t'].astype('category')      
    variantes_DF['clinvar']=variantes_DF['clinvar'].astype('category')

    ### Filtrado de variantes que no interesan
    ##!@@ print('Filtering Non-causal variants')
    ##!@@ variantes_DF=flt.filtradoVariantes(variantes_DF)

    ### Sustituimos nan
    variantes_Nan_transformed=dfmt.transformation_NaN(variantes_DF.copy())
       
    ###
    ### -- Predecimos
    for selectedmodel in modelsforpredict:

        if selectedmodel== 'imbalancednaive':
            print('Imbalanced naive model')
            imbalancednaive = pickle.load(open('MLmodel/imbalancednaive-model.sav' , 'rb'))    
            y_predict=imbalancednaive.predict(variantes_Nan_transformed.drop(columns=['allele_coverage','Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','gene']))
            print_causals(file,variantes_Nan_transformed,y_predict,'imbalancednaive')

        elif selectedmodel== 'balancednaive':
            print('\nBalanced naive model')    
            balancednaive = pickle.load(open('MLmodel/balancednaive-model.sav', 'rb'))
            y_predict=balancednaive.predict(variantes_Nan_transformed.drop(columns=['allele_coverage','Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','gene']))
            print_causals(file,variantes_Nan_transformed,y_predict,'balancednaive')
            
        elif selectedmodel== 'imbalancedDT':
            print('Imbalanced DT model')   
            imbalancedDT = pickle.load(open('MLmodel/imbalancedDecTree-model.sav', 'rb')) 
            y_predict=imbalancedDT.predict(variantes_Nan_transformed.drop(columns=['allele_coverage','Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','gene']))
            print_causals(file,variantes_Nan_transformed,y_predict,'imbalancedDT')
            
        elif selectedmodel== 'balancedDT':    
            print('\nBalanced DT model')  
            balancedDT = pickle.load(open('MLmodel/balancedDecTree-model.sav', 'rb'))
            y_predict=balancedDT.predict(variantes_Nan_transformed.drop(columns=['allele_coverage','Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','gene']))
            print_causals(file,variantes_Nan_transformed,y_predict,'balancedDT')
            
        elif selectedmodel== 'imbalancedGB':  
            print('Imbalanced GB model')  
            imbalancedGB = pickle.load(open('MLmodel/imbalanced_gradboost-model.sav', 'rb'))
            y_predict=imbalancedGB.predict(variantes_Nan_transformed.drop(columns=['allele_coverage','Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','gene']))
            print_causals(file,variantes_Nan_transformed,y_predict,'imbalancedGB')
            
        elif selectedmodel== 'balancedGB':     
            print('\nBalanced GB model')  
            balancedGB = pickle.load(open('MLmodel/balanced_gradboost-model.sav', 'rb'))  
            y_predict=balancedGB.predict(variantes_Nan_transformed.drop(columns=['allele_coverage','Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','gene']))
            print_causals(file,variantes_Nan_transformed,y_predict,'balancedGB')
            
        ##!@@ elif selectedmodel== 'imbalancedCNN':      
            ##!@@ print('\nImbalanced CNN model')   
            ##!@@ imbalancedCNN = pickle.load(open('MLmodel/imbalanced_CNN-model.sav', 'rb'))
            ##!@@ y_predict=imbalancedCNN.predict(variantes_Nan_transformed.drop(columns=['allele_coverage','Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','gene'])).round()
            ##!@@ print_causals(file,variantes_Nan_transformed,y_predict,'imbalancedCNN')
            
        ##!@@ elif selectedmodel== 'balancedCNN':      
            ##!@@ print('\nBalanced CNN model')  
            ##!@@ balancedCNN = pickle.load(open('MLmodel/balanced_CNN-model.sav', 'rb'))  
            ##!@@ y_predict=balancedCNN.predict(variantes_Nan_transformed.drop(columns=['allele_coverage','Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','gene'])).round()
            ##!@@ print_causals(file,variantes_Nan_transformed,y_predict,'balancedCNN')
            
        else :    
            print('Selected model not known, please train before')
    
    
    






