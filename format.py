#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 14:47:58 2019

@author: natalia castejon fernandez
"""
## --   Modulos -- ##
import pandas as pd



## --  Parametros globales -- ##

path='/media/natalia/DATOS/DATA/Natalia/ML/EstudioCorrelaciVariables/TOTALVariants.txt'
header=0
fields=['Chr','Start','Ref','Alt_Annovar','Alt_IR','avsnp147','genotype','maf','gene',
        'causal','SIFT_score','ljb23_sift','sift','Polyphen2_HDIV_score',
        'Polyphen2_HVAR_score','polyphen',
        'Func.refGene','Func.ensGene','ExonicFunc.refGene','ExonicFunc.ensGene','MutationTaster_score',
        'MutationTaster_pred','PROVEAN_score','PROVEAN_pred','CADD_raw','CADD_phred','phyloP20way_mammalian',
        'SiPhy_29way_logOdds','FATHMM_score','FATHMM_pred','CLNREVSTAT','CLNSIG','allele_coverage','function',
        'grantham','5000Exomes','FATHMM','clinvar','phylop','gnomAD_genome_ALL','1000g2015aug_all','1000G_ALL',
        'ExAC_ALL','AF','AF_male','AF_female','PopFreqMax','alelos']


## --  Funciones -- ##

def formateo_min(row,c):
    if row == 'nan':
        return float('NaN') 
    elif row == '':
        return float('NaN')
    elif c not in row:
        return float(row)
    elif c in row:
       return float(min(row.split(c)))
    else:
        return 3.0       

def formateo_sift(row):
    if row == 'nan':
        return float('NaN') 
    elif row == '':
        return float('NaN')
    elif '|' not in row :
        return float(row)
    elif '|' in row:
        return float(row.split('|')[0])
    else :
        return 3.0 
        
####################
## -- int main -- ##
####################

###
### -- Lectura del archivo de variantes
variantes_DF=pd.read_csv(path , sep='\t' , header=header,dtype={'causal': 'bool',},decimal=',', usecols=fields,low_memory=False)


###
### -- Formato de columnas

# Celda a celda separamos por : , convertimos a float y escogemos el valor minimo
ser=[]
[ser.append(formateo_min(str(row),':')) for row in variantes_DF['maf'] ]            
variantes_DF['maf']=pd.to_numeric(ser)
del(ser)
            # print(variantes_DF[variantes_DF['maf']==3.0])   # Ctrl errores. Las filas con este valor no se han controlado. 

# Obtenemos el primer numero como el valor del predictor y casteamos a float
variantes_DF=variantes_DF.join(variantes_DF.pop('ljb23_sift').str.split(',', expand=True))
variantes_DF=variantes_DF.rename(columns={ 0:'ljb23_sift',  1:'resto'})
variantes_DF['ljb23_sift']=variantes_DF['ljb23_sift'].astype(float)

# Eliminamos los | del principio y el final y celda a celda buscamos formato correspondiente: 
#   Transformamos cadena vacía en nan y nos quedamos con el primer valor de la columna, qeu corresponde al predictor sift
variantes_DF['sift']=variantes_DF['sift'].str.strip('|')

ser=[]
[ser.append(formateo_sift(str(row))) for row in variantes_DF['sift'] ]            
variantes_DF['sift']=pd.to_numeric(ser)
del(ser)
            # print(variantes_DF[variantes_DF['sift']==3.0])   # Ctrl errores. Las filas con este valor no se han controlado. 

# Celda a celda separamos por | , convertimos a float y escogemos el valor minimo
variantes_DF['polyphen']=variantes_DF['polyphen'].str.strip('|')
ser=[]
[ser.append(formateo_min(str(row),'|')) for row in variantes_DF['polyphen']]    
variantes_DF['polyphen']=pd.to_numeric(ser)
del(ser)

# Transformamos valores str en numericos según hoja equivalencias.txt
#variantes_DF['Func.refGene']=transformar_equivalencia(variantes_DF['Func.refGene'])



#def transformar_equivalencia(serie):

ser=[]

import re



def cambio (linea):
    dic={'exonic':1,'splicing':1,'ncRNA':2,'UTR5':3,'UTR3':3,'intronic':4,'upstream':5,'downstream':5,'intergenic':6}
    lista=re.sub("|".join(['"', '[', ']', '\'']), "", str(linea.split(';')).split('_')).strip('[').strip(']').split(',')
    num=[]
    for elem in lista:
        num.append(dic[elem])   
    return max(num)
    

[print(cambio(str(row)))  for row in variantes_DF['Func.refGene'] ]


variantes_DF['Func.refGene']=ser
#variantes_DF['Func.refGene']=pd.to_numeric(ser)
del(ser)

print(variantes_DF['Func.refGene'].unique())

#print(variantes_DF.dtypes)

#variantes_DF[variantes_DF['causal']==1].to_csv('Pruebas.txt', sep='\t', index = False)


