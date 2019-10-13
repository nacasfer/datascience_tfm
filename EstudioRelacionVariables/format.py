#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 14:47:58 2019
@author: natalia castejon fernandez
"""

## --   Modulos -- ##
import pandas as pd
import sys



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
      #  return 3.0 
        print ('ERROR: Function formateo_min; value out of range')
        sys.exit()       

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
     #   return 3.0 
        print ('ERROR: Function formateo_sift; value out of range')
        sys.exit() 

def categoria_por_valor (linea,dic):
    if str(linea) in dic.keys():
        return dic[linea]
    elif '_' in str(linea):       
        lista=linea.replace(';','*').replace('_','*').split('*')     
    elif ';' in str(linea):      
        lista=linea.replace(';','*').replace('_','*').split('*')
    elif str(linea) == 'nan':
        return float('NaN')        
    else:
    #   return 7.
        print ('ERROR: Function categoria_por_valor; value out of range')
        sys.exit()   
    return max(list(map(dic.get, lista)))

####################
## -- int main -- ##
####################

###
### -- Lectura del archivo de variantes
variantes_DF=pd.read_csv(path , sep='\t', 
                         header=header,dtype={'causal': 'bool','MutationTaster_score':'float','SIFT_score':'float','Polyphen2_HDIV_score':'float',
                                              'Polyphen2_HVAR_score':'float','MutationTaster_score':'float','PROVEAN_score':'float',
                                              'CADD_phred':'float','phyloP20way_mammalian':'float','SiPhy_29way_logOdds':'float',
                                              'FATHMM_score':'float','gnomAD_genome_ALL':'float','1000g2015aug_all':'float','1000G_ALL':'float',
                                              'ExAC_ALL':'float','AF':'float','AF_male':'float','AF_female':'float','PopFreqMax':'float','alelos':'float'},
                                              decimal=',', usecols=fields,low_memory=False)


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
variantes_DF=variantes_DF.rename(columns={ 0:'ljb23_sift'})
variantes_DF['ljb23_sift']=variantes_DF['ljb23_sift'].astype(float)

# Eliminamos los | del principio y el final y celda a celda buscamos formato correspondiente: 
# Transformamos cadena vacía en nan y nos quedamos con el primer valor de la columna, qeu corresponde al predictor sift
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

# Transformamos valores str en numericos segun tabla de equivalencias.
# Cuando haya dos campos juntos afectando una variante, se toma el campo de mayor grado
dic1={'exonic':1.,'splicing':1.,'ncRNA':2.,'UTR5':3.,'UTR3':3.,'intronic':4.,'upstream':5.,'downstream':5.,'intergenic':6.}

ser=[]
[ser.append(categoria_por_valor(str(row),dic1))  for  row in variantes_DF['Func.refGene'] ]
variantes_DF['Func.refGene']=pd.to_numeric(ser)
del(ser)

ser=[]
[ser.append(categoria_por_valor(str(row),dic1))  for  row in variantes_DF['Func.ensGene'] ]
variantes_DF['Func.ensGene']=pd.to_numeric(ser)
del(ser)

            # print(variantes_DF[variantes_DF['Func.refGene']==7.0])   # Ctrl errores. Las filas con este valor no se han controlado. 
            # print(variantes_DF[variantes_DF['Func.ensGene']==7.0])   # Ctrl errores. Las filas con este valor no se han controlado. 

variantes_DF['Func.refGene']=variantes_DF[['Func.refGene','Func.ensGene']].max(axis=1)
del variantes_DF['Func.ensGene']

dic2={'frameshift insertion':1.,'frameshift deletion':2.,'frameshift block substitution':3.,'stopgain':4.,'stoploss':5.,'nonframeshift insertion':6.,'nonframeshift deletion':7.,'dnonframeshift block substitution':8.,'nonsynonymous SNV':9.,'synonymous SNV':10.,'unknown':11.}
      
ser=[]
[ser.append(categoria_por_valor(str(row),dic2))  for  row in variantes_DF['ExonicFunc.refGene'] ]
variantes_DF['ExonicFunc.refGene']=pd.to_numeric(ser)
del(ser)

ser=[]
[ser.append(categoria_por_valor(str(row),dic2))  for  row in variantes_DF['ExonicFunc.ensGene'] ]
variantes_DF['ExonicFunc.ensGene']=pd.to_numeric(ser)
del(ser)

variantes_DF['ExonicFunc.refGene']=variantes_DF[['ExonicFunc.refGene','ExonicFunc.ensGene']].max(axis=1)
del variantes_DF['ExonicFunc.ensGene']




#CLNSIG                    object   cambiar a ordinal
#allele_coverage           object
#function                  object   cambiar a ordinal
#grantham                  object   buscar info sobre significado de distancia grantham (max o min?)
#5000Exomes                object   buscar AMAF GMMAF EMAF
#FATHMM                    object   buscar info de significado
#clinvar                   object   cambiar a ordinal
#phylop                    object   buscar info de significado


#print(variantes_DF['phylop'].unique())

#print(variantes_DF.dtypes)

#variantes_DF[variantes_DF['causal']==1].to_csv('Pruebas.txt', sep='\t', index = False)


