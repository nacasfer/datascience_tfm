#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 14:47:58 2019
@author: natalia castejon fernandez
"""

## --   Modulos -- ##
import sys


## --  Funciones -- ##

def formateo(row,c,op):
    ''' Funcion para transformar columnas con valores numericos, 
    eliminando separadores extrannos y escogiendo el valor que proceda en cada caso:  
    El primero, el maximo o el minimo) '''
    
    if row == 'nan':
        return float('NaN') 
    elif row == '':
        return float('NaN')
    elif c not in row:
        return float(row)
    elif c in row:
        if op == 'min':
            return float(min(row.split(c)))
        elif op == 'max':
            return float(max(row.split(c)))
        elif op == 'primero':
            return float(row.split('|')[0])
        else:
            print ('ERROR: Function formateo; ilegal operand',op)
    else:
      #  return 3.0 
        print ('ERROR: Function formateo_min; value out of range')
        sys.exit()  


def categoria_por_valor  (lista,dic):
    ''' Funcion para realizar ordinal encoding a las variables con categorias ordinales '''
    
    if len(lista) ==1:
        try:
            return dic[lista[0]]
        except:
            print ('ERROR: Function categoria_por_valor; Value out of range',lista[0])
            sys.exit()
    else:
        try:
            return max(list(map(dic.get, lista)))
        except:
             print ('ERROR: Function categoria_por_valor; One or more values out of range',lista)
             sys.exit()


def formato_datos_origen(variantes_DF):
    ''' Funcion que da formato a las variables de la matriz de datos de entrada'''
    
    # Celda a celda separamos por : , convertimos a float y escogemos el valor minimo
    variantes_DF['maf']=[formateo(str(row),':','min')  for row in variantes_DF['maf'] ] 
    
                # print(variantes_DF[variantes_DF['maf']==3.0])   # Ctrl errores. Las filas con este valor no se han controlado. 
    
    # Obtenemos el primer numero como el valor del predictor y casteamos a float
    variantes_DF=variantes_DF.join(variantes_DF.pop('ljb23_sift').str.split(',', expand=True))
    variantes_DF=variantes_DF.rename(columns={ 0:'ljb23_sift'})
    del variantes_DF[ 1]; del variantes_DF[ 2]
    variantes_DF['ljb23_sift']=variantes_DF['ljb23_sift'].astype(float)
    
    # Eliminamos los | del principio y el final y celda a celda buscamos formato correspondiente: 
    # Transformamos cadena vacía en nan y nos quedamos con el primer valor de la columna, qeu corresponde al predictor sift
    variantes_DF['sift']=variantes_DF['sift'].str.strip('|')
    variantes_DF['sift']=[formateo(str(row),'|','primero')  for row in variantes_DF['sift'] ] 
    
                # print(variantes_DF[variantes_DF['sift']==3.0])   # Ctrl errores. Las filas con este valor no se han controlado. 
    
    # Celda a celda separamos por | , convertimos a float y escogemos el valor minimo
    variantes_DF['polyphen']=variantes_DF['polyphen'].str.strip('|')
    variantes_DF['polyphen']=[formateo(str(row),'|','min')  for row in variantes_DF['polyphen'] ] 
    
    
    # Transformamos cadenas str en numericos ordinales segun tabla de equivalencias. Cuando haya dos campos juntos afectando una variante se toma el campo de mayor grado
    #           Func.refGene y Func.ensGene son equivlentes. se rellena una a partir de la otra y se borra
    dic={'':float('NaN'),'nan':float('NaN'),'exonic':6.,'splicing':6.,'ncRNA':5.,'UTR5':4.,'UTR3':4.,'intronic':3.,'upstream':2.,'downstream':2.,'intergenic':1.}
    
    variantes_DF['Func.refGene']=variantes_DF['Func.refGene'].str.replace(';',',').str.replace('_',',')
    variantes_DF['Func.refGene']=[categoria_por_valor(str(row).split(',') ,dic)  for row in variantes_DF['Func.refGene'] ] 
    
    variantes_DF['Func.ensGene']=variantes_DF['Func.ensGene'].str.replace(';',',').str.replace('_',',')
    variantes_DF['Func.ensGene']=[categoria_por_valor(str(row).split(',') ,dic)  for row in variantes_DF['Func.ensGene'] ] 
    
                # print(variantes_DF[variantes_DF['Func.refGene']==7.0])   # Ctrl errores. Las filas con este valor no se han controlado. 
                # print(variantes_DF[variantes_DF['Func.ensGene']==7.0])   # Ctrl errores. Las filas con este valor no se han controlado. 
    
    variantes_DF['Func.refGene']=variantes_DF[['Func.refGene','Func.ensGene']].max(axis=1)
    del variantes_DF['Func.ensGene']
    
    #           ExonicFunc.refGene y ExonicFunc.ensGene son equivlentes. se rellena una a partir de la otra y se borra
    dic={'':float('NaN'),'nan':float('NaN'),'unknown':float('NaN'),
'synonymous SNV':1.,'nonsynonymous SNV':2.,
'nonframeshift insertion':3.,'nonframeshift deletion':3.,'nonframeshift block substitution':3.,
'frameshift insertion':4.,'frameshift deletion':4.,'frameshift block substitution':4.,
'stopgain':5.,'stoploss':5.}
    
    variantes_DF['ExonicFunc.refGene']=variantes_DF['ExonicFunc.refGene'].str.replace('_','').str.replace(':',',').str.replace(';',',').str.replace('_',',') 
    variantes_DF['ExonicFunc.refGene']=[categoria_por_valor(str(row).split(',') ,dic)  for row in variantes_DF['ExonicFunc.refGene'] ] 
    
    variantes_DF['ExonicFunc.ensGene']=variantes_DF['ExonicFunc.ensGene'].str.replace('_','').str.replace(':',',').str.replace(';',',').str.replace('_',',') 
    variantes_DF['ExonicFunc.ensGene']=[categoria_por_valor(str(row).split(',') ,dic)  for row in variantes_DF['ExonicFunc.ensGene'] ] 
    

    variantes_DF['ExonicFunc.refGene']=variantes_DF[['ExonicFunc.refGene','ExonicFunc.ensGene']].max(axis=1)
    del variantes_DF['ExonicFunc.ensGene']
    
    #           CLNSIG y clinvar son equivalentes. se rellena una a partir de la otra y se borra
    dic={'':float('NaN'),'nan':float('NaN'),'Benign':1.,'Benign/Likelybenign':2,'Likelybenign':2,
         'protective':3,
          'Uncertainsignificance':0,'Conflictinginterpretationsofpathogenicity':0,'other':0,'notprovided':0,'-':0,
          'Affects':4,'association':4,'riskfactor':5,'drugresponse':5,'Likelypathogenic':6,'Pathogenic/Likelypathogenic':7,'Pathogenic':7}
    
    variantes_DF['CLNSIG']=variantes_DF['CLNSIG'].str.replace('_','').str.replace(':',',').str.replace(' ','')
    variantes_DF['CLNSIG']=[categoria_por_valor(str(row).split(',') ,dic)  for row in variantes_DF['CLNSIG'] ] 
    
    variantes_DF['clinvar']=variantes_DF['clinvar'].str.replace('_','').str.replace(':',',').str.replace(' ','')
    variantes_DF['clinvar']=[categoria_por_valor(str(row).split(',') ,dic)  for row in variantes_DF['clinvar'] ] 
    
    variantes_DF['clinvar']=variantes_DF[['clinvar','CLNSIG']].max(axis=1)
    del variantes_DF['CLNSIG']
    
    #           function tiene su propia nomenclatura, con los campos separados por '|'
    
    dic={'':float('NaN'),'nan':float('NaN'),'unknown':float('NaN'),
'synonymous':1.,'missense':2.,
'nonframeshiftInsertion':3.,'nonframeshiftDeletion':3.,'nonframeshift':3.,'nonframeshiftBlockSubstitution':3.,
'frameshiftInsertion':4.,'frameshiftDeletion':4.,'frameshiftBlockSubstitution':4.,
'nonsense':5.,'stoploss':5.,}
    
    variantes_DF['function']=variantes_DF['function'].str.strip('|')
    variantes_DF['function']=[categoria_por_valor(str(row).split('|') ,dic)  for row in variantes_DF['function'] ] 
    
    
    # Eliminamos | y transformamos en float. Escogemos el numero mayor
    variantes_DF['grantham']=variantes_DF['grantham'].str.strip('|')
    variantes_DF['grantham']=[formateo(str(row),'|','max')  for row in variantes_DF['grantham'] ] 
    
    # En 5000genomes nos uqedamos con GMAF
    variantes_DF=variantes_DF.join(variantes_DF.pop('5000Exomes').str.strip(':').str.replace('GMAF=','').str.split(':', expand=True))
    variantes_DF=variantes_DF.rename(columns={ 1:'5000Exomes'})
    variantes_DF['5000Exomes']=variantes_DF['5000Exomes'].astype(float)
    del variantes_DF[ 0]; del variantes_DF[ 2]; del variantes_DF[ 3]; del variantes_DF[ 4]; del variantes_DF[ 5]
    
    # En FATHMM se toma el numero mayor. Por encima de 0,5 patogenico
    variantes_DF['FATHMM']=[formateo(str(row),':','max')  for row in variantes_DF['FATHMM'] ] 
    
    # En polyphen tomamos el mayor. El intervalo es [-20, 30]
    variantes_DF['phylop']=[formateo(str(row),',','max')  for row in variantes_DF['phylop'] ] 


    return variantes_DF


###### futuras implementaciones
# De allele_coverage separamos cobertura del alelo referencia y cobertura del alelo alternativo
#variantes_DF=variantes_DF.join(variantes_DF.pop('allele_coverage').str.split(',', expand=True))
#variantes_DF=variantes_DF.rename(columns={ 0:'Ref_Coverage', 1:'U', 2:'V'})


#variantes_DF['Alt_coverage']=[lambda row: row+','+str(variantes_DF['V'][i]) for i,row in enumerate (variantes_DF['U'])]
#print(variantes_DF[ 2].unique())
#print(variantes_DF.dtypes)
#print (variantes_DF['phylop'])
#variantes_DF[variantes_DF['causal']==1].to_csv('Pruebas.txt', sep='\t', index = False)
