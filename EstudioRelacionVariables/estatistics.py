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
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score

## --  Funciones -- ##

def estudia_coeficientes (variables,df):
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
            pearson,Ppval= stats.pearsonr(X, Y)
            dic_pearson[str(x)+','+str(y)]=pearson
            dic_covar[str(x)+','+str(y)]=np.cov(X, Y)[0,1]      
            corr, Spval = spearmanr(X, Y)
            dic_spearman[str(x)+','+str(y)]=(corr)*10
            print (str(x)+' | '+str(y),
                   "\n    Pearson: %.3f"% pearson,Ppval ,       
                   "\n    Covariate: %.3f"% np.cov(X, Y)[0,1],
                   "\n    Spearman: %.3f"% corr,Spval,)
    del df2

    
def modelo_reg_lineal(x,y):
    # Definimos el modelo lineal de prediccion
    reg=LinearRegression()
    # Entrenamos dicho modelo
    reg.fit(x,y)

    # Hacemos cross validation para conocer las metricas resultantes
    MSE=cross_val_score(reg,x,y,cv=10,scoring='neg_mean_squared_error')
    RMSE= np.sqrt(abs(MSE)).mean()
    MSE=MSE.mean()
    # Devolvemos el modelo entrenado,el valor del coef de regresion, el valor de MSE y el de RMSE
    return reg, reg.coef_, MSE, RMSE


def vale_valor(row):
    return float(row)


def modelo_varias_variables(variables,df,mod):
    principal=str(variables[0])
    #!#! print ('Nan en:',principal,pd.isna(df[principal]).sum() )
    for var1 in variables[1:]:
        df2=df[variables].copy()
        df2=df2[np.isfinite(df2[principal])]
        df2=df2[np.isfinite(df2[str(var1)])]
        modelo, reg_coef, MSE, RMSE=mod(df2[[principal]],df2[str(var1)])    
        
        #modelo, reg_coef, MSE, RMSE=modelo_reg_lineal(df2[[principal]],df2[str(var1)])    
        df[principal]=[modelo.predict([[df.loc[i,str(var1)]]])[0] if (str(row) in ['nan','']) & (np.isfinite(df.loc[i,str(var1)])) else vale_valor(row) for  i,row in enumerate(df[principal])] 
    #!#! print('Tras',var1,':', pd.isna(df[principal]).sum())   
    #!#! print ('MSE:',MSE,'\nRMSE:',RMSE)
    return df[principal]

##########################3

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


### Relacion entre variables: #Eliminamos valores nan y calculamos correlacion
# Variables de frecuencia poblacional
variables=['maf', 'gnomAD_genome_ALL', '1000g2015aug_all', '1000G_ALL', 'ExAC_ALL',
           'AF', 'PopFreqMax', '5000Exomes']
estudia_coeficientes (variables,df)

# Predictores, algunos hay que normalizarlos : 
# PROVEAN [-14, 14], CADD [0, 99], phylop [-20 , 30], grantham [5, 215]
variables=['SIFT_score', 'sift','ljb23_sift',
           'Polyphen2_HDIV_score', 'Polyphen2_HVAR_score', 'polyphen', 
           'phyloP20way_mammalian','phylop',
           'MutationTaster_score', 
           'PROVEAN_score', 
           'CADD_phred', 
           'FATHMM',
           'grantham']

df['PROVEAN_score']=(df['PROVEAN_score']+14) / 28
df['CADD_phred']=(df['CADD_phred']) / 99
df['phylop']=(df['phylop'] +20)/50
df['grantham']=(df['grantham']-5)/ 215

estudia_coeficientes (variables,df)

# 

### Generamos FRECUENCIA_t: A partir de 5000Exomes rellenamos nan de maf:
variables=['maf','5000Exomes']
#!#! print('Nan en maf:',pd.isna(df['1000G_ALL']).sum())
df['FRECUENCIA_t']=modelo_varias_variables(variables,df,modelo_reg_lineal)     
#!#! print('Nan en FRECUENCIA_t:',pd.isna(df['FRECUENCIA_t']).sum())
df.drop(['maf','5000Exomes'], axis=1)


### Generamos POBLACION_t: A partir de 1000g2015aug_all + gnomAD_genome_ALL + AF + ExAC_ALL + PopFreqMax rellenando sobre nan de 1000G_ALL
variables=['1000G_ALL','1000g2015aug_all','gnomAD_genome_ALL','AF','ExAC_ALL','PopFreqMax']
#!#! print('Nan en 1000G_ALL:',pd.isna(df['1000G_ALL']).sum())
df['POBLACION_t']=modelo_varias_variables(variables,df,modelo_reg_lineal)     
#!#! print('Nan en POBLACION_t:',pd.isna(df['POBLACION_t']).sum())
df.drop(['1000g2015aug_all','gnomAD_genome_ALL','AF','ExAC_ALL','PopFreqMax'], axis=1)

### Generamos SIFT_t: A partir de SIFT_score	sift 	 ljb23_sift 
variables=['SIFT_score','sift','ljb23_sift']
#!#! print('Nan en 1000G_ALL:',pd.isna(df['1000G_ALL']).sum())
df['SIFT_t']=modelo_varias_variables(variables,df,modelo_reg_lineal)     
#!#! print('Nan en POBLACION_t:',pd.isna(df['POBLACION_t']).sum())
df.drop(['SIFT_score','sift','ljb23_sift'], axis=1)


### Generamos POLYPHEN_t: A partir de  Polyphen2_HDIV_score  polyphen	Polyphen2_HVAR_score
variables=['Polyphen2_HDIV_score','polyphen','Polyphen2_HVAR_score']
#!#! print('Nan en 1000G_ALL:',pd.isna(df['1000G_ALL']).sum())
df['POLYPHEN_t']=modelo_varias_variables(variables,df,modelo_reg_lineal)     
#!#! print('Nan en POBLACION_t:',pd.isna(df['POBLACION_t']).sum())
df.drop(['Polyphen2_HDIV_score','polyphen','Polyphen2_HVAR_score'], axis=1)


### Generamos PHYLOP_t: A partir de  phyloP20way_mammalian	phylop + FATHMM + FATHMM_score
variables=['phylop','FATHMM','FATHMM_score']
#!#! print('Nan en 1000G_ALL:',pd.isna(df['1000G_ALL']).sum())
df['phylop']=modelo_varias_variables(variables,df,modelo_reg_lineal) 
variables=['phyloP20way_mammalian','phylop']
df['PHYLOP_t']=modelo_varias_variables(variables,df,modelo_reg_lineal) 










#antes=pd.isna(df['1000G_ALL']).sum()






