#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: natalia castejon fernandez
"""

## --   Modulos -- ##
import pandas as pd

## --  Parametros globales -- ##

value=float(0.007) # valor de frecuencia poblacional a filtrar. Se eliminan de la hoja final valores mayores
valueMAF=float(0.005) # Valor de filtrado para el maf. Se eliminan de la hoja final valores mayores


## --  Funciones -- ##


def filtradoVariantes(df):

    # Filtramos por frecuencia y poblacional
    df=pd.concat([df[pd.isnull(df['POBLACION_t'])],df[df['POBLACION_t']<= value]],ignore_index=True)
    #!#! print('After filtering POBLACION_t:\n',df.causal.value_counts())

    df=pd.concat([df[pd.isnull(df['FRECUENCIA_t'])],df[df['FRECUENCIA_t']<= valueMAF]],ignore_index=True)
    #!#! print('After filtering FRECUENCIA_t:\n',df.causal.value_counts())   

    # Filtramos las variables benignas y probablemente benignas (categorias 1, 2)
    df=pd.concat([df[pd.isnull(df['Func.refGene'])],df[df['Func.refGene']== 3],df[df['Func.refGene']== 4],df[df['Func.refGene']== 5],df[df['Func.refGene']== 6]],ignore_index=True)      
    #!#! print('After filtering Func.refGene:\n',df.causal.value_counts())
    
    df=pd.concat([df[pd.isnull(df['FUNCTION_t'])],df[df['FUNCTION_t']== 3],df[df['FUNCTION_t']== 4],df[df['FUNCTION_t']== 5]],ignore_index=True)   
    #!#! print('After filtering FUNCTION__t:\n',df.causal.value_counts()) 
    
    df=pd.concat([df[pd.isnull(df['clinvar'])],df[df['clinvar']== 4],df[df['clinvar']== 5],df[df['clinvar']== 6],df[df['clinvar']== 7]],ignore_index=True)       
    #!#! print('After filtering clinvar:\n',df.causal.value_counts())
    return df




