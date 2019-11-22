#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: natalia castejon fernandez
"""

## --   Modulos -- ##
import pandas as pd

## --  Parametros globales -- ##

value=float(0.005) # valor de frecuencia poblacional a filtrar. Se eliminan de la hoja final valores mayores
valueMAF=float(0.005) # Valor de filtrado para el maf. Se eliminan de la hoja final valores mayores



## --  Funciones -- ##


def filtradoVariantes(df):
    

    # Filtramos por frecuencia y poblacional
 #   print('Numero de variantes pertenecientes a cada categoría:\n',df.causal.value_counts())

    df=pd.concat([df[pd.isnull(df['POBLACION_t'])],df[df['POBLACION_t']<= value]],ignore_index=True)
    df=pd.concat([df[pd.isnull(df['FRECUENCIA_t'])],df[df['FRECUENCIA_t']<= valueMAF]],ignore_index=True)

  #  print('Numero de variantes pertenecientes a cada categoría:\n',df.causal.value_counts())

    # Filtramos las variables benignas y probablemente benignas