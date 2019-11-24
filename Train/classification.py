#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: natalia castejon fernandez
"""

## --   Modulos -- ##

import random
from sklearn.naive_bayes import BernoulliNB
import pickle                                                       
from sklearn.model_selection import cross_validate          
from sklearn.metrics import make_scorer,accuracy_score, precision_score, recall_score
from sklearn.tree import DecisionTreeClasifier
from sklearn.model_selection import GridSearchCV


## --  Parametros globales -- ##

random.seed(7)
scoring = {'accuracy' : make_scorer(accuracy_score), 
           'precision' : make_scorer(precision_score),
           'recall' : make_scorer(recall_score)}
    
## --  Funciones -- ##
                                                              
def modelo_bernouilli(percent,df,typeMod):
    ''' Modelo naive, basado en una bernoilli que usa como probabilidad el porcentaje de exitos del dataset de entrenamiento  ''' 
    clf = BernoulliNB()
    
    X = df.copy()
    X=X.drop(columns=['causal'])
    y=df['causal'].copy()
    
    clf.fit(X, y,sample_weight=percent)
    
    # Guardamos el modelo
    filename = '../MLmodel/'+ typeMod +'naive-model.sav'
    pickle.dump(clf, open(filename, 'wb')) 
    
    # Hacemos cross validation para conocer las metricas resultantes
    result=cross_validate(clf,X,y,cv=10,scoring=scoring)
    # Devolvemos el modelo entrenado,el valor del coef de regresion, el valor de MSE y el de RMSE 
    for key in result.keys():
        result[key]=result[key].mean()

    return clf, result
    
                                          
def modelo_DecTree(df,typeMod):
    ''' Modelo basado en un decision tree ''' 
    DTclf=DecisionTreeClasifier()
    
    # Guardamos el modelo
    filename = '../MLmodel/'+ typeMod +'DecTree-model.sav'
    pickle.dump(DTclf, open(filename, 'wb')) 

    # Hacemos GridSearch para conocer las metricas resultantes
    result=GridSearchCV
    



