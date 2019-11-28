#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: natalia castejon fernandez
"""

## --   Modulos -- ##

import random
import numpy as np
from sklearn.naive_bayes import BernoulliNB
import pickle                                                       
from sklearn.model_selection import cross_validate          
from sklearn.metrics import make_scorer,accuracy_score, precision_score, recall_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import GradientBoostingClassifier
from keras.models import Sequential
from keras.layers import Dense, Conv2D,  Flatten, MaxPooling2D, Reshape



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
    # Devolvemos el modelo entrenado,y metricas de accuracy, precision y recall 

    return clf, result
    

def modelo_arbol(df,typeMod,typeData):
    ''' Modelo basado en arboles de decision. Segun parametros de entrada implementa un Decision Tree simple 
    o un Gradient Boosted''' 
    X = df.copy()
    X=X.drop(columns=['causal'])
    y=df['causal'].copy()
    
    # Hacemos GridSearch para conocer las metricas resultantes
    grid=len(list(df))+1
    
    if typeMod == 'dectree':
        result=GridSearchCV(DecisionTreeClassifier(), param_grid={'max_depth':np.arange(2,grid)},
                        cv=10,scoring='recall')
    
    elif typeMod == 'gradboost':
        result=GridSearchCV(GradientBoostingClassifier(n_estimators=20), param_grid={'max_depth':np.arange(2,grid),
                        'learning_rate':np.arange(0.1,1)},
                        cv=10,scoring='recall')
    
    result.fit(X,y)
    # Guardamos el modelo
    filename = '../MLmodel/'+ typeData +'_'+typeMod+'-model.sav'
    pickle.dump(result, open(filename, 'wb')) 

    result_fin=cross_validate(result,X,y,cv=10,scoring=scoring)
    # Devolvemos el modelo entrenado,y metricas de accuracy, precision y recall 


    return result, result_fin   


def modelo_CNN(df):
    ''' Modelo basado en Red neuronal''' 
    X = df.copy()
    X=X.drop(columns=['causal'])
    y=df['causal'].copy()
    
    model = Sequential()
    model.add(Reshape((1,X.shape[1],1)))
    model.add(Conv2D(filters = 32, kernel_size = (1,5),padding = 'Same',
             activation ='relu', input_shape = (1,X.shape[1],1)))
    model.add(Conv2D(filters = 16, kernel_size = (1,5),padding = 'Same',
             activation ='relu', input_shape = (1,X.shape[1],1)))
    model.add(MaxPooling2D(pool_size = (1,6), strides=(1,2)))
    model.add(Flatten())
    model.add(Dense (500, activation='relu'))
    model.add(Dense (1, activation='relu'))
    model.compile(loss='binary_crossentropy', optimizer='adam',
              metrics=['accuracy'])
    
    model.fit(np.array(X), np.array(y), nb_epoch=4, 
          validation_data=(np.array(X), np.array(y)))


    # Guardamos el modelo
    filename = '../MLmodel/CNN-model.sav'
    pickle.dump(model, open(filename, 'wb')) 
    
    return model




def modelo_CNN2(df):
    ''' Modelo basado en Red neuronal''' 
    X = df.copy()
    X=X.drop(columns=['causal'])
    y=df['causal'].copy()
    
    model = Sequential()
    model.add(Reshape((1,X.shape[1],1)))
    model.add(Conv2D(filters = 32, kernel_size = (1,5),padding = 'Same',
             activation ='relu', input_shape = (1,X.shape[1],1)))
    model.add(Conv2D(filters = 16, kernel_size = (1,5),padding = 'Same',
             activation ='relu', input_shape = (1,X.shape[1],1)))
    model.add(MaxPooling2D(pool_size = (1,6), strides=(1,2)))
    model.add(Flatten())
    model.add(Dense (500, activation='relu'))
    model.add(Dense (1, activation='relu'))
    model.compile(loss='binary_crossentropy', optimizer='adam',
              metrics=['accuracy'])
    
    model.fit(np.array(X), np.array(y), nb_epoch=4, 
          validation_data=(np.array(X), np.array(y)))

    result_fin=cross_validate(model(),X,y,cv=10,scoring=scoring)
    # Devolvemos el modelo entrenado,y metricas de accuracy, precision y recall 

    # Guardamos el modelo
    filename = '../MLmodel/CNN-model.sav'
    pickle.dump(model, open(filename, 'wb')) 
    
    return result_fin, result_fin 




