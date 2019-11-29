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
#from sklearn.model_selection import ShuffleSplit
from keras import backend as K
import matplotlib.pyplot as plt

from tensorflow.keras.metrics import Recall

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

def recall_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    possible_positives = K.sum(K.round(K.clip(y_true, 0, 1)))
    recall = true_positives / (possible_positives + K.epsilon())
    return recall

def precision_m(y_true, y_pred):
    true_positives = K.sum(K.round(K.clip(y_true * y_pred, 0, 1)))
    predicted_positives = K.sum(K.round(K.clip(y_pred, 0, 1)))
    precision = true_positives / (predicted_positives + K.epsilon())
    return precision

    
def modelo_CNN(df,typemod):
    ''' Modelo basado en Red neuronal''' 
    X = df.copy()
    X=X.drop(columns=['causal'])
    y=df['causal'].copy()
   
    
    model = Sequential()
    model.add(Reshape((1,X.shape[1],1)))
    
    model.add(Dense(64, input_dim=64, activation='relu'))
    model.add(Conv2D(filters = 64, kernel_size = (1,5),padding = 'Same',
             activation ='relu', input_shape = (1,X.shape[1],1)))
    model.add(Conv2D(filters = 32, kernel_size = (1,5),padding = 'Same',
             activation ='relu', input_shape = (1,X.shape[1],1)))
    model.add(Conv2D(filters = 16, kernel_size = (1,5),padding = 'Same',
             activation ='relu', input_shape = (1,X.shape[1],1)))  
    model.add(Conv2D(filters = 8, kernel_size = (1,5),padding = 'Same',
             activation ='relu', input_shape = (1,X.shape[1],1)))    
    model.add(Conv2D(filters = 1, kernel_size = (1,5),padding = 'Same',
             activation ='relu', input_shape = (1,X.shape[1],1)))   
    model.add(MaxPooling2D(pool_size = (1,6), strides=(1,2)))  

    model.add(Flatten())

    model.add(Dense (1, activation='sigmoid'))

    model.compile(loss='binary_crossentropy', optimizer='adam',
              metrics=['binary_accuracy',precision_m, recall_m])
        
    history=model.fit(np.array(X), np.array(y), nb_epoch=2, 
          validation_data=(np.array(X), np.array(y)))

    ##!@@ plt_CNN(history,np.array(X), np.array(y))
    # Guardamos el modelo
    filename = '../MLmodel/'+typemod+'_CNN-model.sav'
    pickle.dump(model, open(filename, 'wb')) 

    return model


def metricas_CNN(df,model):
    ''' Metricas para testear la Red neuronal''' 
    long=int((len(df)/10)*9)
    fin=len(df)-long

    df = df.sample(frac=1).reset_index(drop=True)
    y_train=df['causal'].copy()
    y_train=y_train.head(long)
    X_train = df.copy()
    X_train=X_train.drop(columns=['causal'])
    X_train=X_train.head(long)

    y_test=df['causal'].copy()
    y_test=y_test.tail(fin)
    X_test = df.copy()
    X_test=X_test.drop(columns=['causal'])
    X_test=X_test.tail(fin)
    

    loss_train, accuracy_train,  precision_train, recall_train = model.evaluate(X_train, y_train, verbose=0)   
    loss_test, accuracy_test,  precision_test, recall_test = model.evaluate(X_test, y_test, verbose=0)
    
    print('Recall train %.3f and test %.3f:'% (recall_train, recall_test))
    print('Accuracy train %.3f  and test %.3f:'% (accuracy_train, accuracy_test))
    print('Precision train %.3f  and test %.3f:'% (precision_train, precision_test))



def plt_CNN(history,X,y):
    
  #  history=model.fit(X,y, nb_epoch=2, 
   #       validation_data=(X,y))
    
    plt.plot(history.history['recall_m'])
    plt.plot(history.history['val_recall_m'])
    plt.title('Model recall')
    plt.ylabel('Recall')
    plt.xlabel('Epoch')
    plt.legend(['Train', 'Test'], loc='upper left')
    plt.show()   







