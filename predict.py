#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: natalia castejon fernandez
"""



## --  Modulos -- ##

import pickle





## --  Funciones -- ##

def aplica_modelos_fillNaN(X_test,Y_test):
    filename=''
    loaded_model = pickle.load(open(filename, 'rb'))
    result = loaded_model.score(X_test, Y_test)
    print(result)