# Leemos las imagenes y creamos los datos
# Para ello empleamos una funcon desarrollada por EfficientNetB3 y publicada en kaggle

import scipy
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import cv2

from tqdm import tqdm

def Dataset_loader(DIR, RESIZE, sigmaX=10):
    IMG = []
    read = lambda imname: np.asarray(Image.open(imname).convert("RGB"))
    for IMAGE_NAME in tqdm(os.listdir(DIR)):
        PATH = os.path.join(DIR,IMAGE_NAME)
        _, ftype = os.path.splitext(PATH)
        if ftype == ".jpg":
            img = read(PATH)
            img = cv2.resize(img, (RESIZE,RESIZE))
            IMG.append(np.array(img))
    return IMG

benign_train = np.array(Dataset_loader('/home/dsc/Desktop/TFM/Images/train/benign',224))
malign_train = np.array(Dataset_loader('/home/dsc/Desktop/TFM/Images/train/malignant',224))
benign_test = np.array(Dataset_loader('/home/dsc/Desktop/TFM/Images/test/benign',224))
malign_test = np.array(Dataset_loader('/home/dsc/Desktop/TFM/Images/test/malignant',224))

