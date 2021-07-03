# import pymks
import numpy as np
from PIL import Image
# from pymks.stats import autocorrelate
# from pymks import PrimitiveBasis
# from pymks.tools import draw_microstructures
# from pymks.tools import draw_autocorrelations
# from pymks.stats import crosscorrelate
# from pymks.tools import draw_crosscorrelations
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from Scripts.SpatialCorrelations import *


def image_preprocessing(image_as_numpy):
    img = image_as_numpy-0.5
    img = gaussian_filter(img, sigma=4)
    return img


def image_preprocessing2(image_as_numpy):
    img = image_as_numpy-0.5
    for i in range(1000):
        S = img/(np.sqrt(1+img*img))
        G = 1-gradientMagnitude(img)
        deltaTau = 0.01
        img = img+deltaTau*S*G
    
    return img/20

def regionsOfBoundary(levelSetImage):
    '''
    Input : Level Set image
    
    Output : Only boundaries
    '''
    contourMap = np.zeros(levelSetImage.shape)
    for i in range(1,levelSetImage.shape[0]-1):
        for j in range(1,levelSetImage.shape[1]-1):
            if levelSetImage[i][j]<=0:
                if levelSetImage[i+1][j]>0 or levelSetImage[i-1][j]>0 or levelSetImage[i][j+1]>0 or levelSetImage[i][j-1]>0 or levelSetImage[i+1][j+1]>0 or levelSetImage[i-1][j-1]>0 or levelSetImage[i+1][j-1]>0 or levelSetImage[i-1][j+1]>0:
                    contourMap[i][j]=1
            
    return contourMap

def contouring(levelSetImage):
    
    contourMap = np.zeros(levelSetImage.shape)
    for i in range(1,levelSetImage.shape[0]-1):
        for j in range(1,levelSetImage.shape[1]-1):
            if levelSetImage[i][j]<=0:
                if levelSetImage[i+1][j]>0 or levelSetImage[i-1][j]>0 or levelSetImage[i][j+1]>0 or levelSetImage[i][j-1]>0 or levelSetImage[i+1][j+1]>0 or levelSetImage[i-1][j-1]>0 or levelSetImage[i+1][j-1]>0 or levelSetImage[i-1][j+1]>0:
                    contourMap[i][j]=1
            
    return contourMap
                    


def gradientMagnitude(img_):
    
    '''
    Input : numpy image 
    
    Output : Magnitude of gradient at each point
    
    '''
    gradients = np.gradient(img_)
    x_grad = gradients[0]
    y_grad = gradients[1]
    gradient_abs = np.sqrt(x_grad*x_grad+y_grad*y_grad)
    return gradient_abs

def gradientDirection(img_):
    
    '''
    Input : image
    
    Output : gradient unit vector
    
    '''
    gradients = np.gradient(img_)
    x_grad = gradients[0]
    y_grad = gradients[1]
    gradient_abs = np.sqrt(x_grad*x_grad+y_grad*y_grad)
    a = (x_grad/gradient_abs)
    b = (y_grad/gradient_abs)
    return a, b

def dphidt(listOffiles, fileNumber, timestepForDt):
    '''
    Input : String list of files in order, file you want derivative for, time step for dt
    
    Output : dpho by dt of that image
    
    '''
    
    path = listOffiles[fileNumber]
    path_plus = listOffiles[fileNumber+timestepForDt]
    path_minus = listOffiles[fileNumber-timestepForDt]
    path_plus2 = listOffiles[fileNumber+2*timestepForDt]
    path_minus2 = listOffiles[fileNumber-2*timestepForDt]
    
    img_ = dat_to_numpy(path)
    img_ =image_preprocessing(img_)
    
    img_plus = dat_to_numpy(path_plus)
    img_plus =image_preprocessing(img_plus)
    
    
    img_minus = dat_to_numpy(path_minus)
    img_minus =image_preprocessing(img_minus)
    
    img_plus2 = dat_to_numpy(path_plus2)
    img_plus2 =image_preprocessing(img_plus2)    
    
    img_minus2 = dat_to_numpy(path_minus2)
    img_minus2 =image_preprocessing(img_minus2)    
    
    return (-img_plus2+8*img_plus-8*img_minus+img_minus2)/(12*timestepForDt)
    
    
    
def velocityMagnitude(listOffiles, fileNumber, timestepForDt):
    '''
    Input : String list of files in order, file you want derivative for, time step for dt

    Output : Velocity in boundary regions
    '''
    
    path = listOffiles[fileNumber]
    img_ = dat_to_numpy(path)
    img_ = image_preprocessing(img_)
    contour = regionsOfBoundary(img_)

    gradient_magnitude =gradientMagnitude(img_)
    dphidt_ =dphidt(listOffiles, fileNumber, timestepForDt)
    
    
    velocity = (np.divide(dphidt_,gradient_magnitude))
    velocity = np.nan_to_num(velocity)
    
    return velocity*contour

def velocityDirection(listOffiles, fileNumber, timestepForDt):
    
    '''
    Input : image you need velocity unit direction for
    
    Output : velocity unit vector
    
    '''
    path = listOffiles[fileNumber]
    img_ = dat_to_numpy(path)
    img_processed = image_preprocessing(img_)
    
    normalVector = gradientDirection(img_processed)
    contour = regionsOfBoundary(img_processed)
    
    
    normalVector[0][np.isnan(normalVector[0])] = 0
    normalVector[1][np.isnan(normalVector[1])] = 0

    return normalVector[0]*contour*-1, normalVector[1]*contour*-1
    


def radiusCircular(listOffiles, fileNumber, timestepForDt):
    '''
    Input : String list of files in order, file you want derivative for, time step for dt

    Output : Velocity in boundary regions
    '''
    
    path = listOffiles[fileNumber]
    img_ = dat_to_numpy(path)
   
    r = np.sum(img_[256]+img_[257])/4
    r1= np.sum(img_[:][256]+img_[:][257])/4
    r_avg = (r+r1)/2
    return r_avg

  
def velocityCircular(listOffiles, fileNumber, timestepForDt):
    '''
    Input : String list of files in order, file you want derivative for, time step for dt

    Output : Velocity in boundary regions
    '''
    
    r_plus = radiusCircular(listOffiles, fileNumber+timestepForDt, timestepForDt)
    r_minus = radiusCircular(listOffiles, fileNumber-timestepForDt, timestepForDt)
    r_plus2 = radiusCircular(listOffiles, fileNumber+2*timestepForDt, timestepForDt)
    r_minus2 = radiusCircular(listOffiles, fileNumber-2*timestepForDt, timestepForDt)
    
    velocity = (-r_plus2+8*r_plus-8*r_minus+r_minus2)/(12*timestepForDt)
    print('velocity in pixels per frame = ',velocity)
    return velocity

