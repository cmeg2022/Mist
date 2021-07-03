import numpy as np
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from Scripts import SpatialCorrelations as corr
import matplotlib.pyplot as plt
import pandas as pd
import os
from Scripts import velocityCalculations as vel
import math as maths
from numpy import linalg as LA
import time

#. -----------------Helper functions for hoshen kopelman to implement compressed union find algorithm.----------
def root(arr, i):
    j = i
    try:
        while j != arr[j]:
            arr[j] = arr[arr[j]]
            j = arr[j]
    except:
        print(j, ' ', len(arr))
    return j

def find(arr, p, q):
    return root(arr,p)==root(arr,q)

def union(arr1,sz1, p, q):
    arr=arr1
    sz=sz1
    
    i = arr[p]
    j = arr[q]
    if (sz[i] < sz[j]):
        arr[i] = j
        sz[j] =sz[j]+ sz[i]
    else:
        arr[j] = i
        sz[i] =sz[i]+ sz[j]
    
    return arr, sz


#. --------------------------------------------------------------------------------------------------------------------------



def hoshenKoplemanLabels(img_):
    
    '''
    Input : microstructure as numpy array
    
    Output : array with each precipitate separated and numbered, matrix numbered as 0
    
    '''
    
    properLabels = []
    sz = []
    properLabels=[]
    sz=[]
    properLabels.append(0)
    sz.append(1)
        
    
    largestLabel = 0
    labels = np.zeros(img_.shape).astype('int32')
    for j in range(img_.shape[0]):
        for i in range(img_.shape[1]):
            if img_[j][i] == 1:
                if i>0 and j>0:
                    left = img_[j][i-1]
                    above = img_[j-1][i]
        
                    if left==0 and above== 0:
                        largestLabel=largestLabel+1
                        labels[j][i] = largestLabel
                        properLabels.append(largestLabel)
                        sz.append(1)
                    
                    if left==1 and above== 0:
                        labels[j][i] = root(properLabels,labels[j][i-1])
                    if left==0 and above== 1:
                        labels[j][i] = root(properLabels,labels[j-1][i])
                    if left==1 and above== 1:
                    
                        if labels[j][i-1]!= labels[j-1][i]:
                            if labels[j][i-1]> labels[j-1][i]:
                                properLabels, sz= union(properLabels,sz,int(labels[j][i-1]),int(labels[j-1][i]))
                                labels[j][i] = root(properLabels,int(labels[j][i-1]))
                            if labels[j][i-1]< labels[j-1][i]:
                                properLabels , sz=union(properLabels,sz,int(labels[j-1][i]),int(labels[j][i-1]))
                                labels[j][i] = root(properLabels,int(labels[j][i-1]))
                        else:
                            labels[j][i] = root(properLabels,labels[j][i-1])
                        
                if i == 0 and j == 0:
                    if img_[j][i] == 1:
                        largestLabel=largestLabel+1
                        labels[j][i] = largestLabel
                        properLabels.append(largestLabel)
                        sz.append(1)
                if i == 0 and j>0:
                    above = img_[j-1][i]
                    if img_[j][i] == 1:
                        if above==1:
                            labels[j][i] = root(properLabels,labels[j-1][i])
                        if above==0:
                            largestLabel=largestLabel+1
                            labels[j][i] = largestLabel
                            properLabels.append(largestLabel)
                            sz.append(1)
                if i>0 and j==0:
                    left = img_[j][i-1]
                    if img_[j][i] == 1:
                        if left==1:
                            labels[j][i] = root(properLabels,labels[j][i-1])
                    
                        if left==0 :
                            largestLabel=largestLabel+1
                            labels[j][i] = largestLabel
                            properLabels.append(largestLabel)
                            sz.append(1)
                            
    for i in range(labels.shape[0]):
        for j in range(labels.shape[1]):
            if labels[i][j]!=0:
                labels[i][j] = root(properLabels,labels[i][j])
                            

    #------------------ 2nd Loop starts here--------------------------#
    
    for j in range(img_.shape[0]):
        for i in range(img_.shape[1]):
            if img_[j][i] == 1:
                # part 1
                if i==0 or j==0:
         
                    if i==0:
                        iminus1 = img_.shape[1]-1
                    else:
                        iminus1 = i-1
                    if j==0:  
                        jminus1 = img_.shape[0]-1
                    else:
                        jminus1 = j-1
                    left = img_[j][iminus1]
                    above = img_[jminus1][i]
                    if left==1 and above== 0:
                        properLabels,sz=union(properLabels,sz,int(labels[j][i]),int(labels[j][iminus1]))
                        #labels[j][i] = root(properLabels,labels[j][iminus1])
                    if left==0 and above== 1:
                        properLabels,sz=union(properLabels,sz,int(labels[j][i]),int(labels[jminus1][i]))
                        #labels[j][i] = root(properLabels,labels[jminus1][i])
                    if left==1 and above== 1:
                        if labels[j][iminus1]!= labels[jminus1][i]:
                            properLabels,sz=union(properLabels,sz,int(labels[j][iminus1]),int(labels[jminus1][i]))
     
    listy ={}
    num=0
    for i in range(labels.shape[0]):
        for j in range(labels.shape[1]):
            if labels[i][j]!=0:
                labels[i][j] = root(properLabels,labels[i][j])
                
    for i in range(labels.shape[0]):
        for j in range(labels.shape[1]):
            if labels[i][j]!=0:            
                if labels[i][j] in listy.keys():
                    labels[i][j]=listy[labels[i][j]]
                else:
                    num=num+1
                    listy[labels[i][j]]=num
                    labels[i][j]=num
    return labels
    
def precipitateCentres(labelImage, labelNumber):
    
    '''
    Input : labels output from function: hoshenKoplemanLabels() , precipitate number 
    
    Output : cog along x and y
    
    '''
    img = (labelImage==labelNumber)*1
    m_y =[]
    m_x =[]
    x = []
    y = []

    for i in range(img.shape[0]):
        m_y.append(np.sum(img[i]))
        y.append(i/(img.shape[0]-1)*np.pi*2)
    for j in range(img.shape[1]):
        m_x.append(np.sum(img[:,j]))
        x.append(j/(img.shape[1]-1)*np.pi*2)
    m_x = np.array(m_x)
    m_y = np.array(m_y)
    x = np.array(x)
    y = np.array(y)
    

    alpha_x = np.cos(x)
    beta_x = np.sin(x)
    alpha_mean = np.sum(m_x*np.cos(x))/(np.sum(m_x))
    beta_mean = np.sum(m_x*np.sin(x))/(np.sum(m_x))
    theta_avg_x = np.arctan2(-1*beta_mean, -1*alpha_mean)+np.pi
    cog_x = img.shape[1]*theta_avg_x/(2*np.pi)
    
    alpha_y = np.cos(y)
    beta_y = np.sin(y)
    alpha_mean = np.sum(m_y*np.cos(y))/(np.sum(m_y))
    beta_mean = np.sum(m_y*np.sin(y))/(np.sum(m_y))
    theta_avg_y = np.arctan2(-1*beta_mean, -1*alpha_mean)+np.pi
    cog_y = img.shape[0]*theta_avg_y/(2*np.pi)
    
    
        
    return int(cog_y),int(cog_x)

def findAngleMajorMinorEigenvector(labelImage, labelNumber):
    img = (labelImage==labelNumber)*1
    if np.sum(img[0]) +np.sum(img[:,0]) + np.sum(img[img.shape[0]-1]) + np.sum(img[:,img.shape[1]-1])==0:
        t,l,b,e = angleOfInclination(labelImage, labelNumber)
    else:
        t,l,b,e = PBC_angleOfInclination(labelImage, labelNumber)
        
    return t,l,b,e
        
    
    

def angleOfInclination(labelImage, labelNumber):
    
    '''
    Input : labels output from function: hoshenKoplemanLabels() , precipitate number 
    
    Output : angle of inclination in degrees
    
    '''
    
    img = (labelImage==labelNumber)*1
    M00 = 0
    M10 = 0
    M01 = 0
    M20 = 0
    M02 = 0
    M11 = 0
    for i in range(img.shape[0]):
        for j in range(img.shape[1]):
            if img[i][j]==1:
                M00 = M00+1
                M10 = M10+j
                M01 = M01+i
                M20 = M20+j*j
                M02 = M02+i*i
                M11 = M11+i*j
    
    x_avg = M10/M00
    y_avg = M01/M00
    mu20 = M20/M00-x_avg*x_avg
    mu02 = M02/M00-y_avg*y_avg
    mu11 = M11/M00-x_avg*y_avg
        
    theta = 0.5*np.arctan(2*mu11/(mu20-mu02))
    l = np.sqrt(8*(mu02+mu20+np.sqrt(4*mu11*mu11+(mu20-mu02)*(mu20-mu02))))
    b = np.sqrt(8*(mu02+mu20-np.sqrt(4*mu11*mu11+(mu20-mu02)*(mu20-mu02))))
    
    matrix = np.array([[mu20, mu11],[mu11, mu02]])
    eigenVectors = LA.eig(matrix)
    return theta*180/np.pi, l, b, eigenVectors[1]
    

def PBC_angleOfInclination(labelImage, labelNumber):
    
    
    '''
    Input : labels output from function: hoshenKoplemanLabels() , precipitate number 
    
    Output : angle of inclination in degrees for precipitate at boundary
    
    '''
    
    img = (labelImage==labelNumber)*1
    x = img.shape[0]
    y = img.shape[1]
    img_extended = np.zeros((3*img.shape[0],3*img.shape[1]))
    img_extended[0:x,0:y] = img
    img_extended[x:2*x,0:y] = img
    img_extended[2*x:3*x,0:y] = img
    img_extended[0:x,y:2*y] = img
    img_extended[0:x,2*y:3*y] = img
    img_extended[x:2*x,y:2*y] = img
    img_extended[x:2*x,2*y:3*y] = img
    img_extended[2*x:3*x,y:2*y] = img
    img_extended[2*x:3*x,2*y:3*y] = img
    
    new_cogx, new_cogy = precipitateCentres(labelImage, labelNumber)
    new_cogx = new_cogx+x
    new_cogy = new_cogy+y
    
    labelsBig = hoshenKoplemanLabels(img_extended)
    label_for_our_ppt = labelsBig[new_cogx][new_cogy]
    theta, l, b, e = angleOfInclination(labelsBig, label_for_our_ppt)

    return theta, l, b, e
    
    
def areaDistribution(labels):
    
    
    '''
    Input : labels output from function: hoshenKoplemanLabels()
    
    Output : array with area of each pricipitate
    
    '''
    
    A =[]
    noOfLabels = np.max(labels)
    for i in range(1,noOfLabels+1):
        A.append(np.sum((labels==i)*1))
    return A

def precipitateCentresNonCog(labelImage, labelNumber):
    
    '''
    Input : labels output from function: hoshenKoplemanLabels(), precipitate number
    
    Output : centre of gravities along x and y
    
    '''
    img = (labelImage==labelNumber)*1
    m_y =[]
    m_x =[]
    x = []
    y = []

    for i in range(img.shape[0]):
        m_y.append(np.sum(img[i]))
        y.append(i/(img.shape[0]-1)*np.pi*2)
    for j in range(img.shape[1]):
        m_x.append(np.sum(img[:,j]))
        x.append(j/(img.shape[1]-1)*np.pi*2)
    m_x = np.array(m_x)
    m_y = np.array(m_y)
    x = np.array(x)
    y = np.array(y)
    

    alpha_x = np.cos(x)
    beta_x = np.sin(x)
    alpha_mean = np.sum(m_x*np.cos(x))/(np.sum(m_x))
    beta_mean = np.sum(m_x*np.sin(x))/(np.sum(m_x))
    theta_avg_x = np.arctan2(-1*beta_mean, -1*alpha_mean)+np.pi
    cog_x = img.shape[1]*theta_avg_x/(2*np.pi)
    
    alpha_y = np.cos(y)
    beta_y = np.sin(y)
    alpha_mean = np.sum(m_y*np.cos(y))/(np.sum(m_y))
    beta_mean = np.sum(m_y*np.sin(y))/(np.sum(m_y))
    theta_avg_y = np.arctan2(-1*beta_mean, -1*alpha_mean)+np.pi
    cog_y = img.shape[0]*theta_avg_y/(2*np.pi)
    
    
        
    return int(cog_y),int(cog_x)
