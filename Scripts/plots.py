import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
#import seaborn as sns

def plotFromData(x,y,xlabel='xlabel',ylabel='ylabel', title='title'):
    plt.rc('font', family='serif')
    plt.rcParams["axes.edgecolor"] = "black"
    plt.rcParams["axes.linewidth"] = 1
    plt.rc('xtick', labelsize='x-small')
    plt.rc('ytick', labelsize='x-small')
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    ax.set_facecolor('white')
    plt.plot(x,y)
    
    return True

def plotFromDataOnlyY(y,xlabel='xlabel',ylabel='ylabel', title='title'):
    x=[]
    i=0
    for j in y:
        x.append(i)
        i=i+1    
    plt.rc('font', family='serif')
    plt.rcParams["axes.edgecolor"] = "black"
    plt.rcParams["axes.linewidth"] = 1
    plt.rc('xtick', labelsize='x-large')
    plt.rc('ytick', labelsize='x-large')
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    
    ax.set_facecolor('white')
    plt.plot(x,y)
    
    return True