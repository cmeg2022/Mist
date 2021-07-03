
import numpy as np
from PIL import Image
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter



def dat_to_numpy(image_path):
    
    '''
    Input : path to .dat file
    
    Output : image as a numpy array
    
    '''
    
    a = pd.read_csv(image_path,sep=" ", header=None)
    row_img =np.max(a[0])+1
    column_img = np.max(a[1])+1
    img = np.zeros((row_img,column_img))
    a_numpy = a.to_numpy()
    for i in range(a_numpy.shape[0]):
        img[int(a_numpy[i][0])][int(a_numpy[i][1])] = a_numpy[i][2]
    
    return img
def binarize_image_dat(image_as_numpy):
    
    '''
    Input : Image in numpy array format of 0 or 0 cell values
    
    Output : Binary Image as numpy
    
    '''
    return (image_as_numpy>=0.5)*1

def png_to_numpy(image_path, threshold, x_min=72, x_max=328, y_min=139, y_max=395):
    
    '''
    Input : A .png image path for B&W image
            x and y range values of image to be subsetted 
    
    Output: Numpy array of image
    
    
    '''
    
    image_in_numpy = np.array(Image.open(image_path).convert('L'))
    im_bool = image_in_numpy > threshold
    im_bin = (image_in_numpy > threshold) * 255
    
    # Subsetting to remove whitespace
    image_subset =im_bin[x_min:x_max,y_min:y_max]
    
    imagefinal = np.uint8(image_subset)
    
    
    return imagefinal

def binarize_image(image_as_numpy):
    
    '''
    Input : Image in numpy array format of 0 or 255 cell values
    
    Output : Binary Image as numpy
    
    '''
    return (image_as_numpy==255)*1
    
# def auto_corr_from_pymks(img_binary):
#     '''
#     Input : Binary image with shape (X,Y)
    
#     Output : Image of auto Correlation for both states
#     '''
#     immg = np.zeros((1,img_binary.shape[0],img_binary.shape[1]))
#     immg[0] = img_binary
#     p_basis = PrimitiveBasis(n_states=2)
#     X_auto = autocorrelate(immg, p_basis, periodic_axes=(0, 1))
#     correlations = [('black', 'black'), ('white', 'white')]
#     draw_autocorrelations(X_auto[0], autocorrelations=correlations)
    
#     return X_auto

# def cross_corr_from_pymks(img_binary):
#     '''
#     Input : Binary image with shape (X,Y)
    
#     Output : Image of cross correlation
#     '''
#     p_basis = PrimitiveBasis(n_states=2)
#     immg = np.zeros((1,img_binary.shape[0],img_binary.shape[1]))
#     immg[0] = img_binary
#     X_cross = crosscorrelate(immg, p_basis, periodic_axes=(0, 1))
#     correlations = [('black', 'white')]
#     draw_crosscorrelations(X_cross[0], crosscorrelations=correlations)
    
#     return X_cross

def probability_matrix(img_binary):
    '''
    Input : Binary Image
    
    Output : Probability of states on array
    
    
    '''
    
    m_white = np.zeros(img_binary.shape)
    m_black = np.zeros(img_binary.shape)


    for i in range(img_binary.shape[0]):
        for j in range(img_binary.shape[1]):
            ei = i+1
            ej = j+1
            wi = i-1
            wj = j-1
        
            if(ei>img_binary.shape[0]-1):
                ei = ei - img_binary.shape[0]
            if(wi<0):
                wi = wi + img_binary.shape[0]
            if(ej>img_binary.shape[1]-1):
                ej = ej - img_binary.shape[1]
            if(wj<0):
                wj = wj + img_binary.shape[1]
            
            p_white = img_binary[i][j]#+img_binary[ei][j]+img_binary[wi][j]+img_binary[i][ej]+img_binary[i][wj]
            p_white = p_white#/(5.0)
            p_black = 1.0 - p_white
        
            m_white[i][j] = p_white
            m_black[i][j] = p_black
        
    return m_white, m_black

def get_2_point_statistics(m1,m2):
    
    Fourier1 = np.fft.fft2(m1)
    Fourier2 = np.fft.fft2(m2)
    Fourier2 = Fourier2/(m1.shape[0]*m1.shape[1])
    Fourier2_conj =np.conjugate(Fourier2)
    NetFourier = np.multiply(Fourier1,Fourier2_conj)
    
    back_to_time =np.fft.ifft2(NetFourier)
    img_back = np.abs(back_to_time)
    
    image =np.reshape(img_back,(m1.shape[0],m2.shape[1]))
    return np.fft.fftshift(image) 
    
def auto_corr_from_code(img_binary):
    
    '''
    Input : Binary image with shape (X,Y)
    
    Output : Images of auto correlation
    '''
    m_white,m_black =probability_matrix(img_binary)
    auto1 = get_2_point_statistics(m_white,m_white)
    auto2 = get_2_point_statistics(m_black,m_black)
    
    
    return auto1, auto2

def cross_corr_from_code(img_binary):
    '''
    Input : Binary image with shape (X,Y)
    
    Output : Image of cross correlation
    '''
    m_white,m_black =probability_matrix(img_binary)
    auto1 = get_2_point_statistics(m_white,m_black)
    return auto1
def radialDistribution(cross):
    '''
    Input : Correlation Vector
    
    Output : probability radially distributed
    '''   
    
    radiusVector = np.zeros(cross.shape)
    for i in range(cross.shape[0]):
        for j in range(cross.shape[1]):
            radiusVector[i][j] = np.sqrt((i-cross.shape[0]/2)**2 +(j-cross.shape[1]/2)**2)
    r = np.zeros(int(cross.shape[0]/2))
    for i in range(r.shape[0]):
        r1 = (radiusVector>=i)*1
        r2 = (radiusVector<=i+1)*1
        r_net = r1*r2
        sum_val =np.sum(r_net)
        convolved = np.multiply(r_net,cross)
        r[i] =np.sum(convolved)/sum_val
    r[0] = cross[int(cross.shape[0]/2)][int(cross.shape[1]/2)]
    return r

def radialDistribution_hm(cross,radiusVector):
    '''
    Input : Correlation Vector
    
    Output : probability radially distributed hm
    '''   
    r = np.zeros(int(cross.shape[0]/2))
    for i in range(r.shape[0]):
        r1 = (radiusVector>=i)*1
        r2 = (radiusVector<=i+1)*1
        r_net = r1*r2
        convolved = np.multiply(r_net,cross)
        aa = convolved[np.nonzero(convolved)] 
        r[i] = stats.hmean(aa)
        
    r[0] = cross[int(cross.shape[0]/2)][int(cross.shape[1]/2)]
    return r

def angularDistribution(cross, theta1, theta2, angularShift = 0):
    '''
    Input : Correlation Vector, 2 angles between which you need prob distribution
    
    Output : angle/sector of consideration, probability radially distributed between angles theta1 and theta2
    '''
    start = (theta1+angularShift)*np.pi/180
    end = (theta2+angularShift)*np.pi/180
    
    
    angleVector = np.zeros(cross.shape)
    for i in range(cross.shape[0]):
        for j in range(cross.shape[1]):
            if i ==cross.shape[0]/2:
                angleVector[i][j] = 0
            if j == cross.shape[1]/2:
                angleVector[i][j] = np.pi/2
            else:
                angleVector[i][j] = np.arctan((i-cross.shape[0]/2)/(cross.shape[1]/2-j))
    
    for i in range(angleVector.shape[0]):
        for j in range(angleVector.shape[1]):
            if i < cross.shape[0]/2 and j < cross.shape[1]/2 :
                angleVector[i][j] = angleVector[i][j] + np.pi
            if i > cross.shape[0]/2 and j < cross.shape[1]/2:
                angleVector[i][j] = angleVector[i][j] + np.pi
            if i > cross.shape[0]/2 and j > cross.shape[1]/2:
                angleVector[i][j] = angleVector[i][j] + np.pi*2
    
    
    r1 = (angleVector>=start)*1
    r2 = (angleVector<=end)*1
    r_net_angle = r1*r2
    
    
    radiusVector = np.zeros(cross.shape)
    for i in range(cross.shape[0]):
        for j in range(cross.shape[1]):
            radiusVector[i][j] = np.sqrt((i-cross.shape[0]/2)**2 +(j-cross.shape[1]/2)**2)
    r = np.zeros(int(cross.shape[0]/2))
    for i in range(r.shape[0]):
        r1 = (radiusVector>=i)*1
        r2 = (radiusVector<=i+1)*1
        r_net = r1*r2
        r_net = r_net * r_net_angle
        sum_val =np.sum(r_net)
        convolved = np.multiply(r_net,cross)
        r[i] =np.sum(convolved)/sum_val
    r[0] = cross[int(cross.shape[0]/2)][int(cross.shape[1]/2)]
    
        
    return r_net_angle, r





        
