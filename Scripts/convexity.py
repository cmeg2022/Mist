import numpy as np

def monteCarloConvexity(images,number=50000):
    
    '''
    Input : image of precipitate as numpy array, number of monte carly simulations to be done
    
    Output : convexity value of precipitate
    
    '''
    
    counter = 0
    counterPositive12 = 0
    A12=[]

    pointX=[]
    pointY=[]
    for i in range(images.shape[0]):
        for j in range(images.shape[1]):
            if images[i][j]==1:
                pointX.append(i)
                pointY.append(j)
           
    for i in range(number):
        p1 = np.random.randint(0,len(pointX))
        p2 = np.random.randint(0,len(pointX))
        x1 = pointX[p1]
        y1 = pointY[p1]
        x2 = pointX[p2]
        y2 = pointY[p2] 
        x12 = int((x1+x2)/2)
        y12 = int((y1+y2)/2)

        counter=counter+1
        if images[x12][y12]==1:
            counterPositive12 = counterPositive12+1
            

    return (counterPositive12/counter)

def shortRangedAverageConvexity(images, parts=4):
    
    '''
    Input : image of microstructure as numpy array, number of parts to break image into
    
    Output : short range averaged convexity value of microstructure
    
    '''
    
    X_max = images.shape[0]
    Y_max = images.shape[1]
    x_ = int(X_max/parts)
    y_ = int(Y_max/parts)
    B = []
    img_ = (images>0.5)*1
    convexitySum=0
    
    for ii in range(parts):
        for jj in range(parts):
            img_new = img_[x_*ii:x_*(ii+1),y_*jj:y_*(jj+1)]
            convexitySum=convexitySum+monteCarloConvexity(img_new)
    
    return convexitySum/parts/parts

    
    
            
    