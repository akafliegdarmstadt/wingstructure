import numpy as np


def calcarea(outline):
    
    x_i, y_i = outline.T
    x_ip1, y_ip1 = np.roll(outline.T, 1, axis=1)
    
    A = 0.5 * np.sum(y_ip1*x_i-y_i*x_ip1)
    
    return A


def staticmoments(outline):
    
    x_i, y_i = outline.T
    x_ip1, y_ip1 = np.roll(outline.T, 1, axis=1)
    
    S_x = 1/6 * np.sum((y_i+y_ip1)*(y_ip1*x_i-y_i*x_ip1))
    S_y = 1/6 * np.sum((x_i+x_ip1)*(y_ip1*x_i-y_i*x_ip1))
    
    return S_x, S_y


def inertiamoments(outline):
    
    x_i, y_i = outline.T
    x_ip1, y_ip1 = np.roll(outline.T, 1, axis=1)
    
    I_xx = 1/12 * np.sum((y_ip1**2 + (y_i+y_ip1)*y_i)\
                         +(y_ip1*x_i-y_i*x_ip1))
    I_yy = 1/12 * np.sum((x_ip1**2 + (x_i+x_ip1)*x_i)\
                         +(y_ip1*x_i-y_i*x_ip1))
    I_xy = 1/12 * np.sum(0.5*x_ip1**2*y_i**2-0.5*x_i**2*y_ip1**2\
                          -(y_ip1*x_i-y_i*x_ip1)*(x_i*y_i+x_ip1*y_ip1))
    
    return I_xx, I_yy, I_xy


def neutralcenter(outline):
    
    area_ = calcarea(outline)
    
    S_x, S_y = staticmoments(outline)
    
    x_nc = S_y/area_
    y_nc = S_x/area_
    
    return x_nc, y_nc