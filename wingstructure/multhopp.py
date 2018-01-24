# -*- coding: utf-8 -*-
"""
Multhopp Lifting Line Method

Reference: Konstruktionsseminar Strukturentwurf
"""
import numpy as np
import scipy.optimize as optimize


def multhopp(alpha:float, c_li: np.array, y_li: np.array, N_M:int=None, dcl:float=2*np.pi)->dict:
    """Calculates lift distribution with multhopp method.
    
    :param alpha: angle of attack in radians, either whole wing or section wise
    :param c: list or numpy.array of chord lengthes
    :param y: list of corresponding span positions
    :param N_M: number of calculation grid points in spanwise direction
    :param dcl: lift slope of airfoil or list of lift slopes
    :rtype: dict
    """
    
    # prepare calculation 
    alpha, N_M, v_ar, theta_ar, y_ar, chord_ar, dcl, b, AR = prepare_multhopp(alpha, c_li, y_li, N_M, dcl)
    
    # build up euquation system and solve it
    result = solve_multhopp(alpha, y_ar, chord_ar, dcl, b, AR)
    
    return result


def prepare_multhopp(alpha, c_li, y_li, N_M, dcl):
    """Prepares problem for multhopp calculation
    """

    # calculate wingspan
    b = 2*max(y_li)

    # calculate wing area
    S = 2 * np.trapz(y=c_li, x=y_li)

    # calculate aspect ratio
    AR = b ** 2 / S

    # number of grid points
    if N_M is None:
        N_M = int(round(AR)*4-1)  # has to be uneven, not more than 4*aspect ratio

    # grid point indices
    v_ar = np.arange(1,N_M+1)

    # create thetas
    theta_ar = v_ar * np.pi / (N_M + 1)

    # array of grid points
    y_ar = b/2 * np.cos(theta_ar)
    
    # chord depthes
    chord_ar = np.interp(np.abs(y_ar), y_li, c_li)
    
    # distribute lift slope over wing
    if not isinstance(dcl,np.ndarray):
        dcl *= np.ones(N_M)
    elif len(dcl) != N_M:
        dcl = np.interp(y_ar,y_li,dcl)

    # distribute aoa over wing
    if not isinstance(alpha, np.ndarray) or len(alpha)==1:
        alpha = np.ones(N_M)*alpha
    elif len(alpha) == len(y_li):
        # interpolate
        alpha = np.interp(y_ar, y_li, alpha, left=alpha[0], right=0)
    elif len(alpha) == len(y_ar):
        pass
    else:
        print('error') #TODO error handling

    return alpha, N_M, v_ar, theta_ar, y_ar, chord_ar, dcl, b, AR


def solve_multhopp(alpha, y_ar, chord_ar, dcl, b, AR):
    """builds up equation system and solves it"""
    
    N_M = len(y_ar)

    v_ar = np.arange(1, N_M + 1)

    # calculate thetas
    theta_ar = v_ar * np.pi / (N_M + 1)

    # create empyt matrix (N_MxN_M) for multhoppcoefficients
    B = np.zeros((N_M, N_M))

    # calculation of multhopp coefficients
    for v,y_v,theta_v,c,dcl_v in zip(v_ar,y_ar,theta_ar,chord_ar,dcl):
        for n,theta_n in zip(v_ar,theta_ar):
         
            # diagonal elements
            if(v==n):
                B[v-1,v-1]= (N_M+1)/(4*np.sin(theta_v))+2*b/(dcl_v*c)
            # non diagonal elements
            else:
                B[v-1,n-1]=-((1-(-1.)**(v-n))/2*(np.sin(theta_n)/((N_M+1)*(np.cos(theta_n)-np.cos(theta_v))**2)))

    # calculation of local circulation
    gamma_ar = np.dot(np.linalg.inv(B),alpha)

    # lift coefficient for whole wing
    C_A = np.pi*AR/(N_M+1)*sum(gamma_ar*np.sin(v_ar*np.pi/(N_M+1)))
    
    # induced angle
    a_ind = np.zeros(N_M)
    
    for v in range(N_M):
        part1 = (N_M+1)/(4*np.sin((v+1)*np.pi/(N_M+1)))*gamma_ar[v]
        
        part2 = 0
        
        for j in range(N_M):
            if j == v:
                continue
            part2 += B[v,j]*gamma_ar[j]
            
        a_ind[v] = part1 + part2 # should be a subtraction, but only works with addition

    # calculate induced drag
    C_Wi = np.pi*AR/(N_M+1) * sum(gamma_ar*a_ind*np.sin(v_ar*np.pi/(N_M+1)))
    
    # Oswald-factor
    k = C_A**2/(np.pi*AR*C_Wi)
    
    # local lift coefficient
    c_a_li = 2*b/(chord_ar)*gamma_ar

    # calculate effective aoa
    a_eff = alpha-a_ind

    return {'gamma': gamma_ar, 'y': y_ar, 'c_l': c_a_li, 'C_L':C_A,'a_ind': a_ind,
            'a_eff': a_eff,'C_Di': C_Wi,'k': k,'AR': AR,'chords': chord_ar}


def inverse_multhopp(alpha_geo, C_L, c_li, y_li, N_M=None, dcl=2*np.pi):
    """
    Calculates aoa for given lift coefficient
    :param alpha_geo: alphas given by wing geometry
    :param C_L: demanded lift coefficient
    :param c_li: chord list
    :param y_li: span position list
    :param N_M: calculation positions
    :param dcl: lift slope of airfoil[s]
    :return: resultdictionary
    """

    def fun(alpha):
        alpha = alpha_geo+alpha
        result = multhopp(alpha, c_li, y_li, N_M,dcl)
        return abs(C_L-result['C_A'])
    
    res = optimize.fsolve(fun,1/180*np.pi)
    
    res2 = multhopp(res[0],c_li,y_li,N_M,dcl)
    
    res2['alpha'] = res[0]
    
    return res2
