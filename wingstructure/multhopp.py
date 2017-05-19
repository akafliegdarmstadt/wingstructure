# -*- coding: utf-8 -*-
import numpy as np
import scipy.optimize as optimize


def multhopp(alpha, c_li, y_li, N_M=None, dcl=2*np.pi):
    alpha, N_M, v_ar, theta_ar, y_ar, chord_ar, dcl, b, AR = prepare_multhopp(alpha, c_li, y_li, N_M, dcl)
    return solve_multhopp(alpha, y_ar, chord_ar, dcl, b, AR)


def prepare_multhopp(alpha, c_li, y_li, N_M, dcl):
    """Berechnet Auftriebsverteilung nach dem Multhopp-Verfahren.
    
    :param alpha: Anstellwinkel der Fläche
    :type alpha: float
    :param c_li: Liste von Flügeltiefen
    :type c_li: list
    :param y_li: Liste von Stützstellen zu Flügeltiefen
    :type y_li: list
    :param N_M: Anzahl von Stützstellen bei Berechnung
    :type: int
    :param dcl: Auftriebsanstieg Profil
    :type dcl: float
    :rtype: dict
    """

    # Spannweite bestimmen
    b = 2*max(y_li)

    # Flügelfläche bestimmen
    S = 2 * np.trapz(y=c_li, x=y_li)

    # Streckung bestimmen
    AR = b ** 2 / S

    # Stützstellenanzahl
    if N_M is None:
        N_M = int(round(AR)*4-1)  # muss ungerade sein, sollte 4*Streckung nicht überschreiten

    # Vektor mit Indizes der Berechnungspunkte
    v_ar = np.arange(1,N_M+1)

    # Thetas der Stützstellen
    theta_ar = v_ar * np.pi / (N_M + 1)

    # Vektor der Stützstellen definieren
    y_ar = b/2 * np.cos(theta_ar)
    
    # Flügeltiefen
    chord_ar = np.interp(np.abs(y_ar), y_li, c_li)

    if not isinstance(dcl,np.ndarray):
        dcl *= np.ones(N_M)
    elif len(dcl) != N_M:
        dcl = np.interp(y_ar,y_li,dcl)

    # Anstellwinkel an den Stützstellen
    if not isinstance(alpha, np.ndarray) or len(alpha)==1:
        alpha = np.ones(N_M)*alpha
    elif len(alpha) == len(y_li):
        # interpolieren
        alpha = np.interp(y_ar, y_li, alpha, left=alpha[0], right=0)
    elif len(alpha) == len(y_ar):
        pass
    else:
        print('error') #TODO gescheites Fehlerhandling

    return alpha, N_M, v_ar, theta_ar, y_ar, chord_ar, dcl, b, AR


def solve_multhopp(alpha, y_ar, chord_ar, dcl, b, AR):

    N_M = len(y_ar)

    v_ar = np.arange(1, N_M + 1)

    # Thetas der Stützstellen
    theta_ar = v_ar * np.pi / (N_M + 1)

    # leeres N_MxN_M Array (Matrix) für Multhoppkoeffizienten
    B = np.zeros((N_M, N_M))

        # Berechnung der Multhopp Matrix / Multhoppkoeffizienten
    for v,y_v,theta_v,c,dcl_v in zip(v_ar,y_ar,theta_ar,chord_ar,dcl):
        for n,theta_n in zip(v_ar,theta_ar):
         
            # Diagonalelemente
            if(v==n):
                B[v-1,v-1]= (N_M+1)/(4*np.sin(theta_v))+2*b/(dcl_v*c)
            # übrige Elemente
            else:
                B[v-1,n-1]=-((1-(-1.)**(v-n))/2*(np.sin(theta_n)/((N_M+1)*(np.cos(theta_n)-np.cos(theta_v))**2)))

    # Berechnung der lokalen Zirkulation
    gamma_ar = np.dot(np.linalg.inv(B),alpha)

    # Gesamtauftriebsbeiwert
    C_A = np.pi*AR/(N_M+1)*sum(gamma_ar*np.sin(v_ar*np.pi/(N_M+1)))
    
    # Induzierten Anstellwinkel berechne
    a_ind = np.zeros(N_M)
    
    for v in range(N_M):
        part1 = (N_M+1)/(4*np.sin((v+1)*np.pi/(N_M+1)))*gamma_ar[v]
        
        part2 = 0
        
        for j in range(N_M):
            if j == v:
                continue
            part2 += B[v,j]*gamma_ar[j]
            
        a_ind[v] = part1 + part2 # Rätselhafterweise hier ein plus statt minus

    # Induzierten Widerstand bestimmen
    C_Wi = np.pi*AR/(N_M+1) * sum(gamma_ar*a_ind*np.sin(v_ar*np.pi/(N_M+1)))
    
    # Oswald-Faktor
    k = C_A**2/(np.pi*AR*C_Wi)
    
    # lokale Auftriebsbeiwerte
    c_a_li = 2*b/(chord_ar)*gamma_ar

    a_eff = alpha-a_ind

    return {'gamma_li': gamma_ar, 'y_v_li': y_ar, 'c_a_li': c_a_li, 'C_A':C_A,'a_ind': a_ind,
            'a_eff': a_eff,'C_Wi': C_Wi,'k': k,'AR': AR,'chord_li': chord_ar}


def inverse_multhopp(alpha_geo, C_L, c_li, y_li, N_M=None, dcl=2*np.pi):

    def fun(alpha):
        alpha = alpha_geo+alpha
        result = multhopp(alpha, c_li, y_li, N_M,dcl)
        return abs(C_L-result['C_A'])
    
    res = optimize.fsolve(fun,1/180*np.pi)
    
    res2 = multhopp(res[0],c_li,y_li,N_M,dcl)
    
    res2['alpha'] = res[0]
    
    return res2
