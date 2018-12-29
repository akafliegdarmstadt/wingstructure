import numpy as np

def nonlinearLL(a, c_li, y_li, airfoil_li, v_inf, polars, max_iter=200):
    """ Nichtlineares Iterationsverfahren zur Bestimmung der Auftriebsverteilung.
    Erlaubt die direkte Verwendung von Profilpolaren.
    :param a: Anstellwinkel
    :param c_li: Profiltiefen
    :param y_li: Spannweitenpositionen
    :param airfoil_li: Profile
    :param v_inf: Anströmgeschwindigkeit
    :param polars: Polaren für die angegebenen Profile
    :param max_iter: Maximale Anzahl an Iterationen
    :return: """
      
    # Geometriegrößen bestimmen  
    S = PlaneGeom2.calcWingSRef(c_li, y_li)    # Flügelfläche
    b = PlaneGeom2.calcWingBRef(y_li)         # Spannweite
    AR = b**2/S                          # Streckung

    # Stützstellen auf Spannweite verteilt
    y_ar = np.array(y_li)

    # Stützstellen Profiltiefe
    c_ar = np.array(c_li)
    
    # Elliptische Auftriebsverteilung als Startwert für Zirkulationsverteilung

    theta_0 = 0  # werden hier zu hohe Werte definiert werden seltsame Lösungen gefunden
    
    theta_ar = theta_0*np.sqrt(1-(2*y_ar/b)**2)
    
    # Durchführen der Iterationen
    
    convergence = 0
    iteration = 0

    # Nachdem über 5 Iterationen Konvergenzkriterium erreicht wurde oder
    # anzahl maximaler Iterationen erreicht wurde abbrechen
    while convergence < 5 and iteration < max_iter:

        theta_neu = np.zeros(theta_ar.shape)

        # Schleife über Spannweitenpositionen zur Bestimmung der Ableitung
        dThetady = np.zeros(y_ar.shape)
        for jj, y2 in enumerate(y_ar):
            if jj == 0:
                # Vorwärtsdifferenzenquotient
                dThetady[jj] = (theta_ar[1]-theta_ar[0])/(y_ar[1]-y_ar[0])
            elif jj == len(y_ar)-1:
                # Rückwärtsdifferenzenquotient
                dThetady[jj] = (theta_ar[-1]-theta_ar[-2])/(y_ar[-1]-y_ar[-2])
            else:
                # zentraler Differenzenquotient 1. Ordnung
                dThetady[jj] = (theta_ar[jj+1]-theta_ar[jj-1])/(y_ar[jj+1]-y_ar[jj-1])

        # Schleife über Spannweitenpositionen
        for ii, y in enumerate(y_ar):
            # Berechnung der induzierten Anstellwinkel

            deltay = y-y_ar

            if ii == 0:
                deltay[0] = deltay[1]
            elif ii == len(y_ar)-1:
                deltay[-1] = deltay[-2]
            else:
                deltay[ii] = 0.5*(deltay[ii-1]+deltay[ii+1])

            integrand = dThetady/deltay

            if ii == 0:
                integrand[0] = integrand[1]
            elif ii == len(y_ar)-1:
                integrand[-1] = integrand[-2]
            else:
                integrand[ii] = 0.5*(integrand[ii-1]+integrand[ii+1])

            a_ind = 1 / (4 * np.pi * v_inf) * integrate.simps(integrand, y_ar)
            
            a_eff = a - a_ind

            #TODO: Re berechnen
            re = 0

            c_l,c_w = polars.interpolate(airfoil_li[ii],re,a_eff)

            theta_neu[ii] = 0.5 * v_inf * c_li[ii] * c_l

        # Randbedingungen setzen (ohne geht nix..)
        theta_neu[0] = theta_neu[-1] = 0


        theta_d = 0.05*(theta_neu - theta_ar)

        if all(theta_d < 1e-5):
            convergence += 1
        elif convergence > 1:
            convergence = 0

        theta_ar += theta_d

        iteration += 1

        print('Iteration {}, Convergence {}'.format(iteration, convergence))

    c_l = 2 / v_inf * theta_ar/c_ar

    c_L = 2*integrate.simps(c_ar*c_l,y_ar)/S

    return {'Theta':theta_ar,'C_L':c_L,'c_l':c_l}
