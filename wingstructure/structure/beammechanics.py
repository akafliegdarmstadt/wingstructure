import numpy as np

class Crosssection:
    # Crosssection is just for wingstructure (two cells)
    # inputs are tuples (3 arrays per tuple)
    
    def __init__(self, geometries_orig, youngsmoduli, shearmoduli, thickness):
        
        self.youngsmoduli = youngsmoduli
        self.shearmoduli = shearmoduli
        self.thickness = thickness
        
        NC = self._normal_center(geometries_orig)
        
        geometries_nc = [geometry - NC for geometry in geometries_orig]
        
        self.Θ = self._calculate_principalaxis(geometries_nc)
        
        self.geometries = [transform(geometry_nc, self.Θ) for geometry_nc
                           in geometries_nc]
        
    def _normal_center(self, geometries):
        ES_total = np.zeros(2)
        EA_total = 0.0
        for geometry, E, t in zip(geometries, self.youngsmoduli, self.thickness):
            ES_total += ES(geometry, E, t)
            EA_total += EA(geometry, E, t)
        NC = np.flip(ES_total/EA_total, 0)
        return NC
        
    def _calculate_principalaxis(self, geometries_nc):
            
        EI_total = np.zeros(3)
        
        for geometry_nc, E, t in zip(geometries_nc, self.youngsmoduli, self.thickness):
            # translation to normal center
            EI_total += EI(geometry_nc, E, t)

        [EI_y, EI_z, EI_yz] = EI_total

        θ = -0.5*np.arctan(2*EI_yz/(EI_y - EI_z))
        
        return θ
    
    def shearflow_open(self, shearforce):
        # shearflow open profile due to shearforce
        EI_total = 0
        n_xSQz = []
        n_xSQy = []
            
        EI_total = sum([EI(geometry, E, t) for geometry, E, t in
                    zip(self.geometries, self.youngsmoduli, self.thickness)])
            
        for geometry, E, t in zip(self.geometries, self.youngsmoduli, self.thickness):
            ES = ES_circulate(geometry, E, t)
            n_xSQz_i = -(shearforce[1]/EI_total[0])*ES[0]
            n_xSQy_i = -(shearforce[0]/EI_total[1])*ES[1]
            n_xSQz.append(n_xSQz_i)
            n_xSQy.append(n_xSQy_i)
            
        return n_xSQz, n_xSQy
    
    def integrate_sGt(self):
        integrate_sGt_ = []
        for geometry, G, t in zip(self.geometries, self.shearmoduli, self.thickness):
            integrate_sGt_j = 0
            for i in range(0,len(geometry)-1):
                a = np.array(geometry[i])
                b = np.array(geometry[i+1])
                Δs_i = np.linalg.norm(a-b)
                integrate_sGt_j += Δs_i/(G[i]*t[i])
            integrate_sGt_.append(integrate_sGt_j)
            
        return integrate_sGt_
    
    def integrate_nsGt(self, shearforce):
        n_xSQz, n_xSQy = self.shearflow_open(shearforce)
        integrate_nQzsGt = []
        integrate_nQysGt = []
        for geometry, G, t, n_xSQz_, n_xSQy_ in \
            zip(self.geometries, self.shearmoduli, self.thickness, n_xSQz, n_xSQy):
            integrate_nQzsGt_j = 0
            integrate_nQysGt_j = 0
            for i in range(0,len(geometry)-1):
                a = np.array(geometry[i])
                b = np.array(geometry[i+1])
                Δs_i = np.linalg.norm(a-b)
                integrate_nQzsGt_j += n_xSQz_[i]*Δs_i/(G[i]*t[i])
                integrate_nQysGt_j += n_xSQy_[i]*Δs_i/(G[i]*t[i])
            integrate_nQzsGt.append(integrate_nQzsGt_j)
            integrate_nQysGt.append(integrate_nQysGt_j)
            
        return integrate_nQzsGt, integrate_nQysGt
    
    def shearflow_cell(self, shearforce):
        # constant shearflow in cells due to shearforce
        integrate_sGt_ = self.integrate_sGt()
        integrate_nQzsGt, integrate_nQysGt = self.integrate_nsGt(shearforce)
        
        a = integrate_sGt_cell_left = integrate_sGt_[0] + integrate_sGt_[1]
        d = integrate_sGt_cell_right = integrate_sGt_[2] + integrate_sGt_[1]
        b = c = -integrate_sGt_[1]
        
        e = -(-integrate_nQzsGt[0] + integrate_nQzsGt[1])
        f = -(integrate_nQzsGt[2] - integrate_nQzsGt[1])
        
        g = -(-integrate_nQysGt[0] + integrate_nQysGt[1])
        h = -(integrate_nQysGt[2] - integrate_nQysGt[1])
        
        A = np.matrix([[a, b], [c, d]])
        y = np.array([e,f]).T
        z = np.array([g,h]).T

        n_xS0_Qz = np.linalg.solve(A, y) # left, right
        n_xS0_Qy = np.linalg.solve(A, z) # left, right
        
        return n_xS0_Qz, n_xS0_Qy

    def shearflow_shearforce(self, shearforce):
        # shearflow due to shearforce
        # start point: lower intersection left, center, right
        n_xSQz, n_xSQy = self.shearflow_open(shearforce)
        n_xS0_Qz, n_xS0_Qy = self.shearflow_cell(shearforce)
        
        n_xSQz_left = n_xSQz[0] - n_xS0_Qz[0]
        n_xSQz_center = n_xSQz[1] - n_xS0_Qz[1] + n_xS0_Qz[0]
        n_xSQz_right = n_xSQz[2] + n_xS0_Qz[1]

        n_xSQy_left = n_xSQy[0] - n_xS0_Qy[0]
        n_xSQy_center = n_xSQy[1] - n_xS0_Qy[1] + n_xS0_Qy[0]
        n_xSQy_right = n_xSQy[2] + n_xS0_Qy[1]
        
        n_xSQz = n_xSQz_left, n_xSQz_center, n_xSQz_right
        n_xSQy = n_xSQy_left, n_xSQy_center, n_xSQy_right

        return n_xSQz, n_xSQy # left, center, right
        
    def moment_shearflow(self, shearforce):
        # moment due to shearflow
        n_xSQz, n_xSQy = self.shearflow_shearforce(shearforce)
        moment_Qz = []
        moment_Qy = []
        for geometry, n_xSQz_, n_xSQy_ in \
        zip(self.geometries, n_xSQz, n_xSQy):
            moment_Qz_ = 0
            moment_Qy_ = 0
            for i in range(0,len(geometry)-1):
                a = np.array(geometry[i])
                b = np.array(geometry[i+1])
                m = (a+b)/2
                n = b-a
                distance = abs(np.cross(m,n))
                moment_Qz_ += distance*n_xSQz_[i]
                moment_Qy_ += distance*n_xSQy_[i]
                
            moment_Qz.append(moment_Qz_)
            moment_Qy.append(moment_Qy_)
        
        return moment_Qz, moment_Qy # left, center, right
    
    def shearcenter(self):
        # shear center due to shearflow referred to normal center
        shearforce = [1,1]
        moment_Qz, moment_Qy = self.moment_shearflow(shearforce)
        # transform shearforce because moment_shearflow is due to transformed shearforce
        Q_y = shearforce[0]
        Q_z = shearforce[1]
        # counterclockwise positive
        moment_Qz_total = -moment_Qz[0] - moment_Qz[1] + moment_Qz[2]
        moment_Qy_total = -moment_Qy[0] - moment_Qy[1] + moment_Qy[2]
        
        if Q_z == 0:
            y_shearcenter = 0
        else:
            y_shearcenter = moment_Qz_total/Q_z
    
        if Q_y == 0:
            z_shearcenter = 0
        else:
            z_shearcenter = moment_Qy_total/Q_y
    
        return np.array((y_shearcenter, z_shearcenter))

    def shearflow_torsion(self, torsionalmoment):
        # (constant) shearflow due to torsion
        # counterclockwise
        integrate_sGt_ = self.integrate_sGt()
        a = integrate_sGt_left = integrate_sGt_[0]
        b = integrate_sGt_center = integrate_sGt_[1]
        c = integrate_sGt_right = integrate_sGt_[2]
        
        polygon_left = np.concatenate((self.geometries[0], np.flip(self.geometries[1],0)))
        polygon_right = np.concatenate((self.geometries[2], np.flip(self.geometries[1],0)))
        Am_left = polygonarea(polygon_left)
        Am_right = polygonarea(polygon_right)
        
        u = a/Am_left + b/Am_left + b/Am_right
        v = c/Am_right + b/Am_right + b/Am_left
        w = 1 + (Am_right/Am_left)*(u/v)

        n_xS_left = torsionalmoment/(2*Am_left)*(1/w)
        n_xS_right = n_xS_left*(u/v)
        
        return n_xS_left, n_xS_right
        
    def deformation_vec(self, forces):
        # forces = ((N, M_y, M_z))
        EA_total = 0
        EI_total = 0
        for geometry, E, t in zip(self.geometries, self.youngsmoduli, self.thickness):
            EA_total += EA(geometry, E, t)
            EI_total += EI(geometry, E, t)

        EI_y = EI_total[0]
        EI_z = EI_total[1]
        EI_yz = EI_total[2]
        
        matrix = np.array([[EA_total, 0, 0],
                           [0, EI_y, EI_yz],
                           [0, EI_yz, EI_z]])
        
        deformations = np.linalg.solve(matrix, forces) # ϵ, κ_x, κ_y
        
        return deformations
    
    def deformation_element(self, forces):
        # deformation in every element due to normalforce and bending moments
        # deformation towards longitudinal axis of beam
        deformation_vec_ = self.deformation_vec(forces)
        ϵ = deformation_vec_[0]
        κ_y = deformation_vec_[1]
        κ_z = deformation_vec_[2]
        deformations_total = []
    
        for geometry in self.geometries:
            deformations = []
            
            for i in range(0,len(geometry)-1):
                a = np.array(geometry[i])
                b = np.array(geometry[i+1])
                distance = (a+b)/2
                Δy = distance[0]
                Δz = distance[1]
                deformations_i = ϵ + Δy*κ_y + Δz*κ_z
                deformations.append(deformations_i)
                
            deformations_total.append(deformations)
    
        return deformations_total
    
    def stress_element(self, forces):
        # stress in every element due to normalforce and bending moments
        # stress towards longitudinal axis of beam
        deformation_element_ = self.deformation_element(forces)
        stress_total = []
        
        for deformations, E in zip(deformation_element_, self.youngsmoduli):
            stress = []
            
            for i in range(len(deformations)):
                stress_i = deformations[i]*E[i]
                stress.append(stress_i)
                
            stress_total.append(stress)
    
        return stress_total
    
    def c4_point(self):
        # c4_point of profile (on chord line)
        position1 = np.argmin(self.geometries[0][:,0])
        firstpoint = self.geometries[0][position1]
        position2 = np.argmax(self.geometries[2][:,0])
        lastpoint = self.geometries[2][position2]
        c4_point = (lastpoint - firstpoint)/4 + firstpoint
        
        return c4_point
    
    def dist_c4_shearcenter(self):
        # distance between c4_point and shearcenter in y- and z-axis
        distance = self.shearcenter() - self.c4_point()
        return distance
    
    def transform_shearforce(self, shearforce):
        # transform shearforce in c4 into shearforce in shearcenter and torsional moment
        shearforce = transform(np.array((shearforce)), self.Θ)
        distance = self.dist_c4_shearcenter()
        T = np.cross(distance, shearforce)
        
        return shearforce, T
    
    def shearforce_and_torsionalmoment(self, shearforce, torsionalmoment):
        # returns shearforce due to shearcenter and total torsionalmoment
        # input is shearforce due to c4 point and torsionalmoment
        Q, T_trans = self.transform_shearforce(shearforce)
        T = T_trans + torsionalmoment
        
        return Q, T
    
    def shearflow_total(self, shearforce, torsionalmoment):
        # returns shearflow due to shearforce and torsionalmoment superimposed
        # same direction as shearflow_shearforce
        # input is shearforce due to c4 point and torsionalmoment
        
        Q, M_T = self.shearforce_and_torsionalmoment(shearforce, torsionalmoment)
        
        n_xSQz, n_xSQy = self.shearflow_shearforce(Q) # left, center, right
        n_xSQ = np.array(n_xSQz) + np.array(n_xSQy)
        n_xS_tor_left, n_xS_tor_right = self.shearflow_torsion(M_T)
        
        n_xS_left = n_xSQ[0] - n_xS_tor_left
        n_xS_center = n_xSQ[1] + n_xS_tor_left - n_xS_tor_right
        n_xS_right = n_xSQ[2] + n_xS_tor_right
        
        return n_xS_left, n_xS_center, n_xS_right
    
    def twist(self, shearforce, torsionalmoment):
        # twist is just depending on torsionalmoment; shearforce is due to shearcenter
        Q, M_T = self.shearforce_and_torsionalmoment(shearforce, torsionalmoment)
        
        integrate_sGt_ = self.integrate_sGt()
        a = integrate_sGt_left = integrate_sGt_[0]
        b = integrate_sGt_center = integrate_sGt_[1]
        c = integrate_sGt_right = integrate_sGt_[2]
        
        polygon_left = np.concatenate((self.geometries[0], np.flip(self.geometries[1],0)))
        polygon_right = np.concatenate((self.geometries[2], np.flip(self.geometries[1],0)))
        Am_left = polygonarea(polygon_left)
        Am_right = polygonarea(polygon_right)
        
        n_xS_left, n_xS_right = self.shearflow_torsion(M_T)
        
        twist = (1/(2*Am_left))*(n_xS_left*a + (n_xS_left - n_xS_right)*b)
        
        return twist
                
    def test_continuity(self, shearflow):
        # for test: integrals in both cells must be qual to zero
        # just for pure shear force, no torsional moment
        integral_nsGt = []
        for geometry, G, t, n in \
        zip(self.geometries, self.shearmoduli, self.thickness, shearflow):
            integral_nsGt_j = 0
            for i in range(0,len(geometry)-1):
                a = np.array(geometry[i])
                b = np.array(geometry[i+1])
                Δs_i = np.linalg.norm(a-b)
                integral_nsGt_j += n[i]*Δs_i/(G[i]*t[i])
            integral_nsGt.append(integral_nsGt_j)
            
        continuity_leftcell = integral_nsGt[0] - integral_nsGt[1]
        continuity_rightcell = integral_nsGt[2] - integral_nsGt[1]
      
        return continuity_leftcell, continuity_rightcell
    
    def plot_around_profil(self, shearflow):
    
        geometry_left = self.geometries[0]
        geometry_center = self.geometries[1]
        geometry_right = self.geometries[2]
        
        shearflow_left = shearflow[0]
        shearflow_center = shearflow[1]
        shearflow_right = shearflow[2]

        point_profil_left = []
        point_profil_center = []
        point_profil_right = []

        for i in range(0,len(geometry_left)-1):
            a = np.array(geometry_left[i])
            b = np.array(geometry_left[i+1])

            m = (a+b)/2
            n = b-a

            normalvec = [n[1], -n[0]]/np.linalg.norm(n)

            point_profil_left_i = normalvec*shearflow_left[i] + m
            point_profil_left.append(point_profil_left_i)

        for i in range(0,len(geometry_center)-1):
            a = np.array(geometry_center[i])
            b = np.array(geometry_center[i+1])

            m = (a+b)/2
            n = b-a

            normalvec = [n[1], -n[0]]/np.linalg.norm(n)

            point_profil_center_i = normalvec*shearflow_center[i] + m
            point_profil_center.append(point_profil_center_i)
            
        for i in range(0,len(geometry_right)-1):
            a = np.array(geometry_right[i])
            b = np.array(geometry_right[i+1])

            m = (a+b)/2
            n = b-a

            normalvec = [n[1], -n[0]]/np.linalg.norm(n)

            point_profil_right_i = normalvec*shearflow_right[i] + m
            point_profil_right.append(point_profil_right_i)

        x = np.array(point_profil_left)
        y = np.array(point_profil_center)
        z = np.array(point_profil_right)

        return x, y, z

    
def transform(vec, Θ):
    rotmatrix = np.array([[np.cos(Θ), -np.sin(Θ)],
                          [np.sin(Θ),  np.cos(Θ)]])

    vec_transformed = (rotmatrix @ vec.T).T

    return vec_transformed

        
def ES(geometry, youngsmoduli, thickness):
    ES_y = 0
    ES_z = 0
    for i in range(0,len(geometry)-1):
        a = np.array(geometry[i])
        b = np.array(geometry[i+1])
        Δs = np.linalg.norm(a-b)
        distance = (a+b)/2
        Δy = distance[0]
        Δz = distance[1]
        ES_y += youngsmoduli[i]*Δz*Δs*thickness[i]
        ES_z += youngsmoduli[i]*Δy*Δs*thickness[i]
    return np.array((ES_y, ES_z))


def ES_circulate(geometry, youngsmoduli, thickness):
    ES_y = 0
    ES_z = 0
    ES_y_vec = []
    ES_z_vec = []
    for i in range(0,len(geometry)-1):
        a = np.array(geometry[i])
        b = np.array(geometry[i+1])
        Δs = np.linalg.norm(a-b)
        distance = (a+b)/2
        Δy = distance[0]
        Δz = distance[1]
        ES_y += youngsmoduli[i]*Δz*Δs*thickness[i]
        ES_z += youngsmoduli[i]*Δy*Δs*thickness[i]
        ES_y_vec.append(ES_y)
        ES_z_vec.append(ES_z)
    return np.array((ES_y_vec, ES_z_vec))


def EA(geometry, youngsmoduli, thickness):
    EA = 0
    for i in range(0,len(geometry)-1):
        a = np.array(geometry[i])
        b = np.array(geometry[i+1])
        Δs = np.linalg.norm(a-b)
        EA += youngsmoduli[i]*Δs*thickness[i]
    return EA


def EI(geometry, youngsmoduli, thickness):
    EI_y = 0
    EI_z = 0
    EI_yz = 0
    for i in range(0,len(geometry)-1):
        a = np.array(geometry[i])
        b = np.array(geometry[i+1])
        Δs = np.linalg.norm(a-b)
        distance = (a+b)/2
        Δy = distance[0]
        Δz = distance[1]
        EI_y += youngsmoduli[i]*Δz**2*Δs*thickness[i]
        EI_z += youngsmoduli[i]*Δy**2*Δs*thickness[i]
        EI_yz += -youngsmoduli[i]*Δy*Δz*Δs*thickness[i]
    return np.array((EI_y, EI_z, EI_yz))


def polygonarea(polygon):
    A = 0
    for i in range(0,len(polygon)-1):
        a = polygon[i]
        b = polygon[i+1]    
        A += 1/2 * (b[1]*a[0] - a[1]*b[0])
    return abs(A)


def lift(ρ, v, c_l, c):
    L = ρ/2*v**2*c_l*c # lift per length
    return L


def moment(ρ, v, c_m, c):
    m = ρ/2*v**2*c_m*c**2 # torsion moment per length
    return m