
# coding: utf-8

# In[1]:


import numpy as np
from wingstructure.structure2 import SectionBase, Layer, BoxSpar
from wingstructure.material import IsotropicMaterial as Material


# # Werkstoffe definieren

# In[2]:


balsa_schwer = Material(ρ=0.3e3, E=210e3, G=80e3)
balsa_mittel = Material(ρ=0.18e3, E=210e3, G=80e3)
balsa_leicht = Material(ρ=0.16e3, E=210e3, G=80e3)
gfk_rovings  = Material(ρ=1.225, E=210e3, G=80e3)
gfk_gewebe   = Material(ρ=1.225, E=210e3, G=80e3)


# # Geometrie definieren

# ## Profilkoordinaten einlesen

# In[3]:


coords = np.loadtxt('airfoils/ah93157.dat', skiprows=1)*1.2
sektion = SectionBase(coords)
sektion


# In[4]:


außenlage = Layer(sektion, gfk_gewebe, 3e-2)
kastenholm = BoxSpar(außenlage, {'flange': gfk_rovings, 'web': balsa_mittel},
                                0.45, 0.200, 3e-2, 1e-2)
l3 = Layer(kastenholm, balsa_schwer, 3e-3)
sektion


# In[5]:


import triangle

außenlage.geometry.exterior


# In[6]:


t1 = triangle.triangulate({'vertices':np.array(außenlage.geometry.exterior)})


# In[7]:


alle = np.vstack((np.array(außenlage.geometry.exterior), np.array(außenlage.geometry.interiors[0])))
#t2 = triangle.triangulate({'vertices':alle}, 'l')


# In[8]:


#t2.keys()


# In[ ]:


#a = np.array([1,1])
segments = [(ii-1,ii) for ii in range(np.array(außenlage.geometry.exterior).shape[0])]
segments[0] = (119,0)
segments = np.array(segments)
segments2 = [(120+ii-1,120+ii) for ii in range(np.array(außenlage.geometry.interiors[0]).shape[0])]
segments2[0] = (120+95,120)
segments2 = np.array(segments2)
#a =np.array([1,2]) 
t2 = triangle.triangulate({'vertices':alle, 'segments':np.vstack((segments,segments2,))},'p')

#holes = np.array((1,1))
# In[1]:


segments, segments2


# In[12]:



import triangle.plot as plot
from matplotlib import pyplot as plt
fig = plt.figure(figsize=(25,15))
ax = plt.subplot(111)
plot.plot(ax, **t2)
plt.show()
