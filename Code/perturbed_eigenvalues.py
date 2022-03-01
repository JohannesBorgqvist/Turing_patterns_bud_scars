#!/usr/bin/env python
# coding: utf-8

# In[1]:


# Importera nödvändiga moduler
#get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams.update({'font.size': 21})


# In[4]:


#Parametervärden
a = 0.2
b = 1
d = 18
k_squared=2*(2+1)


# In[5]:


#Partiella deviator etc.
f_u = ( (b - a) / (b + a) )
f_v = ( (a + b)**2 )
g_u = ( (-2*b) / (a+b) )
g_v = -f_v
d_c = ( (f_u*g_v) - (2*f_v*g_u) )
d_c +=  np.sqrt( d_c**2 - ((f_u**2)*(g_v**2)) )
d_c = ( (d_c) / (f_u**2) )
gamma_c = ( ( 2 * d_c * k_squared) / ( (d_c*f_u) + g_v ) )


# In[8]:


#Beräkna gränserna
lm = d*(b-a) - (a+b)**3
gL = gamma_c*(lm - np.sqrt(lm**2 -4*d*(a+b)**4))/(2*d*(a+b))
gM = gamma_c*(lm + np.sqrt(lm**2 -4*d*(a+b)**4))/(2*d*(a+b))
print(gL,gM)


# In[75]:


#Plotta perturberade egenvärden
eps_max=0.5 #Största värdet på epsilon
epsilon = np.linspace(0,eps_max)
ev=[2,3,4] #eigenvalues to consider
fig, axes = plt.subplots(1,1,figsize=(12,8))
for n in ev:
    evp = n*(n+1)+4*(2*n+1)/(n**2+n)*epsilon**2 #m=0 egenvärdet
    axes.plot(epsilon,evp)
    axes.text(epsilon[-1], evp[-1],"$\lambda_{" + str(n)+",0}$", fontsize=12)
    for m in range(1,n+1): #loopa igenom övriga egenvärden
        Cnm_nom = (2*n+1)*np.math.factorial(m+n)
        Cnm_den = (4**m)*np.math.factorial(n-m)*np.math.factorial(m)*np.math.factorial(m-1)
        evp = n*(n+1)-Cnm_nom/Cnm_den*epsilon**(2*m) #m=0 egenvärdet
        axes.plot(epsilon,evp)
        axes.text(epsilon[-1], evp[-1],"$\lambda_{" + str(n)+","+ str(m)+"}$", fontsize=14)
axes.plot([0,eps_max],gL*np.array([1,1]),'k--')
axes.plot([0,eps_max],gM*np.array([1,1]),'k--')
axes.set_xlabel("$\epsilon$")
plt.show()


# In[ ]:
