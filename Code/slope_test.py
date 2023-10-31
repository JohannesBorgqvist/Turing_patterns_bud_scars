#!/usr/bin/env python
# coding: utf-8

# In[6]:


import math
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import pandas as pd


# In[48]:


df = pd.read_excel('Fig_5_spectral_n_1_d_20.xlsx', sheet_name='$U_{1}^{0}$')
df['mean'] = df.mean(axis=1)


# In[49]:


fig, axes = plt.subplots(figsize=(14,6))
axes.plot(df['Unnamed: 0'],df['mean'],'.')


# In[27]:


df = pd.read_excel('Fig_5_spectral_n_1_d_20.xlsx', sheet_name='$U_{1}^{0}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for $U_{1}^{0}$: ",p_value)


# In[28]:


df = pd.read_excel('Fig_5_spectral_n_1_d_20.xlsx', sheet_name='$U_{1}^{1}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for $U_{1}^{1}$: ",p_value)


# In[29]:


df = pd.read_excel('Fig_6_spectral_n_2_d_18.xlsx', sheet_name='$U_{2}^{0}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for $U_{2}^{0}$: ",p_value)


# In[30]:


df = pd.read_excel('Fig_6_spectral_n_2_d_18.xlsx', sheet_name='$U_{2}^{1}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for $U_{2}^{1}$: ",p_value)


# In[31]:


df = pd.read_excel('Fig_6_spectral_n_2_d_18.xlsx', sheet_name='$U_{2}^{2}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for $U_{2}^{2}$: ",p_value)


# In[34]:


df = pd.read_excel('Fig_7_metrics_n_1_d_20.xlsx', sheet_name='Minimal distance')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for 'Minimal distance (n,d)=(1,20)': ",p_value)


# In[35]:


df = pd.read_excel('Fig_8_metrics_n_2_d_18.xlsx', sheet_name='Minimal distance')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for 'Minimal distance (n,d)=(2,18)': ",p_value)


# In[36]:


df = pd.read_excel('Fig_S7_spectral_n_3_d_18.xlsx', sheet_name='$U_{3}^{0}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{3}^{0}$': ",p_value)


# In[37]:


df = pd.read_excel('Fig_S7_spectral_n_3_d_18.xlsx', sheet_name='$U_{3}^{1}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{3}^{1}$': ",p_value)


# In[38]:


df = pd.read_excel('Fig_S7_spectral_n_3_d_18.xlsx', sheet_name='$U_{3}^{2}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{3}^{2}$': ",p_value)


# In[39]:


df = pd.read_excel('Fig_S7_spectral_n_3_d_18.xlsx', sheet_name='$U_{3}^{3}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{3}^{3}$': ",p_value)


# In[40]:


df = pd.read_excel('Fig_S10_spectral_n_4_d_18.xlsx', sheet_name='$U_{4}^{0}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{4}^{0}$': ",p_value)


# In[41]:


df = pd.read_excel('Fig_S10_spectral_n_4_d_18.xlsx', sheet_name='$U_{4}^{1}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{4}^{1}$': ",p_value)


# In[42]:


df = pd.read_excel('Fig_S10_spectral_n_4_d_18.xlsx', sheet_name='$U_{4}^{2}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{4}^{2}$': ",p_value)


# In[43]:


df = pd.read_excel('Fig_S10_spectral_n_4_d_18.xlsx', sheet_name='$U_{4}^{3}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{4}^{3}$': ",p_value)


# In[44]:


df = pd.read_excel('Fig_S10_spectral_n_4_d_18.xlsx', sheet_name='$U_{4}^{4}$')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for '$U_{4}^{4}$': ",p_value)


# In[46]:


df = pd.read_excel('Fig_S8_metrics_n_3_d_18.xlsx', sheet_name='Minimal distance')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for 'Minimal distance (n,d)=(3,18)': ",p_value)


# In[47]:


df = pd.read_excel('Fig_S11_metrics_n_4_d_18.xlsx', sheet_name='Minimal distance')
df['mean'] = df.mean(axis=1)
x = df['Unnamed: 0']
y = df['mean']
slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
print("p-value for 'Minimal distance (n,d)=(4,18)': ",p_value)


# In[ ]:




