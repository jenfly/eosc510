
# coding: utf-8

# # EOSC 510 Tutorial 2
# 
# ## Testing multiple linear regression on artificial data

# In[215]:

get_ipython().magic(u'matplotlib inline')

from IPython.display import display 
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.stats
import statsmodels.formula.api as smf

mpl.rcParams['legend.fontsize'] = 'small'


# ### Quick reference of relevant Python commands:
# 
# Operation | Matlab | Python
# --- | --- | ---
# Matrix multiplication | `X * Y` | `np.dot(X, Y)` 
# Matrix transpose | `X'` | `X.transpose()`
# Matrix inverse | `X^(-1)` | `np.linalg.inv(X)`
# Correlation coeffs & p-values | `corrcoef(x1, x2)` | `scipy.stats.linregress(x1, x2)`
# Multiple linear regression | `regress()` | `statsmodels.formula.api.ols()`
# 
# 

# ### Example 1
# 
# Create artificial data set based on:
# $Y = a_0 + a_1X_1 + a_2X_2 + a_3X_3 + a_4X_4$

# In[216]:

def readmat(filename, varname):
    """Read selected variable from .mat and output as array of floats"""
    data_in = sio.loadmat(filename)
    return data_in[varname].astype(float)


# In[217]:

X = readmat('Xdata.mat', 'X')
n, npred = X.shape
ind = np.arange(1, n + 1)
data = pd.DataFrame(X, index=ind, columns=['X1', 'X2', 'X3', 'X4'])
data.head()


# In[218]:

a0, a1, a2, a3, a4 = 0, 1, -2, 3, -4
Y = a0 + a1*data['X1'] + a2*data['X2'] + a3*data['X3'] + a4*data['X4']
data['Y'] = Y
data.head()


# In[219]:

style = {'X1' : 'b', 'X2' : 'r', 'X3' : 'g', 'X4' : 'm', 'Y' : 'k'}
data.plot.line(style=style)


# In[220]:

# Multiple linear regression model
def mlr(formula='Y ~ X', data=None):
    lm = smf.ols(formula=formula, data=data).fit()
    Y_regr = lm.predict()
    display(lm.summary())
    return lm, Y_regr


# In[221]:

lm, Y_regr = mlr(formula='Y ~ X1 + X2 + X3 + X4', data=data)
data['Y_regr'] = Y_regr


# In[222]:

style['Y_regr'] = 'y--'
data.plot.line(style=style)


# Instead of using `smf.ols()` lets derive the coefficients from  minimizing sum of squared errors (SSE) with respect to a 
# $a=(X^T X)^{-1} X^T y$

# In[223]:

def calc_coeffs(X, Y):
    return np.dot(np.linalg.inv(np.dot(X.transpose(), X)), np.dot(X.transpose(), Y))
areg = calc_coeffs(X, Y)
areg


# Can't find a direct equivalent of Matlab's stepwisefit() function in any of the Python stats modules.  Could write a function myself to mimic the Matlab code:
# ```
# % applying stepwise regression on Y and X
# [a SE PVAL INMODEL STATS NEXTSTEP HISTORY]=stepwisefit(X,Y);
# % constant coefficient 
# a0_st=STATS.intercept
# 
# ```

# In[224]:

# Simple linear regression between Y and Y_regr

# lm2 = smf.ols(formula='Y_regr ~ Y', data=data).fit()
# lm2.summary()

def corrcoef(x1, x2):
    """Return correlation coeff and p-value"""
    m, b0, r, p, stderr = scipy.stats.linregress(x1, x2)
    print('Correlation coefficient: %f' % r)
    print('p-value: %f' % p)
    return r, p

def scatterplot(data, xname, yname):
    """Scatter plot and one-to-one line"""
    x, y = data[xname], data[yname]
    r, p = corrcoef(x, y)
    data.plot.scatter(x=xname, y=yname)
    plt.title('r=%.2f' % r)
    combo = np.concatenate((x, y))
    xline = np.arange(min(combo), max(combo) + 1)
    plt.plot(xline, xline, 'k')


# In[225]:

scatterplot(data, 'Y', 'Y_regr')


# ### Example 2
# 
# Create artificial data set based on:
# $Y = a_0 + a_1X_1 + a_2X_2 + a_3X_3 + a_4X_4 + Y_{rand}$

# In[226]:

Yrand = readmat('Yrand.mat', 'Yrand')
Yrand = Yrand.flatten()
Ynew = Y + 5 * Yrand
data = data.drop('Y_regr', axis=1)
data['Ynew'] = Ynew
data.head()


# In[227]:

style['Ynew'] = 'k--'
data.plot.line(style=style)


# In[228]:

# Multiple linear regression model
lm, Y_regr = mlr(formula='Ynew ~ X1 + X2 + X3 + X4', data=data)
data['Y_regr'] = Y_regr


# In[230]:

data.plot.line(style=style)


# In[231]:

scatterplot(data, 'Ynew', 'Y_regr')


# In[ ]:




# In[ ]:




# In[ ]:



