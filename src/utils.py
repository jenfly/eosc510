import itertools
from IPython.display import display
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio
import scipy.stats
from sklearn.decomposition import PCA
import statsmodels.formula.api as smf

def standardize(df, cols=None):
    """Return DataFrame with standardized columns.

    Specify a subset of columns with the `cols` input.  If omitted, all columns
    are standardized.

    Columns are standardized as: y_out = (y - y.mean()) / y.std()
    """
    if cols is None:
        cols = df.columns
    df_out = df.copy()
    for nm in cols:
        df_out[nm] = (df[nm] - df[nm].mean()) / df[nm].std()
    return df_out


def mlr(formula='Y ~ X', data=None, verbose=True):
    """Multiple linear regression model"""
    lm = smf.ols(formula=formula, data=data).fit()
    Y_regr = lm.predict()
    if verbose:
        display(lm.summary())
    return lm, Y_regr

def corrcoef(x1, x2, verbose=False):
    """Return correlation coeff and p-value"""
    m, b0, r, p, stderr = scipy.stats.linregress(x1, x2)
    if verbose:
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

def split_data(data, isplit=25, verbose=True):
    """Split data into calibration and validation sets"""
    isplit = 25
    data_cal = data.iloc[:isplit]
    data_val = data.iloc[isplit:]
    if verbose:
        print('Calibration Data:')
        display(data_cal)
        print('\nValidation Data:')
        display(data_val)
    return data_cal, data_val


def get_formula(pred_list, response='Ynew'):
    """Return regression formula for list of predictor indices"""
    return response + ' ~ ' + ' + '.join([pred for pred in pred_list])


def formula_list(v0=['X1', 'X2', 'X3', 'X4'], response='Ynew'):
    combos = []
    for i, v in enumerate(v0):
        combos = combos + list(itertools.combinations(v0, i + 1))
    formulas = [get_formula(combo, response) for combo in combos]
    return formulas


def cross_validate(data_cal, data_val, formula='', response='Y'):
    """Cross validate a single MLR model"""
    lm = smf.ols(formula=formula, data=data_cal).fit()
    Yval = lm.predict(data_val)
    r_val, p_val = corrcoef(data_val[response], Yval, verbose=False)
    return r_val, p_val


def optimal_mlr(data_cal, data_val, data=None, formulas=[], response='Y', verbose=True):
    """Find best MLR fit using cross validation
    """
    r_vals, p_vals = [], []
    for formula in formulas:
        r, p = cross_validate(data_cal, data_val, formula=formula, response=response)
        r_vals.append(r)
        p_vals.append(p)

    results = pd.DataFrame(index=formulas, columns=['r_val', 'p_val'])
    results['r_val'], results['p_val'] = r_vals, p_vals
    if verbose:
        display(results)

    formula_best = results.idxmin()['p_val']
    if verbose:
        print('Best formula ' + formula_best)
    lm_best, _ = mlr(formula=formula_best, data=data_cal, verbose=verbose)
    if data is not None:
        Yfinal = lm_best.predict(data)
        rfinal, pfinal = corrcoef(data[response], Yfinal, verbose=False)
    else:
        Yfinal, rfinal, pfinal = None, None, None
    output = {'results' : results, 'formula_best' : formula_best,
              'lm_best' : lm_best, 'Yfinal' : Yfinal, 'rfinal' : rfinal,
              'pfinal' : pfinal}
    return output


def princomp(y, kmax=None, eigenrows=True, real=True):
    """Perform principal component analysis on matrix y.

    Uses numpy functions to compute eigenvalues and eigenvectors.

    Parameters
    ----------
    y : np.array (2-dimensional)
        Input data matrix.  PCA is performed along the rows of y.
    kmax : int, optional
        Truncates at first kmax modes (or returns all modes if kmax is None).
    eigenrows : bool, optional
        If eigenrows is True (default), eigenvectors and principal components
        are output as rows of their respective matrices, and the data y can be
        reconstructed as:
          y_rec = np.dot(A, E)
          where A = output['scores'] is the matrix of principal components and
          E = output['eigenvec'] is the matrix of eigenvectors.
        If eigenrows is False, then y_rec = np.dot(A.T, E.T)
    real : bool, optional
        If True, return eigenvectors, eigenvalues, etc. as real, otherwise
        return complex arrays.

    Output
    ------
    pca : dict
        Dictionary of eigenvectors, eigenvalues, principal components (scores),
        fraction of variance for each mode, original data (y_orig), and
        reconstructed data (y_rec).
    """

    # Subtract the mean (along columns) and transpose
    mean_y = np.mean(y.T, axis=1)
    yn = (y - mean_y).T

    # Compute the covariance matrix
    s = np.cov(yn)

    # Compute eigenvalues and eigenvectors and sort descending
    eigval, eigvec = np.linalg.eig(s)
    if real:
        eigval, eigvec = eigval.real, eigvec.real
    idx = np.argsort(eigval)
    idx = idx[::-1]
    eigval = eigval[idx]
    eigvec = eigvec[:, idx]    

    # Transpose eigvec so each row is an eigenvector
    eigvec = eigvec.T

    # Fraction of variance explained by each mode
    variance = eigval / eigval.sum()

    # Truncate at kmax
    if kmax is not None:
        eigval = eigval[:kmax]
        eigvec = eigvec[:kmax]
        variance = variance[:kmax]

    # Compute principal component scores
    a = np.dot(eigvec, yn)
    a = a.T

    # Reconstruct y from the eigenvectors and principal components
    y_rec = np.dot(a, eigvec)

    if not eigenrows:
        eigvec, a = eigvec.T, a.T

    pca = {'eigval' : eigval, 'eigvec' : eigvec, 'varfrac' : variance,
           'scores' : a, 'y_rec' : y_rec, 'y_orig' : y}

    return pca


def pca_skl(y, kmax=None, **kwargs):
    """Wrapper for sklearn.decomposition.PCA()

    Returns output in a more user friendly format.
    **kwargs are optional inputs to sklearn.decomposition.PCA()
    """
    pca_in = PCA(n_components=kmax, **kwargs)
    pca_in.fit(y)
    eigvec = pca_in.components_
    variance = pca_in.explained_variance_ratio_
    eigval = pca_in.explained_variance_

    mean_y = np.mean(y.T, axis=1)
    yn = (y - mean_y).T
    a = np.dot(eigvec, yn)
    a = a.T

    y_rec = pca_in.inverse_transform(a)

    pca = {'eigval' : eigval, 'eigvec' : eigvec, 'varfrac' : variance,
           'scores' : a, 'y_rec' : y_rec, 'y_orig' : y}

    return pca
