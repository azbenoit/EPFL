# -*- coding: utf-8 -*-
"""Exercise 3.

Ridge Regression
"""

import numpy as np

def mse(y, tx, w):
    """Calculate the loss using either MSE or MAE.

    Args:
        y: shape=(N, )
        tx: shape=(N,2)
        w: shape=(2,). The vector of model parameters.

    Returns:
        the value of the loss (a scalar), corresponding to the input parameters w.
    """
    #MSE loss is (1/2*N)(e')(e)
    #where e = y-Xw
    e = y - np.dot(tx, w)
    return (1/(2*y.shape[0])) * (e.T @ e)

def ridge_regression(y, tx, lambda_):
    """implement ridge regression.

    Args:
        y: numpy array of shape (N,), N is the number of samples.
        tx: numpy array of shape (N,D), D is the number of features.
        lambda_: scalar.

    Returns:
        w: optimal weights, numpy array of shape(D,), D is the number of features.

    >>> ridge_regression(np.array([0.1,0.2]), np.array([[2.3, 3.2], [1., 0.1]]), 0)
    array([ 0.21212121, -0.12121212])
    >>> ridge_regression(np.array([0.1,0.2]), np.array([[2.3, 3.2], [1., 0.1]]), 1)
    array([0.03947092, 0.00319628])
    """
    # ***************************************************
    # COPY YOUR CODE FROM EX03 HERE
    # ridge regression: TODO
    best_weights = np.linalg.solve(tx.T @ tx + (2*y.shape[0]*lambda_) * np.eye(tx.shape[1]), tx.T @ y)
    best_loss = mse(y, tx, best_weights)
    return best_weights
    # ***************************************************
#     raise NotImplementedError
