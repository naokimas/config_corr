# -*- coding: utf-8 -*-
import cvxpy as cv
import numpy as np

def max_ent_config_dmcc(C, tolerance = 1e-5, transform_to_corr_mat = True):
    """
     
    DET-MAX algorithm for estimating the configuration model for correlation matrix data 
    
    Input
      C: covariance matrix
      tolerance: tolerance in relative error
      transform_to_corr_mat = True if one transforms the input covariance matrix to the correlation matrix before running the gradient descent method. Otherwise False.
      The default value is True.
    
    Output
      C_con: estimated covariance matrix
    
    """

    if transform_to_corr_mat == True: # Work on the correlation matrix
        # Transform the original covariance matrix to the correlation matrix
        cov = np.asanyarray(C)
        std_ = np.sqrt(np.diag(cov))
        _C = cov / np.outer(std_, std_)
    else: #  work on the covariance matrix
        _C = C

    N = _C.shape[0] # Number of nodes

    # Covariance matrix we will estimate 
    C_con = cv.Variable((N, N), PSD=True) 

    # Objective function to be maximized
    objective = cv.Minimize(-cv.log_det(C_con))  

    # Constraints on C_con 
    constraints = [cv.sum(C_con, axis=0) == _C.sum(axis = 0), cv.diag(C_con) == np.diag(_C)]

    # Optimization 
    prob = cv.Problem(objective, constraints)
    prob.solve(verbose = True, eps = tolerance)

    return C_con.value
