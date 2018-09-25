# -*- coding: utf-8 -*-
import numpy as np
import numpy.matlib

def max_ent_config_naive_gradient_descent(C, tolerance, r, transform_to_corr_mat):
    """
     
    A naive gradient descent algorithm for estimating the configuration model for correlation matrix data, as presented in Naoki Masuda, Sadamori Kojaku, Yukie Sano. Physical Review E, 98, 012312 (2018).
    
    Input
      C: covariance matrix
      tolerance: tolerance in relative error
      r: learning rate
      transform_to_corr_mat = True if one transforms the input covariance matrix to the correlation matrix before running the gradient descent method. Otherwise False.
      The default value is True.
    
    Output
      C_con: estimated covariance matrix
      alpha, beta: parameters of the estimated covariance/precision matrix
      it: iteration number

      The estimated precision matrix K_est is given by
      K_est(i,i) = alpha(i)
      K_est(i,j) = beta(i) + beta(j)
    
    """

    N = C.shape[0] # number of nodes

    if transform_to_corr_mat == True: # work on the correlation matrix
        # transform the original covariance matrix to the correlation matrix
        D_n_05 = np.diag(np.power(np.diag(C),-0.5)) # = D^{-1/2}
        C = np.dot(np.dot(D_n_05, C), D_n_05) # = D^{-1/2} * C * D^{-1/2}
    s = np.sum(C, axis=1) # node strength including the self-loop

    K = np.linalg.inv(C) # precision matrix from data
    alpha = np.diag(K) # initialization
    beta = np.zeros(N) # initialization
    error = 1e+3 # initialization
    it=0 # number of iteration

    while error > tolerance:
        # construct the precision matrix from the current estimate
        K_est = np.matlib.repmat(beta,N,1) + np.transpose(np.matlib.repmat(beta,N,1))
        K_est = K_est + np.diag(alpha)
    
        # covariance matrix from the current estimate
        C_con = np.linalg.inv(K_est)

        # gradient descent on the log likelihood
        alpha = alpha + r * (np.diag(C_con) - np.diag(C))
        beta = beta + r * (1/N) * (np.sum(C_con, axis=1) - s) # gradient descent on the covariance matrix

        if it % 1000 == 0: # Then measure the relative error
            error = (sum(abs((np.diag(C) - np.diag(C_con)) / np.diag(C))) + sum(abs(((s - np.sum(C_con, axis=1)) / s)))) / N;
            print('%f %f %f' % (error, np.mean(abs(alpha)), np.mean(abs(beta))))
        it = it + 1

    return C_con, alpha, beta, it