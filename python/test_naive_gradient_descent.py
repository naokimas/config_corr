# -*- coding: utf-8 -*-
# test max_ent_config_naive_gradient_descent.py
import numpy as np
import numpy.matlib

data_type = 4

curr_dir = '../../../data/'
# curr_dir = '.' # should be modified to be consistent with your folder name

# We provide only motivation.txt (as described in Masuda, Kojaku & Sano, Physical Review E, 2018)
# Matrix C can be either a covariance or correlation matrix
if data_type==1: # academic motivation data
    C = np.loadtxt(curr_dir + 'motivation.txt')
elif data_type==2: # fMRI1
    C = np.loadtxt(curr_dir + 'cov105115_all.txt');
elif data_type==3: # fMRI2
    C = np.loadtxt(curr_dir + 'cov118932_all.txt');
elif data_type==4: # Japan
    C = np.loadtxt(curr_dir + 'JapanCov.txt');
elif data_type==5: # US
    C = np.loadtxt(curr_dir + 'USCov.txt');

N = C.shape[0] # size of the matrix

transform_to_corr_mat = True
# if True, transform the input covariance matrix to the correlation matrix before running the gradient descent method
# default value is True

# eigenvalues of the original covariance or correlation matrix
eigs_org, eig_vector = np.linalg.eig(C)
min_eig_C = min(eigs_org)
if min_eig_C < 1e-6 * max(np.diag(C)):
    print('min eig = %f' % min_eig_C)
    print('Input correlation/covariance matrix must be full rank')

tolerance = 1e-5; # to judge whether the algorithm has converged. In the paper = 1e-5

r = 1e-4; # learning rate. If r is too large, the algorithm would not converge
# default r = 1e-4

C_con, alpha, beta, it = max_ent_naive_gradient_descent(C, tolerance, r, transform_to_corr_mat);

if transform_to_corr_mat == True:
    D_n_05 = np.diag(np.power(np.diag(C_con),-0.5))
    C_con = np.dot(np.dot(D_n_05, C_con), D_n_05) # corr matrix
K_con = np.linalg.inv(C_con); # precision matrix

L = N; # L should be set to the length of the original data for which the Pearson correlation is calculated (e.g., length of the time series)

X = np.random.multivariate_normal(np.zeros(N), C_con, L);
C_sam = X.T * X / L; # sample covariance matrix generated by the estimated configuration model
if transform_to_corr_mat == True:
    D_n_05 = np.diag(np.power(np.diag(C_sam),-0.5))
    C_sam = np.dot(np.dot(D_n_05, C_sam), D_n_05) # sample correlation matrix
K_sam = np.linalg.inv(C_sam); # sample precision matrix

print(K_con)