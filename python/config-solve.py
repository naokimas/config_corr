import numpy as np
import numpy.matlib

def max_ent_K_config(C, tolerance, r, transform_to_corr_mat):
#
# Input
#   C: covariance matrix
#   tolerance: tolerance in relative error
#   r: learning rate
#   transform_to_corr_mat = True if one transforms the input covariance matrix to the correlation matrix before running the gradient descent method. Otherwise False.
#   default value is True
#
# Output
#   C_con: estimated covariance matrix
#   alpha, beta: parameters of the estimated covariance/precision matrix
#   it: iteration number
#
# estimated precision matrix K_est is given by
# K_est(i,i) = alpha(i)
# K_est(i,j) = beta(i) + beta(j)

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
#        if corr_preserve == True:
#            D_n_05 = np.diag(np.power(np.diag(C_con),-0.5))
#            corr_est = np.dot(np.dot(D_n_05, C_con), D_n_05) # correlation matrix
#            beta = beta + r * (1/N) * (np.sum(corr_est, axis=1) - s) # gradient descent on the correlation matrix
#        else:
        beta = beta + r * (1/N) * (np.sum(C_con, axis=1) - s) # gradient descent on the covariance matrix

        if it % 1000 == 0: # Then measure the relative error
#            if corr_preserve == True:
#                error = (sum(abs((np.diag(C) - np.diag(C_con)) / np.diag(C))) + sum(abs(((s - np.sum(corr_est, axis=1)) / s)))) / N;
#            else:
            error = (sum(abs((np.diag(C) - np.diag(C_con)) / np.diag(C))) + sum(abs(((s - np.sum(C_con, axis=1)) / s)))) / N;
            print('%f %f %f' % (error, np.mean(abs(alpha)), np.mean(abs(beta))))
        it = it + 1

    return C_con, alpha, beta, it
# end of def max_ent_K_config

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

C_con, alpha, beta, it = max_ent_K_config(C, tolerance, r, transform_to_corr_mat);

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

# K = np.linalg.inv(C); # empirical precision matrx
print(K_con)