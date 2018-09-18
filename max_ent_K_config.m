
function [C_con, alpha, beta, it] = max_ent_K_config(C, tolerance, r, corr_preserve)
%
% Input
%   C: covariance matrix
%   tolerance: tolerance in relative error
%   r: learning rate
%   corr_preserve = 1 if the node strength in terms of the Pearson corr 
%     should be preserved. = 0 if the node strength in terms of the
%     covariance should be preserved. Should be set to 0, and one should feed
%     a correlation matrix as C.
%
% Output
%   C_con: estimated covariance matrix
%   alpha, beta: parameters of the estimated covariance/precision matrix
%   it: iteration number
%
% estimated precision matrix K_est is given by
% K_est(i,i) = alpha(i)
% K_est(i,j) = beta(i) + beta(j)

N = size(C,1); % number of nodes

if corr_preserve==1 % work on the correlation matrix
    corr = diag(diag(C))^(-1/2) * C * diag(diag(C))^(-1/2); % corr matrix
    s = sum(corr,2); % node strength including the self-loop
else % work on the covariance matrix
    s = sum(C,2); % node strength including the self-loop
end

K = inv(C); % precision matrix from data
% sum(sum(abs(K*C-eye(N))))
alpha = diag(K); % initialization
beta = zeros(N,1); % initialization

% The following initialization is not used.
% beta = mean(mean(K)) / 2 * ones(N,1); % initialization

% The following initialization is not used.
%
% U = ones(N,N);
% for i=1:N
%    U(i,i) = N-1;
% end
% beta = inv(U) * (sum(K,2) - diag(K)); % initialization
% inv U may be analytically derived
    
error = 1e+3; % initialization
it=0; % number of iteration
while (error > tolerance) 

    % construct the precision matrix from the current estimate
    K_est = repmat(beta,1,N) + repmat(beta,1,N)';
    for i=1:N
        K_est(i,i) = K_est(i,i) + alpha(i);
    end
    
    % covariance matrix from the current estimate
    C_con = inv(K_est);

    % gradient descent on the log likelihood
    alpha = alpha + r * (diag(C_con) - diag(C));
    if corr_preserve==1
        corr_est = diag(diag(C_con))^(-1/2) * C_con * diag(diag(C_con))^(-1/2); % correlation matrix
%        for i=1:N
%            corr_est(i,i) = 1.0;
%            for j=1:i-1
%                corr_est(i,j) = C_con(i,j) / sqrt(C_con(i,i)*C_con(j,j));
%                corr_est(j,i) = corr_est(i,j);
%            end
%        end
        beta = beta + r * (1/N) * (sum(corr_est,2) - s); % gradient descent on the correlation matrix
    else
        beta = beta + r * (1/N) * (sum(C_con,2) - s); % gradient descent on the covariance matrix
    end

    if (mod(it,1000) == 0) % Then measure the relative error
        if corr_preserve==1
            error = (sum(abs((diag(C) - diag(C_con)) ./ diag(C))) + ...
                sum(abs(((s - sum(corr_est,2)) ./ s)))) / N;
        else
            error = (sum(abs((diag(C) - diag(C_con)) ./ diag(C))) + ...
                sum(abs(((s - sum(C_con,2)) ./ s)))) / N;
        end
        fprintf('%f %f %f\n', error, mean(abs(alpha)), mean(abs(beta)));
    end
    it = it + 1;
end % end of while loop
end
