% This code was taken from https://github.com/Leron33/Graph-matching
%it corresponds to the implementation for the paper "Graph Matching with Partially-Correct Seeds" by Liren Yu, Jiaming Xu and Xiaojun Lin. See the paper at https://arxiv.org/pdf/2004.03816.pdf.
% comment inn << >> are theirs... 
%<<Improved implementation of our new robust spectral method 
% A and B are the matrices to be matched 
% eta is the tuning parameter comparable to noise level sigma
% Return permutation matrix P so that P*A*P' is matched to B >>

function [P] = matching_robust_spectral(A, B, eta)
    n = size(A, 1);
    [U, Lambda] = eig(A);
    [V, Mu] = eig(B);
    lambda = diag(Lambda);
    mu = diag(Mu);
    Coeff = 1 ./ ((lambda - mu').^2 + eta^2);
    Coeff = Coeff .* (U' * ones(n) * V);
    X = U * Coeff * V';
    
    %% Rounding by linear assignment - better but slower 
    % M = matchpairs(X', -99999, 'max'); % HT: does not work, changed it
    % below
    %[Mr, Mc] = linear_sum_assignment(-X');
    %P = full(sparse(Mr, Mc, 1, n, n));

    %% Greedy matching - faster but worse 
%     P = full(greedy_match(X'));

    %% Greedy rounding
%     [~, ind_max] = max(X);
%     P = full(sparse(1:n, ind_max, 1, n, n));
%%  Rounding with GMWM 
P=GMWM_alg(X',-2000);

