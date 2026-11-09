function [K, grad_l] = Gram_Gibbs(X,Xp,sigma,l)
arguments
    X     (:,1) double
    Xp    (:,1) double
    sigma (1,1) function_handle
    l     (1,1) double
end

sig = sigma(X);
sig_star = sigma(Xp);
assert(numel(sig)==length(X) && isvector(sig), "sigma0(T) did not return an N-length vector.  Where N = length(T) = " + length(X))
assert(numel(sig_star)==length(Xp) && isvector(sig_star), "sigma0(Ts) did not return an M-length vector.  Where M = length(Ts) = " + length(Xp))

% pairwise squared distances
r2 = pdist2(X, Xp, 'squaredeuclidean');

% base kernel
E = exp(-r2/(2*l^2));
K = (sig(:) * sig_star(:).') .* E;

if nargout > 1
    % dK/dl
    grad_l = K .* (r2 / l^3);
end
end