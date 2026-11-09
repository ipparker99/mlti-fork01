function samples = sample_posterior_predictive(T, Ts, mu0, sigma0, mu_post, S_post, N_samp)
arguments
    T       (:,1) double
    Ts      (:,1) double
    mu0     (1,1) function_handle
    sigma0  (:,1) cell
    mu_post (:,1) double
    S_post  (:,:) double
    N_samp  (1,1) int64
end
assert(all(cellfun(@(f) isa(f,'function_handle'), sigma0)));

N = length(T);
M = length(Ts);

MU0 = mu0(T);
MU0s = mu0(Ts);
assert(size(MU0,1)==N, "mu0(T) did not return an N-by-D matrix.  Where N = length(T) = " + N)
assert(size(MU0s,1)==M, "mu0(Ts) did not return an M-by-D matrix.  Where M = length(Ts) = " + M)
assert(size(MU0,2)==size(MU0s,2))
D = size(MU0, 2);
mu0v = MU0(:);
mu0v_star = MU0s(:);
assert(length(sigma0)==D)

assert(length(mu_post)==N*D+D);
assert(all(size(S_post) == [N*D+D, N*D+D]));

samps_post = mvnrnd(mu_post, S_post, N_samp);
mv_samps = samps_post(:,1:N*D);
theta_samps = samps_post(:, N*D+1:end);

samples = zeros(N_samp, M, D);
for s = 1:N_samp
    K = cell(D,1);
    Ks = cell(D,1);
    Kss = cell(D,1);
    for i = 1:D
        l = exp(theta_samps(s,i));
        K{i} = Gram_Gibbs(T,T,sigma0{i},l);
        Ks{i} = Gram_Gibbs(T,Ts,sigma0{i},l);
        Kss{i} = Gram_Gibbs(Ts,Ts,sigma0{i},l);
    end
    K = blkdiag(K{:});
    Ks = blkdiag(Ks{:});
    Kss = blkdiag(Kss{:});
    Vs = (K+1e-8*eye(size(K)))\Ks;

    mu_star = (mu0v_star + Vs.' * (mv_samps(s,:).' - mu0v)).';
    S_star = (Kss+1e-8*eye(size(Kss))) - Vs.'*Ks;
    S_star = 0.5*(S_star+S_star.');
    samp = mvnrnd(mu_star, S_star);
    samples(s,:,:) = reshape(samp, size(MU0s));
    disp("Sampling Progress: " + s + "/" + N_samp)
end
end