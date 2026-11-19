clear; clc; close all;

true_params = load("true_params.mat");

% [kf(T1),...,kf(TN),Cf(T1),...,Cf(TN),ks(T1),...,ks(TN),Cs(T1),...,Cs(TN)]
indx = [1:14,16:29];
T = (300:100:900)';
N = length(T);
Nsamp = 1e4;

% Store the smallest HDI interval of the following posteriors that contains
% the true diffusivities
Dfs_smallest_HDI_joint = zeros(100,1); % p([Df(T1),...,Df(TN),Ds(T1),...,Ds(TN)]|data)
Afs = [eye(N), -eye(N), zeros(N,2*N);
       zeros(N,2*N), eye(N), -eye(N)];
Df_smallest_HDI_joint = zeros(100,1);  % p([Df(T1),...,Df(TN)]|data)
Af = [eye(N), -eye(N), zeros(N,2*N)];
Ds_smallest_HDI_joint = zeros(100,1);  % p([Ds(T1),...,Ds(TN)]|data)
As = [zeros(N,2*N), eye(N), -eye(N)];
Df_smallest_HDI_ptws = zeros(100,7);   % p(Df(Ti)|data)
Afi = num2cell(Af, 2);
Ds_smallest_HDI_ptws = zeros(100,7);   % p(Ds(Ti)|data)
Asi = num2cell(As, 2);

Df_true =  polyval(true_params.gold_film.k, T) ...
        ./ polyval(true_params.gold_film.c, T) ...
        ./ polyval(true_params.gold_film.rho, T);
Ds_true =  polyval(true_params.YSZ.k, T) ...
        ./ polyval(true_params.YSZ.c, T) ...
        ./ polyval(true_params.YSZ.rho, T);

for trialNum = 1:100
    data = load("results\isoisoex_MAP_results_" + sprintf('%03d', trialNum) + ".mat");
    Hp = schur_comp(data.hessian, indx);
    mup = data.x(indx);

    mu_fs = Afs*mup; Cov_fs = Afs/Hp*Afs.'; Cov_fs = 0.5*(Cov_fs+Cov_fs.');
    mu_f = Af*mup; Cov_f = Af/Hp*Af.'; Cov_f = 0.5*(Cov_f+Cov_f.');
    mu_s = As*mup; Cov_s = As/Hp*As.'; Cov_s = 0.5*(Cov_s+Cov_s.');
    mu_fi = cell(N,1); Cov_fi = cell(N,1);
    mu_si = cell(N,1); Cov_si = cell(N,1);
    for i = 1:N
        mu_fi{i} = Afi{i}*mup; Cov_fi{i} = Afi{i}/Hp*Afi{i}.';
        Cov_fi{i} = 0.5*(Cov_fi{i}+Cov_fi{i}.');

        mu_si{i} = Asi{i}*mup; Cov_si{i} = Asi{i}/Hp*Asi{i}.';
        Cov_si{i} = 0.5*(Cov_si{i}+Cov_si{i}.');
    end

    % sample from the logN, calculate the prob densities, and 
    % construct empirical highest density mass functions
    samps_fs = exp(mvnrnd(mu_fs, Cov_fs, Nsamp));
    probs_fs = logmvnpdf(samps_fs, mu_fs(:)', Cov_fs);
    [probs_fs, sort_indx] = sort(probs_fs); samps_fs = samps_fs(sort_indx, :);

    samps_f = samps_fs(:,1:N); probs_f = logmvnpdf(samps_f, mu_f(:)', Cov_f);
    [probs_f, sort_indx] = sort(probs_f); samps_f = samps_f(sort_indx, :);

    samps_s = samps_fs(:,N+1:2*N); probs_s = logmvnpdf(samps_s, mu_s(:)', Cov_s);
    [probs_s, sort_indx] = sort(probs_s); samps_s = samps_s(sort_indx, :);

    samps_fi = cell(N,1); probs_fi = cell(N,1);
    samps_si = cell(N,1); probs_si = cell(N,1);
    for i = 1:N
        samps_fi{i} = samps_f(:,i);
        probs_fi{i} = logmvnpdf(samps_fi{i}, mu_fi{i}(:)', Cov_fi{i});
        [probs_fi{i}, sort_indx] = sort(probs_fi{i}); samps_fi{i} = samps_fi{i}(sort_indx, :);

        samps_si{i} = samps_s(:,i);
        probs_si{i} = logmvnpdf(samps_si{i}, mu_si{i}(:)', Cov_si{i});
        [probs_si{i}, sort_indx] = sort(probs_si{i}); samps_si{i} = samps_si{i}(sort_indx, :);
    end

    HDI_levels = (Nsamp:-1:1)'/Nsamp;
    
    % calculate probs of true diffusivity
    probs_true_fs = logmvnpdf([Df_true; Ds_true].', mu_fs(:).', Cov_fs);
    probs_true_f = logmvnpdf(Df_true.', mu_f(:).', Cov_f);
    probs_true_s = logmvnpdf(Ds_true.', mu_s(:).', Cov_s);
    probs_true_fi = cellfun(@(Df,mu,cov) logmvnpdf(Df, mu, cov), num2cell(Df_true), mu_fi, Cov_fi);
    probs_true_si = cellfun(@(Ds,mu,cov) logmvnpdf(Ds, mu, cov), num2cell(Ds_true), mu_si, Cov_si);

    % assign smallest HDI levels containing the true diffusivities
    Dfs_smallest_HDI_joint(trialNum) = interp_nearest(probs_fs, HDI_levels, probs_true_fs);
    Df_smallest_HDI_joint(trialNum) = interp_nearest(probs_f, HDI_levels, probs_true_f);
    Ds_smallest_HDI_joint(trialNum) = interp_nearest(probs_s, HDI_levels, probs_true_s);
    for i = 1:N
        Df_smallest_HDI_ptws(trialNum,i) = interp_nearest(probs_fi{i}, HDI_levels, probs_true_fi(i));
        Ds_smallest_HDI_ptws(trialNum,i) = interp_nearest(probs_si{i}, HDI_levels, probs_true_si(i));
    end
    
    disp("Progress: " + trialNum + "%")
end

%%
figure; hold on; axis equal;
plot([0;sort(Dfs_smallest_HDI_joint)*100], 0:100)
plot([0;sort(Df_smallest_HDI_joint)*100], 0:100)
plot([0;sort(Ds_smallest_HDI_joint)*100], 0:100)
plot([zeros(1,N);sort(Df_smallest_HDI_ptws)*100], 0:100)
plot([zeros(1,N);sort(Ds_smallest_HDI_ptws)*100], 0:100)
plot([0,100],[0,100],'--k', LineWidth=2)
xlim([0,100]); ylim([0,100])
xlabel("HDI Level [\%]", Interpreter="latex")
ylabel("\% of Infered HDIs that Captured the True Diffusivity", Interpreter="latex")

function p = logmvnpdf(x,mu,Sigma)
    assert(isequal(size(x,2),size(mu,2),size(Sigma,1),size(Sigma,2)));
    lnx = log(x);
    p = exp(log(mvnpdf(lnx,mu,Sigma)) - sum(lnx, 2));
end

function out = interp_nearest(x, y, xq)
    [~, indx] = min(abs(x - xq));
    out = y(indx);
end