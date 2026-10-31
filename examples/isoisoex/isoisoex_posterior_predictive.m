clear; clc; close all;

colors = colororder;
color_true = colors(2,:);
color_post = colors(1,:);
color_samp = [0.5,0.5,0.5];
lw = 1.5; cs = 5; ms = 2.5;

Nsamp = 1000;

T = (300:100:900)';
Ts = (300:10:900)';

load results\isoisoex_MAP_results_100.mat
true_params = load("true_params.mat");

indx = [1:14,16:29,31:34];
Hp = schur_comp(hessian, indx);
mup = x(indx);
Sp = inv(Hp);
Sp = 0.5*(Sp+Sp');
post_samps = mvnrnd(mup, Sp, Nsamp);

consts = load("isoisoex_constants.mat");

pred_Df_samps = zeros(Nsamp, numel(Ts));
pred_Ds_samps = zeros(Nsamp, numel(Ts));

for i = 1:Nsamp
    K = blkdiag( ...
        RBFKernel(T,T,log(consts.sigma_lnkf), log(consts.sigma_lnkf), post_samps(i, 29)), ...
        RBFKernel(T,T,log(consts.sigma_lnCf), log(consts.sigma_lnCf), post_samps(i, 30)), ...
        RBFKernel(T,T,log(consts.sigma_lnks), log(consts.sigma_lnks), post_samps(i, 31)), ...
        RBFKernel(T,T,log(consts.sigma_lnCs), log(consts.sigma_lnCs), post_samps(i, 32)) ...
    );
    Ks = blkdiag( ...
        RBFKernel(T,Ts,log(consts.sigma_lnkf), log(consts.sigma_lnkf), post_samps(i, 29)), ...
        RBFKernel(T,Ts,log(consts.sigma_lnCf), log(consts.sigma_lnCf), post_samps(i, 30)), ...
        RBFKernel(T,Ts,log(consts.sigma_lnks), log(consts.sigma_lnks), post_samps(i, 31)), ...
        RBFKernel(T,Ts,log(consts.sigma_lnCs), log(consts.sigma_lnCs), post_samps(i, 32)) ...
    );
    Kss = blkdiag( ...
        RBFKernel(Ts,Ts,log(consts.sigma_lnkf), log(consts.sigma_lnkf), post_samps(i, 29)), ...
        RBFKernel(Ts,Ts,log(consts.sigma_lnCf), log(consts.sigma_lnCf), post_samps(i, 30)), ...
        RBFKernel(Ts,Ts,log(consts.sigma_lnks), log(consts.sigma_lnks), post_samps(i, 31)), ...
        RBFKernel(Ts,Ts,log(consts.sigma_lnCs), log(consts.sigma_lnCs), post_samps(i, 32)) ...
    );
    VsT = Ks' / K;

    MU0 = [consts.mu_lnkf, consts.mu_lnCf, consts.mu_lnks, consts.mu_lnCs] .* ones(numel(T), 1);
    mu0v = MU0(:);
    MU0s = [consts.mu_lnkf, consts.mu_lnCf, consts.mu_lnks, consts.mu_lnCs] .* ones(numel(Ts), 1);
    mu0vs = MU0s(:);

    mus = mu0vs + VsT * (post_samps(i,1:28).' - mu0v);
    Ss = Kss - VsT * Ks;
    Ss = 0.5 * (Ss + Ss.') + 1e-8*eye(size(Ss));
    pred_samp = mvnrnd(mus, Ss);
    pred_samp = reshape(pred_samp, [], 4);
    pred_Df_samps(i,:) = exp(pred_samp(:,1)-pred_samp(:,2));
    pred_Ds_samps(i,:) = exp(pred_samp(:,3)-pred_samp(:,4));
end

MUp = reshape(mup(1:28), [], 4);
mu_lnDfp = MUp(:,1) - MUp(:,2);
mu_lnDsp = MUp(:,3) - MUp(:,4);
sigma_lnDfp = sqrt(diag(Sp(1:7,1:7) + Sp(8:14, 8:14) - 2*Sp(1:7,8:14)));
sigma_lnDsp = sqrt(diag(Sp(15:21,15:21) + Sp(22:28, 22:28) - 2*Sp(15:21,22:28)));

mode_lnDfp = zeros(size(mu_lnDfp));
mode_lnDsp = zeros(size(mu_lnDsp));
mult_lnDfp = zeros(size(mu_lnDfp));
mult_lnDsp = zeros(size(mu_lnDsp));
for i = 1:length(mode_lnDfp)
    [mode_lnDfp(i), mult_lnDfp(i)] = mu_sigma_2_mode_mult(mu_lnDfp(i), sigma_lnDfp(i), 0.99);
    [mode_lnDsp(i), mult_lnDsp(i)] = mu_sigma_2_mode_mult(mu_lnDsp(i), sigma_lnDsp(i), 0.99);
end

figure(1); hold on;
patch([Ts; NaN]*ones(1,Nsamp), [pred_Df_samps'; nan(1, size(pred_Df_samps,1))], 'k', ...
'FaceColor', 'none', ...
'EdgeColor', 'k', ...
'EdgeAlpha', 0.1);

plot(Ts, polyval(true_params.gold_film.k, Ts)./polyval(true_params.gold_film.c, Ts)./polyval(true_params.gold_film.rho, Ts), '--r', LineWidth=2)

errorbar(T, mode_lnDfp, mode_lnDfp-mode_lnDfp./mult_lnDfp, mode_lnDfp-mode_lnDfp.*mult_lnDfp, LineWidth=lw, CapSize=cs, LineStyle="none", Marker='o', MarkerSize=ms, MarkerFaceColor='w', Color=color_post);

figure(2); hold on;
patch([Ts; NaN]*ones(1,Nsamp), [pred_Ds_samps'; nan(1, size(pred_Ds_samps,1))], 'k', ...
'FaceColor', 'none', ...
'EdgeColor', 'k', ...
'EdgeAlpha', 0.1);

figure(2);
plot(Ts, polyval(true_params.YSZ.k, Ts)./polyval(true_params.YSZ.c, Ts)./polyval(true_params.YSZ.rho, Ts), '--r', LineWidth=2)

errorbar(T, mode_lnDsp, mode_lnDsp-mode_lnDsp./mult_lnDsp, mode_lnDsp-mode_lnDsp.*mult_lnDsp, LineWidth=lw, CapSize=cs, LineStyle="none", Marker='o', MarkerSize=ms, MarkerFaceColor='w', Color=color_post);
