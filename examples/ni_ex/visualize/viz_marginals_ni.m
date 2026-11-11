%% Main GP Results
clear; clc; close all;

%% Read in .mat files
% Read in MAP results
trialNum = "002";
ni_results = load("ni_MAP_results_"+trialNum+".mat");
constants = load("ni_constants.mat");
true_params = load("true_params.mat");
Tstruct = load("ni_data_stats.mat","T");
T = Tstruct.T;

num_T = length(T); % number of temperatures

% Indices
film_kC_idx = 1:2*num_T;
kf_idx = 1:num_T;
Cf_idx = num_T+1:2*num_T;
sub_kC_idx = 2*num_T+2:4*num_T+1;
ks_idx = 2*num_T+2:3*num_T+1;
Cs_idx = 3*num_T+2:4*num_T+1;
lnell_idx = 4*num_T+3: 4*num_T+6;

% Marginalize out film and substrate thicknesses 
indx = false(size(ni_results.hessian,1),1);
indx([film_kC_idx, sub_kC_idx, lnell_idx]) = true; % indices to keep

% Obtain marginal hessian with Schur complement
H_marg = schur_comp(ni_results.hessian, indx);

% Obtain the covariance
Sp = inv(H_marg); % invert hessian to obtain covariance
Sp = 0.5*(Sp+Sp.'); % symmetrize against round-off

% Mean of the distribution
mu = ni_results.x(indx);

% Prior mean of the distribution
mu0 = @(T) [constants.mu_lnkf, constants.mu_lnCf, constants.mu_lnks,...
    constants.mu_lnCs] .* ones(length(T),1);

% Prior standard deviation of the distribution
sigma0 = {...
    @(T) constants.sigma_lnkf .* ones(length(T),1);...
    @(T) constants.sigma_lnCf .* ones(length(T),1);...
    @(T) constants.sigma_lnks .* ones(length(T),1);...
    @(T) constants.sigma_lnCs .* ones(length(T),1);...
};

% Sample points for posterior samples
Ts = T;
N_samp = 500;

% sample posterior
ni_samples = sample_posterior_predictive(T, Ts, mu0, sigma0, mu, Sp, N_samp);
ni_samples = reshape(ni_samples, N_samp, [], 4);


%% Mode and S
% Indices for schur_comp call
indx = false(size(ni_results.hessian,1),1);
indx(film_kC_idx) = true;

% Marginalize out everything except kf and Cf
Hp = schur_comp(ni_results.hessian, indx); % compute marginal hessian
Sp = inv(Hp); % invert hessian to obtain covariance
Sp = 0.5*(Sp+Sp.'); % symmetrize round-off error

% Covariance of the diffusivity
Sp_Df = Sp(kf_idx,kf_idx) + Sp(Cf_idx,Cf_idx) - 2*Sp(kf_idx,Cf_idx);

% lnDf = ni_results.x(kf_idx) - ni_results.x(Cf_idx);

% Get mode and covariance
% [mode,S] = get_mode_S(lnDf, sqrt(diag(Sp_Df)));


%% Film diffusivity inference plot
figure(1); clf
hold on;
ni_D = 14.9;   % thermal diffusivity (mm^2/s)
graphite_D = 15.1; % thermal diffusivity (mm^2/s)

lnDf_samples = ni_samples(:,:,1)-ni_samples(:,:,2);
Df_samples = exp(lnDf_samples);

% Plot samples
histogram(exp(lnDf_samples),DisplayName="Posterior samples");
xline(ni_D,Color="r",LineWidth=2,DisplayName="True value")
legend
title("Ni-graphite synthetic data")
xlabel("Graphite (film) Diffusivity (mm^2/s)")
fontsize(18,"points")


%%
% Indicies for substrate k and C
indx = false(size(ni_results.hessian,1),1);
indx(sub_kC_idx) = true;

% Marginalize out everything except ks and Cs
Hp = schur_comp(ni_results.hessian, indx);
Sp = inv(Hp);
Sp = 0.5*(Sp+Sp.');

% Covariance of the diffusivity
Sp = Sp(kf_idx,kf_idx) + Sp(Cf_idx, Cf_idx)...
    - 2*Sp(kf_idx,Cf_idx);

% Get mode of S
[mode,S] = get_mode_S(ni_results.x(ks_idx)...
    -ni_results.x(Cs_idx), sqrt(diag(Sp)));

Ts_values = ones(size(ni_samples(:,:,3))) * Ts;

figure(2); clf
lnDs_samples = ni_samples(:,:,3) - ni_samples(:,:,4);
Ds_samples = exp(lnDs_samples);
plot_samp = histogram(Ds_samples, DisplayName="Posterior samples");
xline(ni_D,Color="r",LineWidth=2, DisplayName="True value")
legend
title("Ni-graphite synthetic data")
xlabel("Nickel (substrate) Diffusivity (mm^2/s)")
fontsize(18,"points")
%% Correlation Matrices

% Plot settings
bw = 0.05;
sc = 0.7;
gw = 0.25;
blw = 0.75;
grid_color = 0.5;
grid_bold = 0.1;
fig = figure;
fig.Units = 'inch';
fig.Position = [8 1 7 6.45];

% Correlation coefficient imagesc
figure(3)
tiledlayout(1,2,'TileSpacing','compact', "Padding","tight");

% k and C correlations
nexttile(1);

% Indices for schur_comp call
indx = false(size(ni_results.hessian,1),1);
indx([film_kC_idx, sub_kC_idx]) = true;
Sp = inv(schur_comp(ni_results.hessian,indx));
Sp = 0.5*(Sp+Sp.');
imagesc(corrcov(Sp))
clim([-1,1])

% Labels and ticks
D = size(Sp,1);
labels = "$\ln " + ["k_f", "C_f", "k_s", "C_s"]+"(\mathbf{T})$";
labels = labels(:);
yticks(1:num_T:D)
yticklabels(labels)
xticks(1:num_T:D)
xticklabels(labels)

% More plot settings
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0 0];
axis equal
xlim([0.5,D+0.5])
ylim([0.5,D+0.5])
title("Gold-YSZ Experiment", Interpreter="latex", FontSize=12)
ytickangle(90)
axis xy
ax.Clipping = 'off';
ax.XAxis.TickLabelGapMultiplier = sc;
ax.YAxis.TickLabelGapMultiplier = sc;

% Diffusivity correlations
nexttile(2);

% compute S_lnD and correlations
A = [eye(num_T),-eye(num_T),zeros(num_T,num_T),zeros(num_T,num_T);
    zeros(num_T,num_T), zeros(num_T,num_T), eye(num_T), eye(num_T)];
S_lnD = A*Sp*A.';
D = size(S_lnD,1);
ni_corr = corrcov(S_lnD);

% Plot diffusivity correlations
imagesc(ni_corr)
clim([-1,1])
c = colorbar;
ylabel(c, "Correlation Coefficient", Interpreter="latex", FontSize=11)

% Labels and ticks
labels = "$\ln " + ["D_f", "D_s"]+"(\mathbf{T})$";
labels = labels(:);
yticks(1:num_T:D)
yticklabels(labels)
xticks(1:num_T:D)
xticklabels(labels)
ytickangle(90)
ax = gca;
ax.TickLabelInterpreter = 'latex';
ax.TickLength = [0 0];
axis xy
ax.Clipping = 'off';
ax.XAxis.TickLabelGapMultiplier = sc;
ax.YAxis.TickLabelGapMultiplier = sc;
hold on
axis equal
xlim([0.5,D+0.5])
ylim([0.5,D+0.5])
c.Layout.Tile = "east";

function [mode, S] = get_mode_S(mu, sigma)
    mode = exp(mu - sigma.^2);
    S = zeros(size(sigma));
    for i = 1:length(sigma)
        S(i) = get_logN_HDI_param(sigma=sigma(i), alpha=0.95);
    end
end