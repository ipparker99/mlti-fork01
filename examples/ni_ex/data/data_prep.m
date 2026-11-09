clear; clc;

% Read in simulation data
filename = "ni_lockin_ref_data.mat";
data_ref = load(filename, "X", "Y", "phase");

% Simulation data notes
% X, size (1225,1), units (m)
% Y, size (1225,1), units (m)
% phase, size(1225,1), units (deg)

%% Interpolate reference phase
nq = 20;
xq = linspace(min(data_ref.X), max(data_ref.X), nq);
yq = linspace(min(data_ref.Y), max(data_ref.Y), nq);
[Xq,Yq] = meshgrid(xq,yq);

% Interpolate phase
data_ref.phase = deg2rad(data_ref.phase); % convert degrees to radians
phase_interpolant = scatteredInterpolant(data_ref.X, data_ref.Y,...
    data_ref.phase);
phase_q = phase_interpolant(Xq(:),Yq(:));

% set temperature (K)
T = 300; 

% Set frequency (MHz)
f = [15, 15]*1E-6; % scalar frequency causes an error, quick fix

% Combine x-y coordinates into one array
Xprobe = [Xq(:), Yq(:)]; % size (N_probe,2)
Xprobe = Xprobe*1e6; % convert m to um

% Make phase_ref the appropriate size
phase_ref(:,:,1,1) = phase_q;
phase_ref(:,:,1,2) = phase_q;

% Save noise free reference
save("ni_data_noise_free.mat", "T", "f", "Xprobe", "phase_ref")

%% Add noise to reference

% Noise settings
Nrep = 100;
Ntrials = 10;

% Initialize cell arrays
phi_obs     = cell(Ntrials, 1);
sigmaT_true = cell(Ntrials, 1);
sigmaD_low  = cell(Ntrials, 1);
sigmaD_hi   = cell(Ntrials, 1);
sigmaD_true = cell(Ntrials, 1);
subs = repmat({':'}, 1, ndims(phase_ref));

% Noise adding loop
for i = 1:Ntrials
    sigmaT_true{i} = deg2rad(unifrnd(0.1,5));
    sigmaD_low{i} = randi(6);
    sigmaD_hi{i}  = sigmaD_low{i} + randi(4);
    sigmaD_true{i} = deg2rad(unifrnd(sigmaD_low{i},sigmaD_hi{i},size(phase_ref)));

    phi_true = phase_ref + normrnd(0,sigmaT_true{i},size(phase_ref));

    phi_obs{i} = zeros([size(phase_ref), Nrep]);
    for j = 1:Nrep
        phi_obs{i}(subs{:},j) = phi_true + normrnd(0,sigmaD_true{i},size(phi_true));
    end
end

save("ni_data_corrupted.mat", "T", "f", "Xprobe", "phi_obs",...
    "sigmaT_true", "sigmaD_low", "sigmaD_hi", "sigmaD_true");

%% Compute summary statistics from noisy data 
% Ntrials = numel(phi_obs);
mu       = cell(Ntrials,1);
Rbar     = cell(Ntrials,1);
kappa    = cell(Ntrials,1);
kappa_mu = cell(Ntrials,1);
for a = 1:Ntrials
    Nrep = size(phi_obs{a},ndims(phi_obs{a}));
    mu{a} = circular_mean(phi_obs{a}, ndims(phi_obs{a}));
    Rbar{a} = R(phi_obs{a}, ndims(phi_obs{a}))./Nrep;
    kappa{a} = get_kappa(Rbar{a});
    kappa_mu{a} = Nrep.*Rbar{a}.*kappa{a};
    disp("Progress: " + a + "/" + Ntrials)
end

save("ni_data_stats.mat", "mu", "Rbar", "kappa", "kappa_mu", "T", "f", "Xprobe")

function out = circular_mean(x, dim)
    s = sin(x);
    c = cos(x);
    out = atan2(sum(s, dim), sum(c, dim));
end

function out = R(x, dim)
    s = sum(sin(x), dim);
    c = sum(cos(x), dim);
    out = sqrt(c.^2+s.^2);
end

function kappas = get_kappa(Rbar)
    kappas = zeros(numel(Rbar), 1);
    for i = 1:numel(Rbar) 
        kappas(i) = fzero(@(x) besseli(1, x, 1)./besseli(0, x, 1) - Rbar(i), 100);
    end
    kappas = reshape(kappas, size(Rbar));
end