clear; clc; rng(22);

load("ni_data_noise_free.mat")

Nrep = 20;

Ntrials = 10;

phi_obs     = cell(Ntrials, 1);
sigmaT_true = cell(Ntrials, 1);
sigmaD_low  = cell(Ntrials, 1);
sigmaD_hi   = cell(Ntrials, 1);
sigmaD_true = cell(Ntrials, 1);
subs = repmat({':'}, 1, ndims(phi_COMSOL));
for i = 1:Ntrials
    sigmaT_true{i} = deg2rad(unifrnd(0.1,5));
    sigmaD_low{i} = randi(6);
    sigmaD_hi{i}  = sigmaD_low{i} + randi(4);
    sigmaD_true{i} = deg2rad(unifrnd(sigmaD_low{i},sigmaD_hi{i},size(phi_COMSOL)));

    phi_true = phi_COMSOL + normrnd(0,sigmaT_true{i},size(phi_COMSOL));

    phi_obs{i} = zeros([size(phi_COMSOL), Nrep]);
    for j = 1:Nrep
        phi_obs{i}(subs{:},j) = phi_true + normrnd(0,sigmaD_true{i},size(phi_true));
    end
end

save("ni_data_corrupted.mat", "T", "f", "Xprobe", "phi_obs",...
    "sigmaT_true", "sigmaD_low", "sigmaD_hi", "sigmaD_true");