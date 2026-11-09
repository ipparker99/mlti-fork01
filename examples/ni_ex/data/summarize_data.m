clear; clc;

load ni_data_corrupted.mat
Ntrials = numel(phi_obs);
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