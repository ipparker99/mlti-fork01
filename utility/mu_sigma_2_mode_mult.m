function [mode, S] = mu_sigma_2_mode_mult(mu, sigma, alpha)
    mode = exp(mu-sigma.^2);
    S = get_logN_HDI_param(sigma=sigma, alpha=alpha);
end