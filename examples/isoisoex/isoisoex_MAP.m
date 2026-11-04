function isoisoex_MAP(trialNum, x0, options)
if nargin < 3 || isempty(options)
    options = optimoptions("fminunc", Display="iter-detailed", ...
        SpecifyObjectiveGradient=true, ...
        FiniteDifferenceType="central", ...
        StepTolerance=1e-10, FunctionTolerance=1e-10, ...
        MaxFunctionEvaluations=100, Algorithm="quasi-newton", MaxIterations=100);
end
data = load("isoisoex_data_stats.mat");
data.kappa_mu = data.kappa_mu{trialNum};
data.mu = data.mu{trialNum};

N = length(data.T);

checkpointDir = "checkpoints";
if ~exist(checkpointDir, 'dir')
    mkdir(checkpointDir);
end

resultsDir = "results";
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

if nargin < 2 || isempty(x0)
    %% LOAD DATA
    data_ws.f = data.f;
    data_ws.Xprobe = data.Xprobe;
    
    %% Initial Guess: x0 - Sample from Prior
    rng("shuffle");
    consts = load("isoisoex_constants.mat");
    lnkf = normrnd(consts.mu_lnkf, consts.sigma_lnkf, N, 1);
    lnCf = normrnd(consts.mu_lnCf, consts.sigma_lnCf, N, 1);
    lnhf = normrnd(consts.mu_lnhf, consts.sigma_lnhf, N, 1);
    lnks = normrnd(consts.mu_lnks, consts.sigma_lnks, N, 1);
    lnCs = normrnd(consts.mu_lnCs, consts.sigma_lnCs, N, 1);
    lnkappaT = normrnd(consts.mu_lnkappaT, consts.sigma_lnkappaT, 1, 1);
    lnell = normrnd(consts.mu_lnell, consts.sigma_lnell, 4, 1);
    
    %% WARM START: Temperature Independent
    err_ind = cell(N, 1);
    x_ind = cell(N,1);
    fval_ind = cell(N,1);
    exitflag_ind = cell(N,1);
    output_ind = cell(N,1);
    grad_ind = cell(N,1);
    hessian_ind = cell(N,1);
    for i = 1:N
        x0 = [lnkf(i); lnCf(i); lnhf; lnks(i); lnCs(i); lnkappaT];
        data_ws.T = data.T(i);
        data_ws.mu = data.mu(:,i,:,:);
        data_ws.kappa_mu = data.kappa_mu(:,i,:,:);
        obj_fun = @(x, x_max) format4optim( ...
            @(xi) nl_posterior( ...
                xi(1), xi(2), xi(3), ...
                xi(4), xi(5), xi(6), lnell, data_ws, x_max ...
            ), x ...
        );
        x_max = calculate_x_max(x0(4), x0(5), data_ws.f, data_ws.Xprobe);
        [~, err_ind{i}] = checkGradients(@(x) obj_fun(x, x_max), x0, options, "Display","on");
        [ ...
            x_ind{i}, fval_ind{i}, exitflag_ind{i}, ...
            output_ind{i}, grad_ind{i}, hessian_ind{i} ...
        ] = fminunc(@(x) obj_fun(x, []), x0, options);
    end

    %% WARM START: lnell
    x0 = [lnkf; lnCf; lnhf; lnks; lnCs; lnkappaT; lnell];
    x = horzcat(x_ind{:});
    x = x.';

    lnkf = x(:,1);
    lnCf = x(:,2);
    lnhf = mean(x(:,3));
    lnks = x(:,4);
    lnCs = x(:,5);
    lnkappaT = mean(x(:,6));

    obj_fun = @(x, x_max) format4optim( ...
        @(xi) nl_posterior( ...
            lnkf, lnCf, lnhf, ...
            lnks, lnCs, ...
            lnkappaT, xi(1:4), ...
            data, x_max ...
        ), x ...
    );

    x_max = calculate_x_max(lnks, lnCs, data.f, data.Xprobe);
    [~, err_lnell] = checkGradients(@(x) obj_fun(x, x_max), lnell, options, "Display","on");
    [ ...
        lnell, fval_lnell, exitflag_lnell, ...
        output_lnell, grad_lnell, hessian_lnell ...
    ] = fminunc(@(x) obj_fun(x, []), lnell, options);

    x = [lnkf; lnCf; lnhf; lnks; lnCs; lnkappaT; lnell];

    save(fullfile(checkpointDir, "isoisoex_MAP_results_warm_start_" ...
        + sprintf('%03d', trialNum) ...
        + "_" ...
        + string(datetime("now", Format="uuuuMMdd'T'HHmmss")) ...
        + ".mat"), "x0", "x", "err_ind", "x_ind", "fval_ind", ...
        "exitflag_ind", "output_ind", "grad_ind", "hessian_ind", ...
        "err_lnell", "fval_lnell", ...
        "exitflag_lnell", "output_lnell", "grad_lnell", "hessian_lnell" ...
    );
    save(fullfile(resultsDir, "isoisoex_MAP_results_" ...
        + sprintf('%03d', trialNum) ...
        + ".mat"), "x0", "x", "err_ind", "x_ind", "fval_ind", ...
        "exitflag_ind", "output_ind", "grad_ind", "hessian_ind", ...
        "err_lnell", "fval_lnell", ...
        "exitflag_lnell", "output_lnell", "grad_lnell", "hessian_lnell" ...
    );

    x0 = x;
end

obj_fun = @(x, x_max) format4optim( ...
    @(xi) nl_posterior( ...
        xi(1:N), ...
        xi(N+1:2*N), ...
        xi(2*N+1), ...
        xi(2*N+2:3*N+1), ...
        xi(3*N+2:4*N+1), ...
        xi(4*N+2), ...
        xi(4*N+3:4*N+6), ...
        data, x_max ...
    ), x ...
);

x_max = calculate_x_max(x0(2*N+2:3*N+1), x0(3*N+2:4*N+1), data.f, data.Xprobe);
[~, err] = checkGradients(@(x) obj_fun(x, x_max), x0, options, "Display","on");
[x,fval,exitflag,output,grad] = fminunc(@(x) obj_fun(x, []), x0, options);

if exitflag < 1
    hessian = NaN;
else
    x_max = calculate_x_max(x(2*N+2:3*N+1), x(3*N+2:4*N+1), data.f, data.Xprobe);
    hessian = central_diff_hessian(@(x) obj_fun(x, x_max), x, 1e-4);
end

save(fullfile(checkpointDir, "isoisoex_MAP_results_" ...
    + sprintf('%03d', trialNum) ...
    + "_" ...
    + string(datetime("now", Format="uuuuMMdd'T'HHmmss")) ...
    + ".mat"), "x0", "x", "err", "fval", ...
    "exitflag", "output", "grad", "hessian" ...
);
save(fullfile(resultsDir, "isoisoex_MAP_results_" ...
    + sprintf('%03d', trialNum) ...
    + ".mat"), "x0", "x", "err", "fval", ...
    "exitflag", "output", "grad", "hessian" ...
);

if ~isequal(x, x0) && (isnan(hessian) || min(eig(hessian)) < 0 || max(abs(grad)) > 1 || exitflag < 1)
    % if options.Algorithm == "trust-region"
    %     options.Algorithm = "quasi-newton";
    %     options.MaxIterations=100;
    % else
    %     options.Algorithm = "trust-region";
    %     options.MaxIterations=10;
    % end
    isoisoex_MAP(trialNum, x, options)
end
end

function x_max = calculate_x_max(lnks, lnCs, f, Xprobe)
    Ds = exp(get_val(lnks)-get_val(lnCs)); % mm^2/s
    Lths = sqrt(reshape(Ds,1,1,[]) ./ pi ./ reshape(f,1,1,1,1,[])); % um
    x_max = max(5*max(sqrt(sum(Xprobe.^2, 2))), 2*Lths);
end

function psi = nl_posterior(lnkf, lnCf, lnhf, lnks, lnCs, lnkappaT, lnell, data, x_max)
    %% DATA
    T = data.T;
    f = data.f;
    Xprobe = data.Xprobe; 
    Nx = 160;

    %% PRIOR PARAMS
    consts = load("isoisoex_constants.mat");
    mu_lnkf        = consts.mu_lnkf        * ones(size(lnkf));
    sigma_lnkf     = consts.sigma_lnkf     * ones(size(lnkf));
    mu_lnCf        = consts.mu_lnCf        * ones(size(lnCf));
    sigma_lnCf     = consts.sigma_lnCf     * ones(size(lnCf));
    mu_lnhf        = consts.mu_lnhf        * ones(size(lnhf));
    sigma_lnhf     = consts.sigma_lnhf     * ones(size(lnhf));
    mu_lnks        = consts.mu_lnks        * ones(size(lnks));
    sigma_lnks     = consts.sigma_lnks     * ones(size(lnks));
    mu_lnCs        = consts.mu_lnCs        * ones(size(lnCs));
    sigma_lnCs     = consts.sigma_lnCs     * ones(size(lnCs));
    mu_lnkappaT    = consts.mu_lnkappaT    * ones(size(lnkappaT));
    sigma_lnkappaT = consts.sigma_lnkappaT * ones(size(lnkappaT));
    mu_lnell       = consts.mu_lnell       * ones(size(lnell));
    sigma_lnell    = consts.sigma_lnell    * ones(size(lnell));

    %% PRIOR FUNCTIONS
    % Ψ(θ)
    psi_lnell = nln(lnell, mu_lnell, sigma_lnell, true);

    % Ψ(lnkf|θ1)
    psi_lnkf = nl_mvn_cov(lnkf, mu_lnkf, ...
        RBFKernel(T,T,log(sigma_lnkf),log(sigma_lnkf),lnell(1)));
    
    % Ψ(lnCf|θ2)
    psi_lnCf = nl_mvn_cov(lnCf, mu_lnCf, ...
        RBFKernel(T,T,log(sigma_lnCf),log(sigma_lnCf),lnell(2)));
    
    % Ψ(lnks|θ3)
    psi_lnks = nl_mvn_cov(lnks, mu_lnks, ...
        RBFKernel(T,T,log(sigma_lnks),log(sigma_lnks),lnell(3)));
    
    % Ψ(lnCs|θ4)
    psi_lnCs = nl_mvn_cov(lnCs, mu_lnCs, ...
        RBFKernel(T,T,log(sigma_lnCs),log(sigma_lnCs),lnell(4)));
    
    % Ψ(lnhf)
    psi_lnhf = nln(lnhf, mu_lnhf, sigma_lnhf);
    
    % Ψ(lnkappaT)
    psi_lnkappaT = nln(lnkappaT, mu_lnkappaT, sigma_lnkappaT);
    
    %% OTHER DETERMINISTIC PARAMS
    lnaf = consts.lnaf;
    lnas = consts.lnas;
    logitRf = consts.logitRf;
    logitRs = consts.logitRs;
    lnRth = consts.lnRth;

    sx = consts.sx;
    sy = consts.sy;
    P = consts.P;

    if isempty(x_max)
        x_max = calculate_x_max(lnks, lnCs, f, Xprobe);
    end

    %% LIKELIHOOD
    % Ψ(ϕ|x)
    T0tilde = ForwardModel("iso", "iso", true).solve( ...
        {lnkf}, {}, lnCf, lnaf, logitRf, lnhf, ...
        {lnks},  {},  lnCs, lnas, logitRs, [], ...
        lnRth, sx^2, sy^2, 0, P, f, x_max, Nx, Xprobe ...
    );
    psi_nll = nll_phase(data.mu, phase(T0tilde), exp(lnkappaT), data.kappa_mu, true);
    
    %% POSTERIOR
    psi = sum_nl_probs( ...
        psi_lnell, psi_lnkf, psi_lnCf, psi_lnhf, ...
        psi_lnks, psi_lnCs, psi_lnkappaT, psi_nll ...
    );
end