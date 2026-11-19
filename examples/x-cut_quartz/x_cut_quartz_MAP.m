function x_cut_quartz_MAP(x0, options)
data = load("x-cut_quartz_data_stats.mat");

checkpointDir = "checkpoints";
if ~exist(checkpointDir, 'dir')
    mkdir(checkpointDir);
end

resultsDir = "results";
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

if nargin < 1 || isempty(x0)
    rng("shuffle");
    consts = load("x-cut_quartz_constants.mat");
    lnkf = normrnd(consts.mu_lnkf, consts.sigma_lnkf);
    lnCf = normrnd(consts.mu_lnCf, consts.sigma_lnCf);
    lnhf = normrnd(consts.mu_lnhf, consts.sigma_lnhf);
    lnks_perp = normrnd(consts.mu_lnks_perp, consts.sigma_lnks_perp);
    lnks_par = normrnd(consts.mu_lnks_par, consts.sigma_lnks_par);
    Os_theta = normrnd(0, deg2rad(10));
    Os_phi = 2*pi*rand;
    Os1 = cos(Os_theta) .* cos(Os_phi);
    Os2 = cos(Os_theta) .* sin(Os_phi);
    Os3 = sin(Os_theta);
    lnCs = normrnd(consts.mu_lnCs, consts.sigma_lnCs);
    lnkappaT = normrnd(consts.mu_lnkappaT, consts.sigma_lnkappaT);

    x0 = [lnkf; lnCf; lnhf; lnks_perp; lnks_par; lnCs; lnkappaT; Os1; Os2; Os3];
end

x0(8:10) = x0(8:10) / norm(x0(8:10));

if nargin < 2 || isempty(options)
    options = optimoptions("fmincon", Display="iter-detailed", ...
        SpecifyObjectiveGradient=true, ...
        SpecifyConstraintGradient=true, ...
        FiniteDifferenceType="central", ...
        StepTolerance=1e-10, FunctionTolerance=1e-10, ...
        MaxFunctionEvaluations=1e4, Algorithm="sqp");
end

obj_fun = @(x, x_max) format4optim( ...
    @(xi) nl_posterior( ...
        xi(1), xi(2), xi(3), xi(4), xi(5), xi(6), ...
        xi(7), xi(8), xi(9), xi(10), data, x_max ...
    ), x ...
);
nonlcon = @(x) format4constraint([], ...
    @(xi) hyper_sphere(xi(8), xi(9), xi(10)), ...
    x ...
);
nonlcon_check = @(x) format4optim( ...
    @(xi) hyper_sphere(xi(8), xi(9), xi(10)), ...
    x ...
);

x_max = calculate_x_max(x0(4), x0(5), x0(6), data.f, data.Xprobe);
[~, err] = checkGradients(@(x) obj_fun(x, x_max), x0, options, "Display","on");
[~, constraint_err] = checkGradients(nonlcon_check, x0, options, "Display","on");
[x,fval,exitflag,output,lambda] = fmincon( ...
    @(x) obj_fun(x, []), x0, [], [], [], [], [], [], nonlcon, options);

save(fullfile(checkpointDir, "x-cut_quartz_MAP_results_" ...
    + string(datetime("now", Format="uuuuMMdd'T'HHmmss")) ...
    + ".mat"), "x0", "x", "err", "fval", ...
    "exitflag", "output", "lambda", "constraint_err" ...
);
save(fullfile(resultsDir, "x-cut_quartz_MAP_results.mat"), ...
    "x0", "x", "err", "fval", ...
    "exitflag", "output", "lambda", "constraint_err" ...
);

if exitflag < 1
    gradM = NaN;
    B = NaN;
    hessM = NaN;
else
    % normalize x
    x(8:10) = x(8:10)/norm(x(8:10));

    x_max = calculate_x_max(x0(4), x0(5), x0(6), data.f, data.Xprobe);
    [fval, grad] = obj_fun(x, x_max);
    
    [~, Jc] = nonlcon_check(x);
    lam = (Jc*Jc') \ Jc * grad(:);
    B = blkdiag(eye(7), null(full(Jc(1, 8:10))));
    Hc = zeros(10,10);
    Hc(8:10,8:10) = 2*lam*eye(3);

    try
        assert(max(abs(Jc*B), [], "all") < 1e-8);
        assert(max(abs(B'*B-eye(52)), [], "all") < 1e-8);
        assert(max(abs(lam+lambda.eqnonlin)) < 1e-8);
        disp("PASS")
    catch exception
        disp("FAIL")
    end
    gradM = B'*grad(:);

    if max(abs(gradM)) > 1 && datetime('now') < datetime(2026, 11, 6, 8, 0, 0)
        hessM = NaN;
    else
        hessian = central_diff_hessian(@(x) obj_fun(x, x_max), x, 1e-4);
        hessM = B'*(hessian-Hc)*B;
    end
end

save(fullfile(checkpointDir, "isouniex_MAP_results_" ...
    + sprintf('%03d', trialNum) ...
    + "_" ...
    + string(datetime("now", Format="uuuuMMdd'T'HHmmss")) ...
    + ".mat"), "x0", "x", "err", "fval", ...
    "exitflag", "output", "lambda", "gradM", "hessM", "B", "constraint_err" ...
);
save(fullfile(resultsDir, "isouniex_MAP_results_" ...
    + sprintf('%03d', trialNum) ...
    + ".mat"), "x0", "x", "err", "fval", ...
    "exitflag", "output", "lambda", "gradM", "hessM", "B", "constraint_err" ...
);

if ~isequal(x, x0) && (isnan(hessM) || min(eig(hessM)) < 0 || max(abs(gradM)) > 1 || exitflag < 1)
    if exitflag < 1 || options.Algorithm == "interior-point"
        options.Algorithm = "sqp";
        options.MaxFunctionEvaluations=500;
    else
        options.Algorithm = "interior-point";
        options.MaxFunctionEvaluations=75;
    end
    isouniex_MAP(trialNum, x, options)
end
end

function x_max = calculate_x_max(lnks_perp, lnks_par, lnCs, f, Xprobe)
    Ds_perp = exp(get_val(lnks_perp)-get_val(lnCs)); % mm^2/s
    Ds_par = exp(get_val(lnks_par)-get_val(lnCs)); % mm^2/s
    Lths = sqrt(reshape(max(Ds_perp, Ds_par),1,1,[]) ./ pi ./ reshape(f,1,1,1,1,[])); % um
    x_max = max(5*max(sqrt(sum(Xprobe.^2, 2))), 2*Lths);
end

function psi = nl_posterior(lnkf, lnCf, lnhf, lnks_perp, lnks_par, lnCs, lnkappaT, Os1, Os2, Os3, data, x_max)
    %% DATA
    f = data.f;
    Xprobe = data.Xprobe;   
    Nx = 160;

    %% PRIOR PARAMS
    consts = load("x-cut_quartz_constants.mat");
    mu_lnkf         = consts.mu_lnkf         * ones(size(lnkf));
    sigma_lnkf      = consts.sigma_lnkf      * ones(size(lnkf));
    mu_lnCf         = consts.mu_lnCf         * ones(size(lnCf));
    sigma_lnCf      = consts.sigma_lnCf      * ones(size(lnCf));
    mu_lnhf         = consts.mu_lnhf         * ones(size(lnhf));
    sigma_lnhf      = consts.sigma_lnhf      * ones(size(lnhf));
    mu_lnks_perp    = consts.mu_lnks_perp    * ones(size(lnks_perp));
    sigma_lnks_perp = consts.sigma_lnks_perp * ones(size(lnks_perp));
    mu_lnks_par     = consts.mu_lnks_par     * ones(size(lnks_par));
    sigma_lnks_par  = consts.sigma_lnks_par  * ones(size(lnks_par));
    mu_lnCs         = consts.mu_lnCs         * ones(size(lnCs));
    sigma_lnCs      = consts.sigma_lnCs      * ones(size(lnCs));
    mu_lnkappaT     = consts.mu_lnkappaT     * ones(size(lnkappaT));
    sigma_lnkappaT  = consts.sigma_lnkappaT  * ones(size(lnkappaT));
    Qs              = consts.Qs;
    Okappas         = consts.Okappas;

    %% PRIOR FUNCTIONS

    % Ψ(lnkf)
    psi_lnkf = nln(lnkf, mu_lnkf, sigma_lnkf);
    
    % Ψ(lnCf)
    psi_lnCf = nln(lnCf, mu_lnCf, sigma_lnCf);
    
    % Ψ(lnks_perp)
    psi_lnks_perp = nln(lnks_perp, mu_lnks_perp, sigma_lnks_perp);

    % Ψ(lnks_par)
    psi_lnks_par = nln(lnks_par, mu_lnks_par, sigma_lnks_par);
    
    % Ψ(lnCs)
    psi_lnCs = nln(lnCs, mu_lnCs, sigma_lnCs);
    
    % Ψ(lnhf)
    psi_lnhf = nln(lnhf, mu_lnhf, sigma_lnhf);
    
    % Ψ(lnkappaT)
    psi_lnkappaT = nln(lnkappaT, mu_lnkappaT, sigma_lnkappaT);

    % Ψ(O)
    psi_Os = nl_Bingham_unnormalized([Os1;Os2;Os3], Qs, Okappas);

    %% OTHER DETERMINISTIC PARAMS
    lnaf = consts.lnaf;
    lnas = consts.lnas;
    logitRf = consts.logitRf;
    logitRs = consts.logitRs;
    lnRth = consts.lnRth;

    Sxx = data.Sxx;
    Syy = data.Syy;
    Sxy = data.Sxy;
    P = consts.P;

    if isempty(x_max)
        x_max = calculate_x_max(lnks_perp, lnks_par, lnCs, f, Xprobe);
    end

    %% LIKELIHOOD
    % Ψ(ϕ|x)
    psi_nll = cell(numel(data.Sxx),1);
    for i = 1:numel(data.Sxx)
        T0tilde = ForwardModel("iso", "uni", true).solve( ...
            {lnkf}, {}, lnCf, lnaf, logitRf, lnhf, ...
            {lnks_perp, lnks_par},  {Os1, Os2, Os3},  lnCs, lnas, logitRs, [], ...
            lnRth, Sxx(i), Syy(i), Sxy(i), P, f, x_max, Nx, Xprobe(i,:) ...
        );
        psi_nll{i} = nll_phase(data.mu(i,:,:,:), phase(T0tilde), exp(lnkappaT), data.kappa_mu(i,:,:,:), true);
    end
    
    %% POSTERIOR
    psi = sum_nl_probs( ...
        psi_lnkf, psi_lnCf, psi_lnhf, ...
        psi_lnks_perp, psi_lnks_par, psi_lnCs, psi_lnkappaT, ...
        psi_Os, psi_nll{:} ...
    );
end