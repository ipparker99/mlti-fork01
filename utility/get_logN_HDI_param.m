function out = get_logN_HDI_param(options)
arguments (Input)
    options.alpha (1,1) double
    options.S     (1,1) double
    options.sigma (1,1) double
end
arguments (Output)
    out (1,1) double {mustBePositive(out)}
end
is_alpha = isfield(options,'alpha');
is_S     = isfield(options,'S');
is_sigma = isfield(options,'sigma');
if sum([is_alpha, is_S, is_sigma])~=2
    error("Must provide exactly two of the three options (alpha, S, sigma).")
end
if ~is_alpha
    syms alph real positive
    v = alph;
    range = [eps,1];
    S = options.S;
    sigma_variable = options.sigma;
elseif ~is_S
    syms S real positive
    v = S;
    range = [1,1e300];
    alph = options.alpha;
    sigma_variable = options.sigma;
elseif ~is_sigma
    % Error message suggested changing syms to sym
    syms sigma_variable real positive
    v = sigma_variable;
    range = [eps, 1e100];
    alph = options.alpha;
    S = options.S;
else
    error("DEVELOPER ERROR: SHOULD NOT BE ABLE TO REACH THIS!")
end
equation = 2*alph == erf((log(S)-sigma_variable^2)/(sigma_variable*sqrt(2))) + erf((sigma_variable^2+log(S))/(sigma_variable*sqrt(2)));
out = vpasolve(equation, v, range);
end