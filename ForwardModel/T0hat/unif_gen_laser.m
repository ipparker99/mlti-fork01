function unif_gen_laser()
syms Kf11 Kf12 Kf13 Kf22 Kf23 Kf33 lnCf lnaf logitRf lnhf
syms Ks11 Ks12 Ks13 Ks22 Ks23 Ks33 lnCs lnas logitRs lnhs
syms lnRth sxx syy sxy P f u v
syms AA BB [6,1]
syms LAMBDA [2,1]
syms E [4,1]

input_vars = [Kf11 Kf12 Kf13 Kf22 Kf23 Kf33 lnCf lnaf logitRf lnhf ...
              Ks11 Ks12 Ks13 Ks22 Ks23 Ks33 lnCs lnas logitRs lnhs ...
              lnRth sxx syy sxy P f u v];

Cf = exp(lnCf); Cs = exp(lnCs);
af = exp(lnaf); as = exp(lnas);
Rf = 1/(1+exp(-logitRf)); Rs = 1/(1+exp(-logitRs));
hf = exp(lnhf); hs = exp(lnhs);
Rth = exp(lnRth);

FF = @(AAi, BBi) [AAi;  AAi * BBi;  (Rth * AAi + 1) * BBi;  0];
SS = @(AAi, BBi) [0;  -AAi;  -1;  AAi * BBi];

M4 = [FF(AA1, BB1), FF(AA2, BB2), FF(AA3, BB3), ...
      SS(AA4, BB4), SS(AA5, BB5), SS(AA6, BB6)];
C4 = [E(1:2);LAMBDA(1);E(3:4);LAMBDA(2)];
sol4 = solve(M4*C4==0, E);
Ep_finite = matlabFunction(sol4.E1, Vars=[AA; BB; LAMBDA; lnRth]);
Em_finite = matlabFunction(sol4.E2, Vars=[AA; BB; LAMBDA; lnRth]);

M3 = M4(1:3,[1:4,6]);
C3 = C4([1:4,6]);
sol3 = solve(M3*C3==0, E(1:3));
Ep_infinite = matlabFunction(sol3.E1, Vars=[AA; BB; LAMBDA; lnRth]);
Em_infinite = matlabFunction(sol3.E2, Vars=[AA; BB; LAMBDA; lnRth]);

% Laser distribution
% Change gaussian distribution
%phuv = exp(-2.*pi.^2.*(u.^2.*sxx + 2.*u.*v.*sxy + v.^2.*syy));

% Uniform (top-hat distribution)
a = sxx; % laser radius in x-direction
b = syy; % laser radius in y-direction
tau = sqrt((a*u).^2 + (b*v).^2);
phuv = a * b * besselj(1,2*pi*tau) / (2 * pi * tau);

psi  = 0.25;

I0f = P .* (1 - Rf);
I0s = P .* (1 - Rf) .* (1 - Rs) .* exp(-af .* hf);

lf = 2 .* pi .* (u .* Kf13 + v .* Kf23);
ls = 2 .* pi .* (u .* Ks13 + v .* Ks23);

gf = 4 .* pi.^2 .* ...
    (u.^2 .* Kf11 + 2 .* u .* v .* Kf12 + v.^2 .* Kf22) ...
    + pi * 2i .* f .* Cf;
gs = 4 .* pi.^2 .* ...
    (u.^2 .* Ks11 + 2 .* u .* v .* Ks12 + v.^2 .* Ks22) ...
    + pi * 2i .* f .* Cs;

epf = (1i .* lf + sqrt(gf .* Kf33 - lf.^2)) ./ -Kf33;
emf = (1i .* lf - sqrt(gf .* Kf33 - lf.^2)) ./ -Kf33;
eps = (1i .* ls + sqrt(gs .* Ks33 - ls.^2)) ./ -Ks33;
ems = (1i .* ls - sqrt(gs .* Ks33 - ls.^2)) ./ -Ks33;

Lf = phuv .* psi .* af .* I0f ...
    ./ (gf + 2i .* lf .* af - Kf33 .* af.^2);
Ls = phuv .* psi .* as .* I0s ...
    ./ (gs + 2i .* ls .* as - Ks33 .* as.^2);

AA_f_ep = 1i .* lf + Kf33 .* epf;
AA_f_em = 1i .* lf + Kf33 .* emf;
AA_f_ma = 1i .* lf - Kf33 .* af;
AA_s_ep = 1i .* ls + Ks33 .* eps;
AA_s_em = 1i .* ls + Ks33 .* ems;
AA_s_ma = 1i .* ls - Ks33 .* as;

BB_f_ep = exp(hf .* epf);
BB_f_em = exp(hf .* emf);
BB_f_ma = exp(hf .* -af);
BB_s_ep = exp(hs .* eps);
BB_s_em = exp(hs .* ems);
BB_s_ma = exp(hs .* -as);

AA = {AA_f_ep;AA_f_em;AA_f_ma;AA_s_ep;AA_s_em;AA_s_ma};
BB = {BB_f_ep;BB_f_em;BB_f_ma;BB_s_ep;BB_s_em;BB_s_ma};
LAMBDA = {Lf,Ls};

T0hat_finite   = Ep_finite(AA{:}, BB{:}, LAMBDA{:}, lnRth) ...
    + Em_finite(AA{:}, BB{:}, LAMBDA{:}, lnRth) + Lf;
T0hat_infinite = Ep_infinite(AA{:}, BB{:}, LAMBDA{:}, lnRth) ...
    + Em_infinite(AA{:}, BB{:}, LAMBDA{:}, lnRth) + Lf;

grad_T0hat_finite = cell(length(input_vars),1);
grad_T0hat_infinite = cell(length(input_vars),1);

for i = 1:length(input_vars)
    grad_T0hat_finite{i} = gradient(T0hat_finite, input_vars(i));
    grad_T0hat_infinite{i} = gradient(T0hat_infinite, input_vars(i));
end

matlabFunction(T0hat_finite, grad_T0hat_finite{:}, ...
    Vars=input_vars, ...
    Optimize=true, ...
    File='unif_finite_auto_optimized', ...
    Outputs=["T0hat", "grad_" + string(input_vars)], ...
    Comments="Created by solve_system_symbolically()");
matlabFunction(T0hat_infinite, grad_T0hat_infinite{:}, ...
    Vars=input_vars, ...
    Optimize=true, ...
    File='unif_infinite_auto_optimized', ...
    Outputs=["T0hat", "grad_" + string(input_vars)], ...
    Comments="Created by solve_system_symbolically()");
end