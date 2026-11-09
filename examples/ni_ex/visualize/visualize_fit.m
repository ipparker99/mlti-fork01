clear
load("tiny_MAP_results_001.mat")
consts = load("isoisoex_constants.mat");
load("ni_data_stats.mat","f","Xprobe")

f = 15;
% forward_model()
% 
% obj_fun = @(x, x_max) format4optim( ...
%             @(xi) nl_posterior( ...
%                 xi(1), xi(2), xi(3), ...
%                 xi(4), xi(5), xi(6), lnell, data_ws, x_max ...
%             ), x ...
%         );
% 


lnkf = x(1);
lnCf = x(2);
lnhf = x(3);
lnks = x(4);
lnCs = x(5);

% OTHER DETERMINISTIC PARAMS
lnaf = consts.lnaf;
lnas = consts.lnas;
logitRf = consts.logitRf;
logitRs = consts.logitRs;
lnRth = consts.lnRth;

sx = consts.sx;
sy = consts.sy;
P = consts.P;

Nx = 160;

x_max = calculate_x_max(lnks, lnCs, f, Xprobe);

T0tilde = ForwardModel("iso", "iso", true).solve( ...
        {lnkf}, {}, lnCf, lnaf, logitRf, lnhf, ...
        {lnks},  {},  lnCs, lnas, logitRs, [], ...
        lnRth, sx^2, sy^2, 0, P, f, x_max, Nx, Xprobe ...
    );

function x_max = calculate_x_max(lnks, lnCs, f, Xprobe)
    Ds = exp(get_val(lnks)-get_val(lnCs)); % mm^2/s
    Lths = sqrt(reshape(Ds,1,1,[]) ./ pi ./ reshape(f,1,1,1,1,[])); % um
    x_max = max(5*max(sqrt(sum(Xprobe.^2, 2))), 2*Lths);
end

% 
% function x_max = calculate_x_max(lnks, lnCs, f)
%     Ds = exp(lnks-lnCs); % mm^2/s
%     Lths = sqrt(Ds ./ pi ./ f); % um
%     x_max = 30*Lths;
% end