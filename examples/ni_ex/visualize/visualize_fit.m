clear

%% Read in .mat files

% Map estimate
map = load("ni_MAP_results_002.mat","x");
lnkf = map.x(1);
lnCf = map.x(2);
lnhf = map.x(3);
lnks = map.x(4);
lnCs = map.x(5);

% Constants
consts = load("ni_constants.mat");
lnaf = consts.lnaf;
lnas = consts.lnas;
logitRf = consts.logitRf;
logitRs = consts.logitRs;
lnRth = consts.lnRth;
sx = consts.sx;
sy = consts.sy;
P = consts.P;

% Experiment (f and Xprobe)
experiment = load("ni_data_stats.mat","f","Xprobe");
f = experiment.f;
Xprobe = experiment.Xprobe; % units of um

% Reference data
ref = load("ni_lockin_ref_data.mat","X","Y","phase");
X = ref.X * 1e6; % convert m to um
Y = ref.Y * 1e6; % convert m to um
ref.interpolant = scatteredInterpolant(X, Y, ref.phase);
ref.phase_interp = ref.interpolant(Xprobe(:,1),Xprobe(:,2));
%% Forward model call
Nx = 160;
x_max = calculate_x_max(lnks, lnCs, f, Xprobe);

T0tilde = ForwardModel("iso", "iso", true).solve( ...
        {lnkf}, {}, lnCf, lnaf, logitRf, lnhf, ...
        {lnks},  {},  lnCs, lnas, logitRs, [], ...
        lnRth, sx^2, sy^2, 0, P, f, x_max, Nx, Xprobe ...
    );

phase = angle(T0tilde); % obtain phase of temperature
phase = rad2deg(phase); % convert phase from rad to deg

%% Figures
% Plot difference between reference and predicted
figure(1); clf
diff = ref.phase_interp - phase(:,1,1,1);
scatter3(Xprobe(:,1),Xprobe(:,2),diff)
xlabel("X (um)")
ylabel("Y (um)")
zlabel("Phase--ref minus predicted (deg)")
fontsize(18,"points")

% Plot parity
figure(2); clf; hold on
plot(ref.phase_interp, phase(:,1,1,1),'.')
plot([-80,0],[-80,0])
xlabel("Reference phase")
ylabel("Predicted phase")
fontsize(18,"points")

%% Helper code
% Method for choosing x_max
function x_max = calculate_x_max(lnks, lnCs, f, Xprobe)
    Ds = exp(get_val(lnks)-get_val(lnCs)); % mm^2/s
    Lths = sqrt(reshape(Ds,1,1,[]) ./ pi ./ reshape(f,1,1,1,1,[])); % um
    x_max = max(5*max(sqrt(sum(Xprobe.^2, 2))), 2*Lths);
end