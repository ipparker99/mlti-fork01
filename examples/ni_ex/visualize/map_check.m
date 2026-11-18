clear

%% Read in .mat files

% Map estimate
map = load("ni_MAP_results_002.mat","x");
lnkf = map.x(1);
lnCf = map.x(2);
lnhf = map.x(3);
lnks = map.x(4);
lnCs = map.x(5);

% % Known inputs
% lnkf = log(24.0);
% lnCf = log((24.0 / 15.1E-06) / 100^3);
% lnhf = log(2);
% lnks = log(60.7);
% lnCs = log((60.7 / 14.9E-06) / 100^3);
% lnhs = log(1500);

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
Nx = 1000;
x_max = calculate_x_max(lnks, lnCs, f, Xprobe);

T0tilde = ForwardModel("iso", "iso", false).solve( ...
        {lnkf}, {}, lnCf, lnaf, logitRf, lnhf, ...
        {lnks},  {},  lnCs, lnas, logitRs, lnhs, ...
        lnRth, sx, sy, 0, P, f, x_max, Nx, Xprobe ...
    );

phase = angle(T0tilde); % obtain phase of temperature
phase = rad2deg(phase); % convert phase from rad to deg

%% Figures
% Plot difference between reference and predicted
figure(1); clf
diff = ref.phase_interp - phase(:,1,1,1);
scatter3(Xprobe(:,1),Xprobe(:,2),diff,[],diff,'filled',MarkerFaceAlpha=0.5)
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
colorbar
fontsize(18,"points")

%% Helper code
% Method for choosing x_max
function x_max = calculate_x_max(lnks, lnCs, f, Xprobe)
    Ds = exp(get_val(lnks)-get_val(lnCs)); % mm^2/s
    f = f * 1E+06; % MHz to Hz
    Lths = sqrt(reshape(Ds,1,1,[]) ./ pi ./ reshape(f,1,1,1,1,[])); % m
    x_max = max(30*Lths);
    x_max = x_max * 1E+03; % mm to um
end