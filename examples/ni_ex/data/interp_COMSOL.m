clear; clc; clf;

% Set up ------------
% T0tilde, 4D complex double, (4766, 7, 1, 17)
filename = "data_COMSOL.mat";



% Read in .mat file
data_COMSOL = load(filename, "X", "Y", "phase");
% --------------------

% Interpolation points
nq = 20;
xq = linspace(min(data_COMSOL.X), max(data_COMSOL.X), nq);
yq = linspace(min(data_COMSOL.Y), max(data_COMSOL.Y), nq);
[Xq,Yq] = meshgrid(xq,yq);
% Interpolate phase
phase_interpolant = scatteredInterpolant(data_COMSOL.X, data_COMSOL.Y,data_COMSOL.phase);
phase_q = phase_interpolant(Xq(:),Yq(:));

% Extract first temperature and first frequency
T = 300;
f = [15, 15];


% Combine x-y coordinates into one array
Xprobe = [Xq(:), Yq(:)];

phi_COMSOL(:,:,1,1) = phase_q;
phi_COMSOL(:,:,1,2) = phase_q;

save("ni_data_noise_free.mat", "T", "f", "Xprobe", "phi_COMSOL")