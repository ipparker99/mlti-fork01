clear

%addpath(fullfile("..","fm"))

% model configuration -----------------------------------------------------------
% Film properties
Kf = 24;         % thermal conductivity (W/m/K)
Df = 15.1E-06;   % thermal diffusivity (m^2/s)
Cf = Kf/Df;      % volumetric heat capacity (J/m^3/K)
af = 3.4385E+07; % optical absorption coefficient (1/m)
Rf = 0.32;       % optical reflectivity, R in [0,1]
hf = 2e-6;       % thickness (m)


% Film thermal conductivity tensor (isotropic configuration)
Kf11 = Kf; Kf12 =  0; Kf13 =  0;
Kf21 =  0; Kf22 = Kf; Kf23 =  0;
Kf31 =  0; Kf32 =  0; Kf33 = Kf;


% Substrate properties
Ks = 60.7;       % thermal conductivity (W/m/K)
Ds = 14.9E-06;   % thermal diffusivity (m^2/s)
Cs = Ks/Ds;      % volumetric heat capacity (J/m^3/K)
as = 8.2547E+07; % optical absorption coefficient (1/m)
Rs = 0.63;       % optical reflectivity
hs = 1.5e-3;     % thickness (m)


% Substrate thermal condutivity tensor (isotropic configuration)
Ks11 = Ks; Ks12 =  0; Ks13 =  0;
Ks21 =  0; Ks22 = Ks; Ks23 =  0;
Ks31 =  0; Ks32 =  0; Ks33 = Ks;


% Heat source properties
lzang = deg2rad(125); % laser angle with the x-axis
ry = 100E-06; % laser radius in y direction (m)
rx = 122E-06; % laser radius in x-direction due to lzang
P = 1; % laser power (W)
f0 = 15; % modulation frequency (Hz)

% Thermal contact resistance
Rth = 0; % no contact resistance
% End of configuration ---------------------------------------------------


% Thermal diffusion length calculated from substrate (m)
Lth  = sqrt(Ds/(pi*f0));


[X, Y, phase, mag] = forward_model(Kf11,Kf12,Kf13,Kf22,Kf23,Kf33,Cf,af,Rf,hf,...
    Ks11,Ks12,Ks13,Ks22,Ks23,Ks33,Cs,as,Rs,hs,...
    Rth,rx,ry,P,f0);


figure(1)
clf
scatter3(X,Y,phase,80,phase,"filled")
xlabel X
ylabel Y
zlabel("Phase (deg)")
colorbar
fontsize(18,"points")
savefig(fullfile("data","phase_plot.fig"))

figure(2)
clf
scatter3(X,Y,mag,80,mag,"filled")
xlabel X
ylabel Y
zlabel("Magnitude")
fontsize(18,"points")
colorbar
savefig(fullfile("data","mag_plot.fig"))

save(fullfile("data","data.mat"),"X","Y","phase","mag")

% function [X, Y, phase, mag] = forward_model(Kf11, Kf12, Kf13, Kf22,...
%     Kf23, Kf33, Cf, af, Rf, hf, Ks11, Ks12, Ks13, Ks22, Ks23, Ks33,...
%     Cs, as, Rs, hs, Rth, rx, ry, P, f0)
% 
%     Lth = sqrt(Ks11 / (pi * Cs * f0));    
% 
%     % Grid set up
%     factor = 30;
%     x_max = factor * Lth;
%     Nx = 2^10;
%     [dx, X, Y, U, V] = build_grid(x_max, Nx);
% 
%     % compute T0hat(U,V)
%     T0hat = T0hat_finite(Kf11, Kf12, Kf13, Kf22, Kf23, Kf33, Cf, ...
%         af, Rf, hf, Ks11, Ks12, Ks13, Ks22, Ks23, Ks33, Cs, ...
%         as, Rs, hs, Rth, rx, ry, P, f0, U, V);
% 
%     % convert T0hat(U,V) -> T0(X,Y)
%     [X,Y,phase,mag] = T0hat_to_T0(X,Y,dx,Lth,T0hat);
% end