function [X,Y,phase,mag] = T0hat_to_T0(X,Y,dx,Lth,T0hat)
%FUNC_T0HAT_TO_T0 Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    X       % X-grid
    Y       % Y-grid
    dx      % Sampling interval 
    Lth     % thermal diffusion length
    T0hat   % Fourier space temperature
end

arguments (Output)
    X % X-grid for -Lth < X < Lth
    Y % Y-grid for -Lth < Y < Lth
    phase % Phase lag (deg)
    mag % Temperature magnitude
end


% ifftshift
T0hat_shift = ifftshift(T0hat);

% remove nans
T0hat_shift(isnan(T0hat_shift)) = 0;

% ifft2 and fftshift
T0 = fftshift(ifft2(T0hat_shift))/dx^2;

% Check for nans
if any(isnan(T0))
    warning("nan value(s) in T0")
end

% Obtain indices withing region of interest
idx = X < Lth & X > -Lth & Y < Lth & Y > -Lth;

% Crop arrays to region of interest
X = X(idx);
Y = Y(idx);
T0 = T0(idx);

% Compute phase lag and magnitude
phase = rad2deg(angle(T0));
mag = abs(T0);
end