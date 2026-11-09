function [X, Y, phase, mag] = forward_model(Kf11, Kf12, Kf13, Kf22,...
    Kf23, Kf33, Cf, af, Rf, hf, Ks11, Ks12, Ks13, Ks22, Ks23, Ks33,...
    Cs, as, Rs, hs, Rth, rx, ry, P, f0)
    
    Lth = sqrt(Ks11 / (pi * Cs * f0));    

    % Grid set up    
    factor = 30;
    x_max = factor * Lth;
    Nx = 2^10;
    [dx, X, Y, U, V] = build_grid(x_max, Nx);
    
    % compute T0hat(U,V)
    T0hat = T0hat_finite(Kf11, Kf12, Kf13, Kf22, Kf23, Kf33, Cf, ...
        af, Rf, hf, Ks11, Ks12, Ks13, Ks22, Ks23, Ks33, Cs, ...
        as, Rs, hs, Rth, rx, ry, P, f0, U, V);
    
    % convert T0hat(U,V) -> T0(X,Y)
    [X,Y,phase,mag] = T0hat_to_T0(X,Y,dx,Lth,T0hat);
end

