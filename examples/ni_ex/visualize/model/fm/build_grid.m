function [dx,X,Y,U,V] = build_grid(x_max,Nx)

    arguments (Input)
        x_max   % maxium x position at which temperature is calculated
        Nx      % number of grid points along 1 dimension
    end
    
    arguments (Output)
        dx      % grid spacing between points
        X       % X values in meshgrid format
        Y       % Y values in meshgrid format
        U       % U values in meshgrid format
        V       % V value in meshgrid format
    end
    
    % calculate dx
    dx = x_max / floor(Nx/2);

    % build x, y, u, and v sample points
    x = dx * (-floor(Nx/2):ceil(Nx/2) - 1);
    y = dx * (-floor(Nx/2):ceil(Nx/2) - 1);
    u = (-floor(Nx/2):ceil(Nx/2) - 1) / (Nx * dx);
    v = (-floor(Nx/2):ceil(Nx/2) - 1) / (Nx * dx);
    
    
    % generate meshgrid of sample points
    [X, Y] = meshgrid(x, y); % size (Nx,Nx)
    [U, V] = meshgrid(u, v); % size (Nx,Nx)
end