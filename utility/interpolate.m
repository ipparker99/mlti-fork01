function out = interpolate(x,y,Fin,Xprobe,method)

arguments
    x
    y
    Fin
    Xprobe (:,2) double
    method = "linear"
end

F = get_val(Fin);

% Ensure each dimension of x and y either matches F or is singleton
for d = 1:ndims(F)
    assert(size(x,d) == size(F,d) || size(x,d) == 1, ...
        "the length of each dimension of x must equal " + ...
        "that of F or be a singleton dimension");

    assert(size(y,d) == size(F,d) || size(y,d) == 1, ...
        "the length of each dimension of y must equal " + ...
        "that of F or be a singleton dimension");
end

assert(size(x,2) == 1, ...
    "the second dimension of x must be a singleton dimension")

assert(size(y,1) == 1, ...
    "the first dimension of y must be a singleton dimension")

sheetSize = size(F);
sheetSize = sheetSize(3:end);

out = zeros([size(Xprobe,1), sheetSize], 'like', F);

if isa(Fin, "DesignVariable")
    outJac = zeros([size(Xprobe,1), sheetSize, Fin.rootLen], 'like', F);
    Jac = reshape(Fin.Jac, [size(F), Fin.rootLen]);
end

subs = cell(1,numel(sheetSize));

for i = 1:prod(sheetSize)

    if isscalar(sheetSize)
        subs = {i};
    else
        [subs{:}] = ind2sub(sheetSize,i);
    end

    % x and y may be singleton along sheet dimensions
    subsx = subs;
    subsy = subs;

    for j = 1:numel(subs)
        subsx{j} = min(size(x,j+2),subsx{j});
        subsy{j} = min(size(y,j+2),subsy{j});
    end

    Fi = F(:,:,subs{:});
    xi = x(:,1,subsx{:});
    yi = y(1,:,subsy{:});

    if ~all(isfinite(Fi),"all") || ...
       ~all(isfinite(xi),"all") || ...
       ~all(isfinite(yi),"all")
        out(:,subs{:}) = NaN;
    else
        out(:,subs{:}) = interp2( ...
            xi,yi,Fi.', ...
            Xprobe(:,1),Xprobe(:,2),method, NaN);
    end

    if isa(Fin,"DesignVariable")
        for j = 1:Fin.rootLen
            Jacij = Jac(:,:,subs{:},j);
            if ~all(isfinite(Jacij), "all") || ...
               ~all(isfinite(xi),"all") || ...
               ~all(isfinite(yi),"all")
                outJac(:, subs{:}, j) = NaN;
            else
                outJac(:,subs{:},j) = interp2( ...
                    xi,yi,Jacij.', ...
                    Xprobe(:,1),Xprobe(:,2),method, NaN);
            end
        end
    end
end

if isa(Fin,"DesignVariable")
    outJac = reshape(outJac,[],Fin.rootLen);
    out = DesignVariable(out,[],Fin.rootLen,outJac);
end

end