function vec = PBC2D(vec,L)


% Vector should be in the range 0 -> L in all dimensions
% Therefore, we need to apply the following changes if it's not in this range:

for dim=1:2
    if (vec(dim) > L)
        vec(dim) = vec(dim)-L;
    elseif (vec(dim) < 0)
        vec(dim) = vec(dim)+L;
    end
    
end