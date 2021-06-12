function vec = WallBC2D(vec,L)

% Vector should be in the range 0 -> L in all dimensions except for the
% rigid walls. One wall is at y=0, the other one is the tip.
% Therefore, we need to apply the following changes if it's not in this range:

% Correct horizontal location
if (vec(1) > L)
    vec(1) = vec(1)-L;
elseif (vec(1) < 0)
    vec(1) = vec(1)+L;
end

% Correct vertical location
if (vec(2) > L)
    vec(2) = vec(2) - L;
end

end