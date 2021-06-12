function [coords, angles, L] = initSquareGrid(nPart, density)

% Initialize with zeroes
coords = zeros(2, nPart);
angles = zeros(1, nPart);

% Get the cooresponding box size
L = sqrt(nPart/density);

% Find the lowest perfect square greater than or equal to the number of
% particles
nSquare = 2;

while (nSquare^2 < nPart)
    nSquare = nSquare + 1;
end


% Start positioning - use a 2D index for counting the spots
index = [0,0]';

% Assign particle positions
for part=1:nPart
    % Set coordinate
    coords(:,part) = (index+[0.5,0.5]')*(L/nSquare);
    
    % Advance the index
    index(1) = index(1) + 1;
    if (index(1) == nSquare)
        index(1) = 0;
        index(2) = index(2) + 1;
    end
end

end