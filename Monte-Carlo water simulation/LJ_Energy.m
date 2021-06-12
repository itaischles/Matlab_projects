function energy = LJ_Energy(coords,L)

energy = 0;

% Get the number of particles
nPart = size(coords,2);

% Loop over all distinct particle pairs
for partA = 1:nPart-1
    for partB = (partA+1):nPart
        
        % Calculate particle-particle distance
        dr = coords(:,partA) - coords(:,partB);
        % Fix according to periodic boundary conditions
        dr = distPBC2D(dr,L);
        % Get the distance squared
        dr2 = sum(dot(dr,dr));
        
        % Lennard-Jones potential:
        % U(r) = 4*epsilon* [(sigma/r)^12 - (sigma/r)^6]
        %
        % Here, we set sigma = 1, epsilon = 1 (reduced distance and
        % energy units). Therefore:
        %
        % U(r) = 4 * [(1/r)^12 - (1/r)^6]
        %
        % For efficiency, we will multiply by 4 only after summing
        % up all the energies.
        
        invDr6 = 1.0/(dr2^3); % 1/r^6
        energy = energy + (invDr6 * (invDr6 - 1));
        
    end
end

% Multiply energy by 4
energy = energy*4;

end