function deltaE = LJ_EnergyChange(coords, trialPos, part, L)

deltaE = 0;

% Get the number of particles
nPart = size(coords,2);

% Loop over all particles and calculate interaction with particle
% 'part'.
for otherPart = 1:nPart
    
    % Make sure to skip particle 'part' so that we don't calculate self
    % interaction
    if (otherPart == part)
        continue
    end
    
    % Calculate particle-particle distance for both the old and new
    % configurations
    drNew = coords(:,otherPart) - trialPos;
    drOld = coords(:,otherPart) - coords(:,part);
    
    % Fix according to periodic boundary conditions
    drNew = distPBC2D(drNew,L);
    drOld = distPBC2D(drOld,L);
    
    % Get the distance squared
    dr2_New = sum(dot(drNew,drNew));
    dr2_Old = sum(dot(drOld,drOld));
    
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
    
    invDr6_New = 1.0/(dr2_New^3); % 1/r^6
    invDr6_Old = 1.0/(dr2_Old^3); % 1/r^6
    
    % Calculate the potential energy
    eNew = (invDr6_New * (invDr6_New - 1));
    eOld = (invDr6_Old * (invDr6_Old - 1));
    
    deltaE = deltaE + eNew - eOld;
end

% Multiply energy by 4
deltaE = deltaE*4;

end