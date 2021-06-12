function [pairEnergyMat, wallEnergyMat, HBvisMat] = LJ_HB_Energy(coords, angles, L, EnStruct, HBvisMat)

% Extract energy parameters
eps_LJpair = EnStruct.eps_LJpair;           % Pair-particle Lennard-Jones interaction strength
eps_LJwall = EnStruct.eps_LJwall;           % Particle-wall Lennard-Jones interaction strength
eps_HB = EnStruct.eps_HB;                   % Hydrogen-bond interaction strength
sigma_LJpair = EnStruct.sigma_LJpair;       % Pair-particle Lennard-Jones distance scaling
sigma_LJwall = EnStruct.sigma_LJwall;       % Particle-wall Lennard-Jones distance scaling
sigma_HB_r = EnStruct.sigma_HB_r;           % Hydrogen-Bond radial Gaussian width
sigma_HB_t = EnStruct.sigma_HB_t;           % Hydrogen-Bond angular Gaussian width
r_HB = EnStruct.r_HB;                       % Hydrogen-bond optimal bond distance
pairEnergyMat = EnStruct.pairEnergyMat;     % Energy matrix of the pair-particle interactions
wallEnergyMat = EnStruct.wallEnergyMat;     % Energy matrix of the particle-wall interactions


% Get the number of particles
nPart = size(coords,2);

% Loop over all particles and calculate interaction with the walls
for part = 1:nPart
    
    % Calculate distance to the walls
    rWall = coords(2,part);
    
    % Get the sigma_LJ over distance
    sigma_over_rWall = sigma_LJwall/rWall;
    
    % calculate interaction energy with the wall for particle 'part'
    wallEnergyMat(part) = 4*eps_LJwall*(sigma_over_rWall^12 - sigma_over_rWall^6);
    
end

% Loop over all distinct particle pairs
for partA = 1:nPart-1
    for partB = (partA+1):nPart
        
        % Calculate particle-particle distance vector
        r = coords(:, partA) - coords(:, partB);
        % Fix according to periodic boundary conditions
        r = distPBC2D(r, L);
        % Calculate the absolute value of the distance
        rAbs = sqrt(r'*r);
        % Get the sigma_LJ over distance
        sigma_over_r = sigma_LJpair/rAbs;
        % Get angles for each particle
        tA = degtorad(angles(:, partA));
        tB = degtorad(angles(:, partB));
        % Find unit vector connecting two particles
        uR = r/rAbs;
        % Find unit vectors of each arm of the molecules.
        % Top row is x component, bottom row is y component.
        % Each column (out of the three) is a single unit vector.
        uA = [[cos(tA), cos(tA + 2*pi/3), cos(tA + 4*pi/3)]; [sin(tA), sin(tA + 2*pi/3), sin(tA + 4*pi/3)]];
        uB = [[cos(tB), cos(tB + 2*pi/3), cos(tB + 4*pi/3)]; [sin(tB), sin(tB + 2*pi/3), sin(tB + 4*pi/3)]];
        
        
        % Lennard-Jones potential:
        % U(r) = 4*eps_LJ* [(sigma_LJ/r)^12 - (sigma_LJ/r)^6]
        %      = 4*eps_LJ*(sigma_LJ/r)^6* [(sigma_LJ/r)^2 - 1]
        
        LJ_term = 4*eps_LJpair*(sigma_over_r^12 - sigma_over_r^6);
        
        
        % Hydrogen-bond potential:
        % U(r, angle1, angle2) =~
        % eps_HB* sum( Gaussian(r)*Gaussian(angle1)*Gaussian(angle2) )
        
        HB_term = eps_HB*gaussian(rAbs, r_HB, sigma_HB_r)*sum(gaussian(uR'*uA-1, 0, sigma_HB_t))*sum(gaussian(uR'*uB+1, 0, sigma_HB_t));
        
        
        % Place the energy of the specific particle pair in the correct
        % places in the energy matrix (symmetric matrix)
        pairEnergyMat(partA, partB) = LJ_term + HB_term;
        pairEnergyMat(partB, partA) = LJ_term + HB_term;
        
        if (abs(HB_term) > abs(LJ_term))
            HBvisMat(partA, partB) = 1;
            HBvisMat(partB, partA) = 1;
        else
            HBvisMat(partA, partB) = 0;
            HBvisMat(partB, partA) = 0;
        end

    end
end

end