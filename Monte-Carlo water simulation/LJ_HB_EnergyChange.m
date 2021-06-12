function [deltaE, newPairInteractionRow, newWallEnergy, newHBvisMat] = LJ_HB_EnergyChange(coords, angles, trialPos, trialAng, part, L, EnStruct, HBvisMat)

% Get the number of particles
nPart = size(coords,2);

% Extract energy parameters
eps_LJpair = EnStruct.eps_LJpair;           % Pair-particle Lennard-Jones interaction strength
eps_LJwall = EnStruct.eps_LJwall;           % Particle-wall Lennard-Jones interaction strength
eps_HB = EnStruct.eps_HB;                   % Hydrogen-bond interaction strength
sigma_LJpair = EnStruct.sigma_LJpair;       % Pair-particle Lennard-Jones distance scaling
sigma_LJwall = EnStruct.sigma_LJwall;       % Particle-wall Lennard-Jones distance scaling
sigma_HB_r = EnStruct.sigma_HB_r;           % Hydrogen-Bond radial Gaussian width
sigma_HB_t = EnStruct.sigma_HB_t;           % Hydrogen-Bond angular Gaussian width
r_HB = EnStruct.r_HB;                       % Hydrogen-bond optimal bond distance
pairEnergyMatOld = EnStruct.pairEnergyMat;  % Old Energy matrix of the pair-particle interactions (before trial move)
wallEnergyMatOld = EnStruct.wallEnergyMat;  % Old Energy matrix of the pair-particle interactions (before trial move)

% Initialize pair interaction row for selected particle
newPairInteractionRow = zeros(1, nPart);

% Initialize the new hydrogen bond visualization matrix to the old matrix
newHBvisMat = HBvisMat;

% Get the number of particles
nPart = size(coords,2);

% Calculate new energy of particle part with the wall:

    % Calculate distance to the walls
    rWallnew = trialPos(2);
    
    % Get the sigma_LJ over distance
    sigma_over_rWall_new = sigma_LJwall/rWallnew;
    
    % calculate interaction energy with the wall for particle 'part'
    newWallEnergy = 4*eps_LJwall*(sigma_over_rWall_new^12 - sigma_over_rWall_new^6);

    
    
% Loop over all particles and calculate interaction with particle
% 'part'.
for otherPart = 1:nPart
    
    % Make sure to skip particle 'part' so that we don't calculate self
    % interaction
    if (otherPart == part)
        continue
    end
    
    % Calculate particle-particle distance for the new configurations
    rNew = coords(:, otherPart) - trialPos;
    
    % Fix according to periodic boundary conditions
    rNew = distPBC2D(rNew, L);
    
    % Calculate the absolute value of the distance
    rAbsNew = sqrt(rNew'*rNew);
    
    % Get the sigma_LJ over distance
    sigma_over_r_New = sigma_LJpair/rAbsNew;
    
    % Get angles for the particle
    tNew = degtorad(trialAng);
    tOther = degtorad(angles(:, otherPart));
    
    % Find unit vector connecting two particles
    uRNew = rNew/rAbsNew;
    
    % Find unit vectors of each arm of the molecules.
    % Top row is x component, bottom row is y component.
    % Each column (out of the three) is a single unit vector.
    uANew = [[cos(tNew), cos(tNew + 2*pi/3), cos(tNew + 4*pi/3)]; [sin(tNew), sin(tNew + 2*pi/3), sin(tNew + 4*pi/3)]];
    uOther = [[cos(tOther), cos(tOther + 2*pi/3), cos(tOther + 4*pi/3)]; [sin(tOther), sin(tOther + 2*pi/3), sin(tOther + 4*pi/3)]];
    
    
    % Lennard-Jones potential:
    % U(r) = 4*eps_LJ* [(sigma_LJ/r)^12 - (sigma_LJ/r)^6]
    %      = 4*eps_LJ*(sigma_LJ/r)^6* [(sigma_LJ/r)^2 - 1]
    
    LJ_New = 4*eps_LJpair*(sigma_over_r_New^12 - sigma_over_r_New^6);
    
    % Hydrogen-bond potential:
    % U(r, angle1, angle2) =~
    % eps_HB* sum( Gaussian(r)*Gaussian(angle1)*Gaussian(angle2) )
    
    HB_New = eps_HB*gaussian(rAbsNew, r_HB, sigma_HB_r)*sum(gaussian(uRNew'*uANew, 1, sigma_HB_t))*sum(gaussian(uRNew'*uOther, -1, sigma_HB_t));
   
    
    newPairInteractionRow(otherPart) = LJ_New + HB_New;
    
    if (abs(HB_New) > abs(LJ_New))
        newHBvisMat(part, otherPart) = 1;
        newHBvisMat(otherPart, part) = 1;
    else
        newHBvisMat(part, otherPart) = 0;
        newHBvisMat(otherPart, part) = 0;
    end
    
end

deltaE = sum(newPairInteractionRow - pairEnergyMatOld(part,:)) + (newWallEnergy - wallEnergyMatOld(part));

end