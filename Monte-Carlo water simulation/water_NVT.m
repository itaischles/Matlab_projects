clear all; close all;

% ===================
%     Initialize
% ===================

% Set configuration parameters
nPart = 7^2;                                    % Number of particles
density = 0.9;                                  % Density of particles

nSteps = 1e3;                                   % Total simulation time
printFreq = nPart;                              % Printing frequency

% Set energy struct parameters
EnStruct.eps_LJpair = 0.1;                      % Pair-particle Lennard-Jones interaction strength
EnStruct.eps_LJwall = 0.1;                      % Particle-wall Lennard-Jones interaction strength
EnStruct.eps_HB = -1;                           % Hydrogen-bond interaction strength
EnStruct.sigma_LJpair = 0.7;                    % Pair-particle Lennard-Jones distance scaling
EnStruct.sigma_LJwall = 0.7;                    % Particle-wall Lennard-Jones distance scaling
EnStruct.sigma_HB_r = 0.08;                     % Hydrogen-Bond radial Gaussian width
EnStruct.sigma_HB_t = 0.08;                     % Hydrogen-Bond angular Gaussian width
EnStruct.r_HB = 1.0;                            % Hydrogen-bond optimal bond distance
EnStruct.pairEnergyMat = zeros(nPart, nPart);   % Pair interaction energy matrix (each row\column corresponds to specific particle)
EnStruct.wallEnergyMat = zeros(1,nPart);        % Particle-wall interaction energy matrix (each column corresponds to specific particle)

% Set initial configuration
[coords, angles, L] = initSquareGrid(nPart, density);


% Set simulation parameters
Temp = 0.2;                                     % Simulation temperature
beta = 1.0/Temp;                                % Inverse temperature
maxDr = EnStruct.r_HB * 1.5;                    % Maximal displacement
maxDt = 45;                                     % maximal rotation (in degrees)
acc_rej = [0; 0];                               % no. of [accepted; rejected] moves

% Initialize the hydrogen bond visualization matrix. It has '1' for a pair
% where its hydrogen bond is stronger than its LJ term.
HBvisMat = zeros(nPart, nPart);

% Calculate initial energy
[EnStruct.pairEnergyMat, EnStruct.wallEnergyMat, HBvisMat] = LJ_HB_Energy(coords, angles, L, EnStruct, HBvisMat);


% ==============
% MC Simulation
% ==============

for i=1:nSteps
    
    % Choose a particle at random
    chosen = floor(rand*nPart)+1;
    
    % Suggest a trial move for particle i
    if (rand < 0.5)
        rTrial = coords(:, chosen) + maxDr*(2*rand(2,1)-1);
        tTrial = angles(chosen);
    else
        rTrial = coords(:, chosen);
        tTrial = wrapTo360(angles(chosen) + maxDt*(2*rand-1));
    end
    
    % Apply boundary conditions
%     rTrial = PBC2D(rTrial, L);
    rTrial = WallBC2D(rTrial, L);
    
    % Calculate the change in energy due to this trial move and the
    % corresponding row+column of interaction energies to put in the
    % energy matrix
    [deltaE, newPairInteractionRow, newWallEnergy, newHBvisMat] = LJ_HB_EnergyChange(coords, angles, rTrial, tTrial, chosen, L, EnStruct, HBvisMat);
    
    % Move according to Metropolis' algorithm
    if (rand < exp(-beta*deltaE))
        % Accept displacement move
        acc_rej = acc_rej + [1; 0];
        coords(:, chosen) = rTrial;                                 % Update positions
        angles(:, chosen) = tTrial;                                 % Update angles
        EnStruct.pairEnergyMat(nPart,:) = newPairInteractionRow;    % Update pair-interaction row
        EnStruct.pairEnergyMat(:,nPart) = newPairInteractionRow;    % Update pair-interaction column
        EnStruct.wallEnergyMat(nPart) = newWallEnergy;              % Update wall-interaction
        HBvisMat = newHBvisMat;                                     % Update hydrogen bond visualization matrix
    else
        acc_rej = acc_rej + [0; 1];
    end
    
    % print step number and draw particle configuration
    if (mod(i,printFreq)==0)
        
        disp(strcat({'Step '}, num2str(i), {' of '}, num2str(nSteps)))
        disp(strcat({'Accept/reject = '}, num2str(acc_rej(1)/acc_rej(2))))
        % set drawn particle radius to equal Lennard-Jones potential
        % minimum
        drawConfig(coords, angles, L, L, 2^(1/6)/2*EnStruct.sigma_LJpair, HBvisMat)

    end
    
end