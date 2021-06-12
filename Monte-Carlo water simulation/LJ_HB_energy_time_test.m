clear all; close all;

% ===================
%     Initialize
% ===================

% Set configuration parameters
nPart = round(linspace(10,500,10));           % Number of particles
density = 0.85;                                % Density of particles
EnStruct.eps_LJ = 0.1;                      % Lennard-Jones interaction strength
EnStruct.eps_HB = 1;                        % Hydrogen-bond interaction strength
EnStruct.sigma_LJ = 0.7;                    % Lennard-Jones distance scaling
EnStruct.sigma_HB = 0.085;                  % Hydrogen-Bond Gaussian width
EnStruct.r_HB = 1;                          % Hydrogen-bond optimal bond distance
EnStruct.energyTot = 0;                     % Sum of all interaction pairs


t = zeros(1, length(nPart));
for i=1:length(nPart)
    
    EnStruct.energyMat = zeros(nPart(i), nPart(i));   % Interaction energy matrix
    
    % Set initial configuration
    [coords, angles, L] = initSquareGrid(nPart(i), density);
    [EnStruct.energyTot, EnStruct.energyMat]  = LJ_HB_Energy(coords, angles, L, EnStruct);
    
    % Calculate energy
    tic
    deltaE_OldFunc(i) = LJ_HB_EnergyChange(coords, angles, [1;1], 45, 1, L, EnStruct);
    tOld(i) = toc;
    
    tic
    deltaE_NewFunc(i) = LJ_HB_EnergyChange_NEW(coords, angles, [1;1], 45, 1, L, EnStruct);
    tNew(i) = toc;
    
    disp(num2str(i))
end

figure()
plot(nPart, tOld, nPart, tNew)
xlabel('Particle number')
ylabel('Energy calculation time [sec]')
legend('Old Function','New Function')

figure()
plot(nPart, deltaE_NewFunc-deltaE_OldFunc)
xlabel('Particle number')
ylabel('Calculation difference')