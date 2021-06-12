clear all; close all;

% ===================
%     Initialize
% ===================

% Set configuration parameters
nPart = 14^2;       % Number of particles
density = 0.85;     % Density of particles

% Set simulation parameters
Temp = 4.0;         % Simulation temperature
beta = 1.0/Temp;    % Inverse temperature
press = 2.0;        % Simulation pressure
maxDr = 0.4;        % Maximal displacement
maxDv = 0.01;       % Maximal ln volume displacement

nSteps = 100;       % Total simulation time (in integration steps)
printFreq = 1;      % Printing frequency

% Set initial configuration
[coords, L] = initSquareGrid(nPart,density);

% Start up an empty histogram for radial distribution function
h.count = 0;
h.range = [0, sqrt(2)*L/2];
h.increment = L/100;

% Calculate initial potential energy
energy = LJ_Energy(coords,L);


% ==============
% MC Simulation
% ==============
for step = 1:nSteps
    
    % Choose whether to do a Volume move or a Displacement move
    
    if (rand*(nPart+1) + 1 < nPart)
        
        % Preform particle displacement
        
        for i=1:nPart
            % Suggest a trial move for particle i
            rTrial = coords(:,i) + maxDr*(rand(2,1)-0.5);
            
            % Apply periodic boundary conditions
            rTrial = PBC2D(rTrial,L);
            
            % Calculate the change in energy due to this trial move
            deltaE = LJ_EnergyChange(coords,rTrial,i,L);
            
            if (rand < exp(-beta*deltaE))
                % Accept displacement move
                coords(:,i) = rTrial;       % Update positions
                energy = energy + deltaE;   % Update energy
            end

        end
        
    else
        
        % Preform volume change move
        
        oldV = L^2;
        
        % Suggest a random volume by preforming a random walk in lnV.
        % See Frenkel & Smit 5.4.1 & 5.4.2 for reference.
        lnvTrial = log(oldV) + (rand - 0.5)*maxDv;
        vTrial = exp(lnvTrial);     % Trial volume
        newL = vTrial^(1.0/2);      % Trial box length
        
        % Rescale all coordinates
        coordsTrial = coords*(newL/L);
        
        % Calculate the energy of the scaled coordinates
        eTrial = LJ_Energy(coordsTrial, newL);
        
        % Calculate the weight function for the acceptance
        weight = (eTrial - energy) + press*(vTrial - oldV) - (nPart+1)*Temp*log(vTrial/oldV);
        
        if (rand < exp(-beta*weight))
            % Volume move accepted
            coords = coordsTrial;   % Update coordinates
            energy = eTrial;        % Update energy
            L = newL;               % Update box size
        end
        
    end
    
    h = RDF(h, coords, L);
    
    if (mod(step,printFreq)==0)
        disp(strcat({'Step '}, num2str(step), {' of '}, num2str(nSteps)))
        drawConfig(coords, L, L, 0.5)
    end
    
% %     % plot the radial distribution function
% %     cla;
% %     bar(h.values, h.histo)
% %     drawnow;
    
end