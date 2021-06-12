function [ h ] = RDF_bulk( h, coords, L )

nPart = size(coords,2); % number of particles

% calculate RDF for each particle, 'part'
for part = 1:nPart-1
    
    % iterate over surrounding particles counting each pair once
    for surroundPart = part+1:nPart
        
%         % do not count the particle for which the RDF is calculated
%         if (surroundPart == part)
%             continue
%         end
        
        % Calculate particle-particle distance vector
        rvec = coords(:,surroundPart) - coords(:,part);
        % Fix according to periodic boundary conditions
        rvec = distPBC2D(rvec,L);
        % Get the distance vector's absolute value
        r = sqrt(rvec'*rvec); 
        
        %insert distance to histogram (normalize by the shell volume later)
        h = histogram(h, r);
        
    end
    
end

% normalize histogram by shells' volumes and by the average density and by
% the number of particles
% Notice: because of periodic boundary conditions the apperant system volume is smaller by a factor of 2
density = nPart/(2*L^2);
h.histo = (h.histo / density) ./ (2*pi*h.values*h.increment) / nPart;

end

