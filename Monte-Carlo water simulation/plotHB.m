function plotHB( coords, HBvisMat, L)
% Plot the hydrogen bonds of the pairs which have stronger HB than LJ



% find the number of particels
nPart = size(coords,2);

% Loop over pairs of particles
for i = 1:nPart-1
    for j = i+1:nPart
        
        xi = coords(1,i);
        xj = coords(1,j);
        yi = coords(2,i);
        yj = coords(2,j);
        
        if (HBvisMat(i,j) == 1)
            % Check if the distance between the particles is less than L/2
            if ( sqrt((xi-xj)^2+(yi-yj)^2) < L/2 )
                line([xi, xj], [yi, yj], 'Color', 'g');
            end
        end
        
    end
end

end

