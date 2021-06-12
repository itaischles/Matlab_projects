function plotTripod( coords, angle, radius )
% Gets the coordinates of a single particle and its angle of rotation and
% plots a tripod (3 legs) coordinated inside the circle denoting the
% particle at the desired rotation angle

% calculate the x and y vectors for each leg of the tripod
x1 = [coords(1), coords(1)+radius*cos(degtorad(angle))];
y1 = [coords(2), coords(2)+radius*sin(degtorad(angle))];
x2 = [coords(1), coords(1)+radius*cos(degtorad(angle+120))];
y2 = [coords(2), coords(2)+radius*sin(degtorad(angle+120))];
x3 = [coords(1), coords(1)+radius*cos(degtorad(angle+240))];
y3 = [coords(2), coords(2)+radius*sin(degtorad(angle+240))];

% plot the tripod
line(x1, y1, 'Color', 'r');
line(x2, y2, 'Color', 'r');
line(x3, y3, 'Color', 'r');

end

