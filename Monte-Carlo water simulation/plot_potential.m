eps_LJ = 0.1;
eps_HB = -1;
sigma_LJ = 0.7;
sigma_HB_r = 0.08;
sigma_HB_t = 0.08;
r_HB = 1.0;

tA = 0;
tB = linspace(0, 2*pi, 100);
uR = [1; 0];
uA = [[cos(tA), cos(tA + 2*pi/3), cos(tA + 4*pi/3)]; [sin(tA), sin(tA + 2*pi/3), sin(tA + 4*pi/3)]];
uB = cell(1,length(tB));
angular_gaussian = cell(1,length(tB));
HB_term = cell(1,length(tB));
energy = cell(1,length(tB));

r = linspace(0.68,3,100);

LJ_term = 4*eps_LJ* ((sigma_LJ./r).^12-(sigma_LJ./r).^6);
 

for i=1:length(tB)
        uB{i} = [[cos(tB(i)), cos(tB(i) + 2*pi/3), cos(tB(i) + 4*pi/3)]; [sin(tB(i)), sin(tB(i) + 2*pi/3), sin(tB(i) + 4*pi/3)]];
        angular_gaussian{i} = sum(gaussian(uR'*uA-1, 0, sigma_HB_t))*sum(gaussian(uR'*uB{i}+1, 0, sigma_HB_t));
end
angular_gaussian = cell2mat(angular_gaussian);

for j = 1:length(angular_gaussian)
    HB_term{j} = eps_HB*gaussian(r, r_HB, sigma_HB_r)*angular_gaussian(j);
    energy{j} = HB_term{j} + LJ_term;
end

energy = cell2mat(energy(:));

[R,TB] = meshgrid(r,tB);

X = R.*cos(TB);
Y = R.*sin(TB);

contour(X, Y, energy)
% surf(X, Y, energy)
axis square