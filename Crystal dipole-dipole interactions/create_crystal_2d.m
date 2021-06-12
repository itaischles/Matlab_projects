function [XA,YA,XB,YB,muXA,muYA,muXB,muYB] = create_crystal_2d(input_params)

crystal_size = input_params{1};
Ax = input_params{2};
Ay=input_params{3};
dA=input_params{4};
dB=input_params{5};
muA=input_params{6};
muB=input_params{7};
thetaA=input_params{8};
thetaB=input_params{9};

xa = linspace(0,Ax*(crystal_size(1)),crystal_size(1)+1);
ya = linspace(0,Ay*(crystal_size(2)),crystal_size(2)+1);
xb = linspace(0,Ax*(crystal_size(1)),crystal_size(1)+1);
yb = linspace(0,Ay*(crystal_size(2)),crystal_size(2)+1);

xa = xa-max(xa)/2+dA(1);
ya = ya-max(ya)/2+dA(2);
xb = xb-max(xb)/2+dB(1);
yb = yb-max(yb)/2+dB(2);

[XA,YA] = meshgrid(xa,ya);
[XB,YB] = meshgrid(xb,yb);

muXA = muA*cos(thetaA)*ones(crystal_size+1)';
muYA = muA*sin(thetaA)*ones(crystal_size+1)';
muXB = muB*cos(thetaB)*ones(crystal_size+1)';
muYB = muB*sin(thetaB)*ones(crystal_size+1)';

% add random noise
noise = 0;
muXA = muXA+noise*randn(size(muXA));
muYA = muYA+noise*randn(size(muYA));
muXB = muXB+noise*randn(size(muXB));
muYB = muYB+noise*randn(size(muYB));

end

