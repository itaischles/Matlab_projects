function [XA,YA,ZA,XB,YB,ZB,muXA,muYA,muZA,muXB,muYB,muZB] = create_crystal_3d(input_params)

tic

crystal_size = input_params{1};
Ax = input_params{2};
Ay=input_params{3};
Az=input_params{4};
dA=input_params{5};
dB=input_params{6};
muA=input_params{7};
muB=input_params{8};
angleA=input_params{9};
angleB=input_params{10};

xa = linspace(0, Ax*crystal_size(1), crystal_size(1)+1);
ya = linspace(0, Ay*crystal_size(2), crystal_size(2)+1);
za = linspace(0, Az*crystal_size(3), crystal_size(3)+1);
xb = linspace(0, Ax*crystal_size(1), crystal_size(1)+1);
yb = linspace(0, Ay*crystal_size(2), crystal_size(2)+1);
zb = linspace(0, Az*crystal_size(3), crystal_size(3)+1);

xa = xa-max(xa)/2+dA(1);
ya = ya-max(ya)/2+dA(2);
za = za-max(za)/2+dA(3);
xb = xb-max(xb)/2+dB(1);
yb = yb-max(yb)/2+dB(2);
zb = zb-max(zb)/2+dB(3);

[XA,YA,ZA] = meshgrid(xa,ya,za);
[XB,YB,ZB] = meshgrid(xb,yb,zb);

muxa = muA*sin(angleA(2))*cos(angleA(1))*ones(1,crystal_size(1)+1);
muya = muA*sin(angleA(2))*sin(angleA(1))*ones(1,crystal_size(2)+1);
muza = muA*cos(angleA(2))*ones(1,crystal_size(3)+1);
muxb = muB*sin(angleB(2))*cos(angleB(1))*ones(1,crystal_size(1)+1);
muyb = muB*sin(angleB(2))*sin(angleB(1))*ones(1,crystal_size(2)+1);
muzb = muB*cos(angleB(2))*ones(1,crystal_size(3)+1);


[muXA,muYA,muZA] = meshgrid(muxa,muya,muza);
[muXB,muYB,muZB] = meshgrid(muxb,muyb,muzb);

% add random noise
noise = 0.0;
muXA = muXA+noise*randn(size(muXA));
muYA = muYA+noise*randn(size(muYA));
muZA = muZA+noise*randn(size(muZA));
muXB = muXB+noise*randn(size(muXB));
muYB = muYB+noise*randn(size(muYB));
muZB = muZB+noise*randn(size(muZB));

elapsedTime = toc;
disp(strcat("Crystal initialization completed in ", num2str(elapsedTime), " seconds"));

end

