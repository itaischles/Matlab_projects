clearvars

% see: Dazzi, 2010, Theory of infrared nanospectroscopy by photothermal
% induced resonance

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% cantilever
L = 450e-6; % length
H = 17e-6; % tip height
w = 50e-6; % width
thickness = 2.7e-6; % thickness
theta = 15 * pi/180; % angle relative to surface
rho_cantilever = 2330; % silicon density
E_cantilever = 110e9; % Young modulus
Rtip = 50e-9; % tip radius
numModes = 4; % number of eigenmodes to include in calculation

% sample
material = 'organic dye'; % see available options in "get_material_properties.m"
sphere_radius = 40e-9;

% cantilever + sample in contact
kx = 10; % lateral tip+sample spring constant (high = high friction, low = low friction)
deltax = 10e-6; % position from cantilever end where force is applied
Gamma = 1/(50e-6); % oscillations damping

% pulse train
Pavg = 100e-6; % average power in beam
reprate = 98e3; % repetition rate
spotradius = 2e-6; % focused beam radius

% AFM detection system
foclen = 40e-3; % focal length of lens which focuses red LD onto cantilever
WtoI = 0.45; % watt to amps conversion of 4QPD
ItoV = 1e4 * 100; % amps to volts of electronics card (transimpedance resistor * voltage gain on ac output)
PLD = 3e-3; % power of laser diode used for cantilever movement detection
LDbeamsize = 3e-3; % beam size on 4QPD
noise_dens = 10e-6; % V/sqrt(Hz)
LIBW = 300; % lock in bandwidth (Hz)

% time
tmax = 2e-3; % maximum simulation time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% CALCULATIED PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = linspace(0, L, 1000); % position along cantilever
t = linspace(0, tmax, 20000); % time
f = linspace(0, 0.5./(t(2)-t(1)), numel(t)/2+1);
material_properties = get_material_properties(material);
alphaT = material_properties(1); % linear thermal expansion coefficient of sample (m/(m*degC))
C = material_properties(2); % specific heat capacity of sample (J/(kg*degC))
rho = material_properties(3); % mass density of sample (kg/m^3)
kappa = material_properties(4); % thermal conductivity of sample (W/(m*degC))
E_sample = material_properties(5); % Young's modulus (Pa)
nreal = real(material_properties(6)); % real part of refractive index at 520 nm
nimag = imag(material_properties(6)); % imaginary part of refractive index at 520 nm
E = 1/(1/E_cantilever + 1/E_sample); % combined Young's modulus (sample+cantilever)
tau_relax = rho*C/(3*kappa)*sphere_radius^2; % thermal relaxation time
crossSectionArea = w*thickness;
I_cant = w*thickness^3 / 12; % cantilever moment of inertia
Vol = 4/3*pi*sphere_radius^3; % sample volume
kc = E_cantilever*I_cant*3/L^3; % free cantilever spring constant
kz = 16/15*E*sqrt(Rtip); % perpendicular tip+sample spring constant
% % % kz = 1e-8; % perpendicular tip+sample spring constant for free cantilever (should be =0)
U = kc*L^2/(3*kx*H^2); % see Eq. A1 in paper
V = kc/(3*kz); % see Eq. A1 in paper
Epulse = Pavg * (1/reprate); % energy in single pulse 
fluence = Epulse/(pi*spotradius^2); % pulse fluence
Eabssquared = 2*fluence/(3e8*8.85e-12); % magnitude squared of electric field of incident beam times pulse duration
Eabs = 2*pi/520e-9 * 3e8*8.85e-12 * 9*nimag*nreal/(nreal^2+2)^2 * Eabssquared * Vol; % absorbed energy in sample after each pulse
Tmax = Eabs/(rho*C*Vol); % maximum temperature change in sample after pulse turns off
noise = sqrt(noise_dens^2*LIBW); % integrated noise

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% CALCULATE TEMPERATURE VS. TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_reprate = 1/reprate;
dt = t(2)-t(1);
T = zeros(1,numel(t));
numpulses = floor(tmax/t_reprate);
for i=1:numpulses
    shifted_T_pulse = Tmax*circshift(exp(-t/tau_relax), (i-1)*round(t_reprate/dt));
    shifted_T_pulse(1:(i-1)*round(t_reprate/dt)) = 0;
    T = T + shifted_T_pulse;
end
avgT = ones(1, numel(t))*mean(T); % average temperature increase
sphere_expansion = sphere_radius*alphaT*T/3; % sample thermal expansion
avg_expansion = ones(1, numel(t))*mean(sphere_expansion); % average thermal expansion

% % % figure('Position', [78 355 438 286])
% % % yyaxis right
% % % ylim([0, max(sphere_expansion)]);
% % % ylabel('Thermal expansion (m)')
% % % title(material)
% % % yyaxis left
% % % hold on
% % % plot(t*1e6, T)
% % % plot(t*1e6, avgT, 'r--')
% % % xlabel('Time (\mus)')
% % % ylabel('Temperature increase (K)')
% % % ylim([0, max(T*1.1)]);
% % % xlim([0,200])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% CALCULATE SOLUTIONS TO EIGENVALUE EQ.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eigEq = @(x) -1 + cos(x).*cosh(x) - V*x.^3.*(cos(x).*sinh(x)+sin(x).*cosh(x)) - U*x.*(sin(x).*cosh(x)-cos(x).*sinh(x)) - U*V*x.^4.*(1+cos(x).*cosh(x));

X = linspace(0,15,5000);
mechFreq = sqrt(E_cantilever*I_cant/(rho_cantilever*crossSectionArea)).*(X/L).^2;
Y = smooth(abs(gradient(log(abs(eigEq(X)))))); % make it easier to find solutions of eigenvalue equation
[~, inds] = findpeaks(Y(Y<max(Y)), 'MinPeakProminence', 0.3);
Xn = X(inds(1:numModes));
clear inds

% % % plot(X, Y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% CALCULATE BENDING EIGENMODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modes = zeros(numel(x), numModes);
labels = cell(1, numModes);
beta = zeros(numModes, 1);
for n=1:numModes
    beta(n) = Xn(n)/L;
    modes(:,n) = (cos(beta(n)*x)-cosh(beta(n)*x)) - ((cos(beta(n)*L)-cosh(beta(n)*L))./(sin(beta(n)*L)-sinh(beta(n)*L))).*(sin(beta(n)*x)-sinh(beta(n)*x));
end

% % % figure; hold on
% % % for n=1:numModes
% % %     plot(x/L, modes(:,n))
% % %     labels{n} = strcat('Mode ', num2str(n));
% % % end
% % % hold off
% % % xlabel('x/L')
% % % ylabel('Normalized mode amplitude (a.u.)')
% % % legend(labels, 'location', 'best')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% CANTILEVER BENDING DUE TO THERMAL STRESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltax_ind = round((L-deltax)/L * numel(x)); % index of position where LD shines on cantilever
B = (cos(theta)+(H/deltax)*sin(theta))*sphere_radius*kz*alphaT/3; % see Eq. 11
Pn = zeros(numModes, 1);
wn = sqrt(E_cantilever*I_cant/(rho_cantilever*crossSectionArea)).*beta.^2;
zn = zeros(numel(x), numel(t), numModes); % respose of nth eigenmode
z = zeros(numel(x), numel(t));
Sn = zeros(numel(t), numModes); % signal on 4QPD - contribution from mode n
S = zeros(numel(t), 1); % total signal on 4QPD
for n=1:numModes
    dgdx = gradient(modes(:,n));
    Pn(n) = B*deltax/(rho*crossSectionArea*L) * dgdx(deltax_ind);
    temporalPart = conv(sin(wn(n)*t).*exp(-Gamma/2 * t), T, 'full');
    temporalPart = temporalPart(1:numel(t));
    zn(:, :, n) = Pn(n)/wn(n)*modes(:,n) .* temporalPart;
    z = z + zn(:,:,n);
    Sn(:, n) = Pn(n)/wn(n)*tan(2*atan(dgdx(deltax_ind))) .* temporalPart * foclen * 2*WtoI*ItoV*PLD/LDbeamsize;
    S = S + Sn(:, n);
end
S_noise = S + noise*randn(size(S));
S_FFT = fft(S-mean(S))/numel(S);
S_FFT = S_FFT(1:numel(S_FFT)/2+1);
S_FFT(2:end-1) = 2*S_FFT(2:end-1);
S_FFT_noise = fft(S_noise-mean(S_noise))/numel(S_noise);
S_FFT_noise = S_FFT_noise(1:numel(S_FFT_noise)/2+1);
S_FFT_noise(2:end-1) = 2*S_FFT_noise(2:end-1);

figure
hold on
% % % % % % plot(t*1e6, zn(deltax_ind, :, 1))
% % % % % % plot(t*1e6, zn(deltax_ind, :, 2))
% % % % % % plot(t*1e6, zn(deltax_ind, :, 3))
% % % % % % plot(t*1e6, zn(deltax_ind, :, 4))
% % % plot(t*1e6, z(deltax_ind, :))
% % % plot(t*1e6, T/max(T)*max(abs(z(deltax_ind, :))))
% % % xlabel('Time (\mus)')
plot(t*1e6, S)
plot(t*1e6, S_noise)
xlabel('Time (\mus)')

% % % figure('Position', [764 348 560 288])
% % % hold on
% % % % % % plot(t*1e6, Sn(:, 1))
% % % % % % plot(t*1e6, Sn(:, 2))
% % % % % % plot(t*1e6, Sn(:, 3))
% % % % % % plot(t*1e6, Sn(:, 4))
% % % plot(t(t<600e-6)*1e6, S(t<600e-6), 'Color', [0.3,0.3,0.3])
% % % % % % plot(t*1e6, T/max(T)*max(abs(S)))
% % % xlabel('Time (\mus)')
% % % ylabel('PTE signal (V)')

figure('Position', [382 504 560 288])
plot(f, abs(S_FFT_noise), 'Color', [0.3,0.3,0.3])
plot(f, abs(S_FFT_noise))
xlim([0,1e6])
xlabel('Frequency (Hz)')
ylabel('PTE signal (V)')

%%%%%%%%%%%%%%%%%%%
% MOVIE
%%%%%%%%%%%%%%%%%%%
% % % WriterObj = VideoWriter('movie.mp4', 'MPEG-4');
% % % WriterObj.FrameRate = 30;
% % % open(WriterObj);
% % % figure('Position', [403 485 560 181], 'Color', 'white');
% % % ax = axes;
% % % xlabel('x (\mum)')
% % % ylabel('amplitude (\mum)')
% % % for i=1:1:round(numel(t)*(1/50))
% % %     hold on
% % %     plot(x*1e6, zn(:,i,1)*1e6, 'LineWidth', 5, 'Color', [255, 153, 153]/255);
% % %     plot(x*1e6, zn(:,i,2)*1e6, 'LineWidth', 5, 'Color', [153, 255, 153]/255);
% % %     plot(x*1e6, zn(:,i,3)*1e6, 'LineWidth', 5, 'Color', [153, 153, 255]/255);
% % %     plot(x*1e6, zn(:,i,4)*1e6, 'LineWidth', 5, 'Color', [255, 255, 153]/255);
% % %     plot(x*1e6, z(:,i)*1e6, 'LineWidth', 5, 'Color', 'k');
% % %     plot(x(deltax_ind)*1e6, z(deltax_ind,i)*1e6, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
% % %     ylim([-H*1e6, H*1e6])
% % %     title(strcat('t=', num2str(round(1e6*(t(i)))), '\mus'))
% % %     frame = getframe(gcf);
% % %     writeVideo(WriterObj,frame);
% % %     cla
% % % end
% % % close(WriterObj);