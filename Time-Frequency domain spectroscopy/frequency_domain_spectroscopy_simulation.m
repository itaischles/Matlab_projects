clearvars

% main pulse
lambda0 = 550e-9; % central wavelength [m]
f0 = 3e8/lambda0; % central frequency of main pulse [Hz]
num_cycles_in_pulse = 6; % (approximately) the number of optical cycles in the main pulse (must be > 1 and preferably < 40)
dt = 1/(15*f0); % time resolution of time vector
tau_pulse = (num_cycles_in_pulse/3)/f0; % temporal width of main pulse (assumed Gaussian)
t = linspace(-300*tau_pulse, 300*tau_pulse, round(300*tau_pulse/dt)); % time vector for pulse
E0 = @(dt) exp(-0.5*(t-dt).^2/tau_pulse^2); % pulse envelope function
Emain = E0(0).*cos(2*pi*f0*t); % main pulse electric field
main_pulse_spectrum = fft(Emain); % the spectrum of the main pulse
main_pulse_spectrum = main_pulse_spectrum/max(abs(main_pulse_spectrum)); % normalized pulse spectrum

% molecular response pulse
tau_relax = 15*tau_pulse; % relaxation time of molecular vibration
t_delay = 50*tau_pulse;
lambda_mol1 = 580e-9; % wavelength of molecular coherence 1
lambda_mol2 = 530e-9; % wavelength of molecular coherence 2
f_mol1 = 3e8/lambda_mol1; % frequency of molecular coherence 1
f_mol2 = 3e8/lambda_mol2; % frequency of molecular coherence 1
Hmol = exp(-t/tau_relax).*sin(2*pi*f_mol1*t) + exp(-t/tau_relax).*sin(2*pi*f_mol2*t); % molecular response function
Hmol(t<0) = 0; % making the molecular response function causal
Emol = conv(Emain, Hmol, 'same'); % the convolution between the molecular response function and the pulse is proportional to the output signal field
Emol = Emol/max(Emol); % normalized signal field

% reference pulse (local oscillator)
Eref = E0(t_delay).*cos(2*pi*f0*(t-t_delay)); % delayed reference pulse (delayed using movable mirror)

% frequency domain (what the spectrometer sees)
Edet = Emol+Eref;
det_spectrum = abs(fft(Edet)).^2; % the spectrum that the spectrometer sees
det_spectrum = det_spectrum/max(det_spectrum); % normalized detector spectrum
Fmax = 1/(2*dt);
freq = linspace(0,Fmax/2,round(numel(t)/2));
hold on
plot1 = plot(freq, abs(main_pulse_spectrum(1:round(numel(t)/2))), 'LineWidth', 0.5, 'Color', 'k', 'LineStyle', '--');
plot2 = plot(freq, abs(det_spectrum(1:round(numel(t)/2))), 'LineWidth', 0.5, 'Color', 'r', 'LineStyle', '-');

