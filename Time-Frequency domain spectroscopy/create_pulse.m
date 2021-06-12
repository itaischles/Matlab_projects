function [Epulse] = create_pulse(t, t0, lambda0, a, tau)
% input parameters:
% t: time vector
% t0: time displacement
% lambda0: central wavelength (m)
% a: chirp parameter (unitless)
% tau: pulse temporal width

f0 = 3e8/lambda0; % central frequency of pulse (Hz)
dt = 1/(20*f0); % time resolution of time vector (sec) (for plotting only)

f = linspace(0, 1/dt, numel(t)); % frequency vector for plotting pulse

Gaussian = @(t0, a, tau) exp(-((1-1i*a)/tau^2)*(t-t0).^2); % generic complex Gaussian (time = vector;  a,tau = scalar parameters) (see Eq. 22.1-11 in Saleh & Teich)
expPhase = @(t0) exp(-1i*2*pi*f0*(t-t0)); % fast component of wave (depends on central wavelength)

Epulse = real(Gaussian(t0, a, tau).*expPhase(t0)); % pulse electric field


%%%%%%%%%%%%%%%%%
% pulse spectrum
%%%%%%%%%%%%%%%%%

Epulse_fft = fft(Epulse); % the FT of the pulse
Epulse_fft = Epulse_fft/max(Epulse_fft); % normalized pulse FT


%%%%%%%%%%%
% plotting
%%%%%%%%%%%

% % % figure;
% % % % plot time domain
% % % subplot(2,1,1);
% % % hold on
% % % plot(t*1e15, real(Epulse), 'LineWidth', 1, 'Color', 'b');
% % % plot(t*1e15, abs(Gaussian(t0, a, tau)).^2, 'LineWidth', 2, 'Color', 'r');
% % % % % % plot(t*1e15, a*((t-t0)/tau).^2, 'LineWidth', 2, 'Color', 'm');
% % % xlim([-3*tau*1e15,3*tau*1e15]);
% % % ylim([-1,1]);
% % % xlabel('time (fs)');
% % % 
% % % % % % % plot frequency domain
% % % % % % subplot(2,1,2);
% % % % % % plot(f, real(Epulse_fft));
% % % % % % xlim([0, 1/(2*dt)]);
% % % 
% % % % plot spectrum (wavelength)
% % % subplot(2,1,2);
% % % plot(3e8./f * 1e9, abs(Epulse_fft));
% % % xlabel('wavelength (nm)');
% % % 
% % % end
% % % 
