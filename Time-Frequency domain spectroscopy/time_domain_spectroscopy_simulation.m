% parameters that reproduce experimental WL probe:
%   lambda0 = 595e-9
%   tau_pulse = 100e-15
%   a = 22

clearvars

% refresh plots every X steps for speeding up simulation and display
refreshEveryXSteps = 50;

% main pulse
lambda0 = 595e-9; % central wavelength [m]
f0 = 3e8/lambda0; % central frequency of main pulse [Hz]
dt = 1/(15*f0); % time resolution of time vector
tau_pulse = 100e-15; % temporal width of main pulse (sec)
a = 22; % unitless number. When 0: no chirp
a2 = 80;
t = linspace(-100*tau_pulse, 100*tau_pulse, round(100*tau_pulse/dt)); % time vector for pulse
Emain = create_pulse(t, 0, lambda0, a, tau_pulse); % create main pulse centered at time 0
pulse_spectrum = abs(fft(Emain)).^2; % the spectrum of the main pulse
pulse_spectrum = pulse_spectrum/max(pulse_spectrum); % normalized pulse spectrum

% molecular response pulse
tau_relax = 0.2*tau_pulse; % relaxation time of molecular transition (vibration/electronic...)
t_delay_min = -3*tau_pulse;
t_delay_max = 20*tau_relax;
t_delay = linspace(t_delay_min, t_delay_max, round((t_delay_max-t_delay_min)/(dt))); % delay-time vector
lambda_mol1 = 575e-9; % wavelength of molecular coherence 1
lambda_mol2 = 500e-9; % wavelength of molecular coherence 2
f_mol1 = 3e8/lambda_mol1; % frequency of molecular coherence 1
f_mol2 = 3e8/lambda_mol2; % frequency of molecular coherence 1
Hmol = -exp(-t/tau_relax).*sin(2*pi*f_mol1*t)-exp(-t/tau_relax).*sin(2*pi*f_mol2*t); % molecular response function
% % % Hmol = -exp(-t/tau_relax).*sin(2*pi*f_mol1*t); % molecular response function
Hmol(t<0) = 0; % making the molecular response function causal
Emol = conv(Emain, Hmol, 'same'); % the convolution between the molecular response function and the pulse is proportional to the output signal field
Emol = 0.0*Emol/max(Emol); % signal field

% frequency vector
Fmax = 1/abs(t_delay(2)-t_delay(1));
freq = linspace(0,Fmax/2,round(numel(t_delay)/2));

% wavelength vector for plotting (zoomed)
f0LowFrac = 0.6;
f0HighFrac = 1.5;
[~, ind1] = min(abs(freq-f0LowFrac*f0));
[~, ind2] = min(abs(freq-f0HighFrac*f0));
wavelength_zoom = linspace(f0LowFrac*3e8/f0*1e9, f0HighFrac*3e8/f0*1e9, ind2-ind1+1);

% initialize figure and axes
fig = figure;
fig.Position = [416 87 600 800];
ax1 = subplot(3,1,1);
ax2 = subplot(3,1,2);
ax3 = subplot(3,1,3);
plotvec = round(numel(t)/8):round(numel(t)*7/8)-1; % index vector for plotting purposes

% initialize data for first plots
Eref = create_pulse(t, t_delay(1), lambda0, a, tau_pulse); % delayed reference pulse (delayed using movable mirror)
Eref2 = create_pulse(t, t_delay(1), lambda0, a2, tau_pulse); % delayed reference pulse (delayed using movable mirror)
Intensity = abs(Emain+Eref+Eref2).^2; % the total intensity on the detector vs. time
I_det_sampled_first = trapz(t, Intensity); % the integrated detector intensity (since the detector is slow it basically integrates over the intensity of the pulses)
I_det_sampled = ones(1,numel(t_delay))*I_det_sampled_first; % initializing the detector intensity vector (as a function of delay-time)
ErefPlotData = real(Eref(plotvec)); % Y source variable for plotting reference pulse
Eref2PlotData = real(Eref2(plotvec)); % Y source variable for plotting reference pulse
IntensityPlotData = Intensity(plotvec); % Y source variable for plotting intensity
spectrum = fft(I_det_sampled);
[~, peak_ind] = max(abs(spectrum(ind1:ind2)));
abs_spectrum = abs(spectrum)/abs(spectrum(peak_ind+(ind1-1))); % normalize spectrum such that the peak has height 1
Fmax = 1/(2*dt);
freq_pulse_spectrum = linspace(0,Fmax/2,round(numel(t)/2));
measSpectPlotData = abs_spectrum(1:round(numel(t_delay)/2)); % Y source variable for plotting measured spectrum

% initialize movie 
WriterObj = VideoWriter('movie.mp4', 'MPEG-4');
WriterObj.FrameRate = 30;
open(WriterObj);

% plot on ax1
hold (ax1, 'on')
mainPulsePlot = plot(ax1, t(plotvec)*1e15, real(Emain(plotvec)), 'LineWidth', 0.1, 'Color', 'b', 'LineStyle', '-');
refPulsePlot = plot(ax1, t(plotvec)*1e15, ErefPlotData, 'LineWidth', 0.1, 'Color', 'r', 'LineStyle', '-');
refPulse2Plot = plot(ax1, t(plotvec)*1e15, Eref2PlotData, 'LineWidth', 0.1, 'Color', 'r', 'LineStyle', '-');
intensityPlot = plot(ax1, t(plotvec)*1e15, IntensityPlotData, 'LineWidth', 0.1, 'Color', 'k', 'LineStyle', '-');
molecRespPlot = plot(ax1, t(plotvec)*1e15, real(Emol(plotvec)), 'LineWidth', 0.1, 'Color', 'm', 'LineStyle', '-');
hold (ax1, 'off')
xlim(ax1, [min(t_delay*1e15), max(t_delay*1e15)]);
ylim(ax1, [-2, 2]);
xlabel(ax1, 'time [fs]');
title(ax1, 'impinging pulses')
legend(ax1, 'Main pulse','Ref. pulse','Intensity','Molecular response');
% % % legend(ax1, 'Main pulse','Intensity','Ref. pulse');
refPulsePlot.YDataSource = 'ErefPlotData'; % link Y data of line to source
refPulse2Plot.YDataSource = 'Eref2PlotData'; % link Y data of line to source
intensityPlot.YDataSource = 'IntensityPlotData'; % link Y data of line to source

% plot on ax2
integratedIntensityPlot = plot(ax2, t_delay*3e8*1e6, I_det_sampled);
xlabel(ax2, 'mirror displacement [um]');
title(ax2, 'detector output (integrated intensity)')
xlim(ax2, [min(t_delay*3e8*1e6), max(t_delay*3e8*1e6)]);
ylim(ax2, [0,5e-13]);
integratedIntensityPlot.YDataSource = 'I_det_sampled'; % link Y data of line to source

% plot on ax3
hold (ax3, 'on')
measSpectPlot = plot(ax3, 3e8./freq*1e9, abs(spectrum(1:round(numel(t_delay)/2))), 'LineWidth', 1, 'Color', 'k', 'LineStyle', '-');
probeSpectPlot = plot(ax3, 3e8./freq_pulse_spectrum*1e9, pulse_spectrum(1:round(numel(t)/2)), 'LineWidth', 0.5, 'Color', 'r', 'LineStyle', '--');
xlim(ax3, [min(wavelength_zoom), max(wavelength_zoom)]);
ylim(ax3, [0, 1.2]);
xlabel(ax3, 'wavelength [nm]')
title(ax3, 'measured spectrum & probe spectrum')
measSpectPlot.YDataSource = 'measSpectPlotData'; % link Y data of line to source

% plot refreshing loop
for i=1:numel(t_delay)
    
    % refresh data
    Eref = create_pulse(t, t_delay(i), lambda0, a, tau_pulse); % delayed reference pulse (delayed using movable mirror)
    Eref2 = create_pulse(t, t_delay(i), lambda0, a2, tau_pulse); % delayed reference pulse (delayed using movable mirror)
    Intensity = abs(Emain+Eref+Eref2).^2; % the total intensity on the detector vs. time
    I_det_sampled(i) = trapz(t, Intensity); % the integrated detector intensity (since the detector is slow it basically integrates over the intensity of the pulses)
    spectrum = fft(I_det_sampled);
    [~, peak_ind] = max(abs(spectrum(ind1:ind2)));
    abs_spectrum = abs(spectrum)/abs(spectrum(peak_ind+(ind1-1))); % normalize spectrum such that the peak has height 1
    measSpectPlotData = abs_spectrum(1:round(numel(t_delay)/2)); % Y source variable for plotting measured spectrum
    
    % do every 'refreshEveryXSteps' frames
    if mod(i, refreshEveryXSteps)==0
        
        % refresh plots on ax1
        ErefPlotData = real(Eref(plotvec)); % refresh Y source variable for plotting reference pulse
        Eref2PlotData = real(Eref2(plotvec)); % refresh Y source variable for plotting reference pulse
        IntensityPlotData = Intensity(plotvec); % refresh Y source variable for plotting intensity
        refreshdata(refPulsePlot)
        refreshdata(refPulse2Plot)
        refreshdata(intensityPlot)

        % refresh plot on ax2
        refreshdata(integratedIntensityPlot)

        % plot fourier transform of short pulse and the intensity on the
        % detector
        refreshdata(measSpectPlot)

        % get frame for movie
% % %         frame = getframe(fig, [0,0,fig.Position(3),fig.Position(4)]);
        frame = getframe(fig);
        writeVideo(WriterObj,frame);
    end
end

close(WriterObj);
disp('done!');