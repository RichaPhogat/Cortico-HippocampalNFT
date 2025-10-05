clc
clear all

load('UncoupledCH.mat');
load('hippoEigen200.mat');
load('hippocampusMesh.mat');

fs = 1 / (dt * PointsToJump);

for i = 1:1:25
    [psdHippoi] = calculatePowerSpectrum(modeCoeffEtaF(i,end/4:end),fs);
    PSDHippo(i,:) = smooth(psdHippoi,5);
end
    nPoints = size(psdHippoi, 2);
    freqAxisHippo = (0:nPoints-1) * (fs / (2*nPoints)); 

fig = figure('Position', [100, 100, 600, 420]);

smoothedPSDHippo = (PSDHippo);
imagesc(2:25, freqAxisHippo, smoothedPSDHippo');
ylim([1 15]);
xlim([1.5 25.5]);
xticks([5, 15, 25]); 
xticklabels({'5', '15', '25'});
yticks([10, 25]); 
yticklabels({'10', '25'});

colormap(bluewhiteredEigen(256,0));
c = colorbar;
c.Ticks = [0 10]; 
c.TickLabels = {'0', '10'}; 
clim([0 15])
set(gca, 'FontSize', 28, 'LineWidth', 1, 'Layer', 'top');
axis xy;

function powerSpectrum = calculatePowerSpectrum(signal,Fs)
    n = length(signal);
    fftResult = fft(signal); 
    powerSpectrum = abs(fftResult).^2; % Square the magnitude for power
    powerSpectrum = powerSpectrum(1:floor(n/2)+1); % Keep only positive frequencies
    powerSpectrum(2:end-1) = 2 * powerSpectrum(2:end-1)/(n*Fs); % Adjust for single-sided spectrum
end
