clc
clear all
close all

load('hippocampusMesh.mat');
load('hippoEigen200.mat');
load('UncoupledCH.mat');


plotColor = [1.0, 0.4, 0.3];

    Fs = 1 / (dt * PointsToJump);  

    numModesC = paramAll(1, end-2);
    numModesH = paramAll(1, end-1);
 
    FFTHippo = [];

    for jk = 2:numModesH
        yy = modeCoeffEtaF(jk, :);
        
        fftValues = calculatePowerSpectrum(yy(end/4:end),Fs);  
        FFTHippo(jk, :) = fftValues;  
    end

    avgFFT = mean(FFTHippo, 1);
    stdFFT = std(FFTHippo, 0, 1);

    nPoints = size(FFTHippo, 2);
    freqAxis = (0:nPoints-1) * (Fs / (2*nPoints));
    freqAxis = freqAxis(2:end);
    avgFFT = avgFFT(2:end);
    stdFFT = stdFFT(2:end);

lowerBound = max(0.00000001, avgFFT - stdFFT); 
upperBound = avgFFT + stdFFT;

fig = figure('Position', [100, 100, 1200, 420]);
plot(freqAxis, avgFFT, 'Color', plotColor, 'LineWidth', 1.5); 
hold on;
fill(([freqAxis, fliplr(freqAxis)]), ...
     ([lowerBound, fliplr(upperBound)]), ...
     plotColor, 'FaceAlpha', 0.5, 'EdgeColor', 'none'); 
hold off
% Plot the mean FFT in log-log scale
set(gca, 'YScale', 'log');

xlim([0.5 40])
ylim([0.0001 100])
xticks([5:10:50]);
yticks([0.0001, 10, 1000]);

set(gca, 'FontSize', 34);

function powerSpectrum = calculatePowerSpectrum(signal,Fs)
    n = length(signal);
    fftResult = fft(signal); % Normalize the FFT
    powerSpectrum = abs(fftResult).^2; % Square the magnitude for power
    powerSpectrum = powerSpectrum(1:floor(n/2)+1); % Keep only positive frequencies
    powerSpectrum(2:end-1) = 2 * powerSpectrum(2:end-1)/(n*Fs); % Adjust for single-sided spectrum
end
