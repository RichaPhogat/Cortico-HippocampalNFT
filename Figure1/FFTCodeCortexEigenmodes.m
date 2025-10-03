clc
clear all
close all

load('UncoupledCH.mat');
load('cortexMesh.mat');
load('hippocampusMesh.mat');
load('corticalEigen200.mat');
load('hippoEigen200.mat');

plotColor = [0.5, 0.6, 1.0];


    Fs = 1 / (dt * PointsToJump);  

    numModesC = paramAll(1, end-2);
    numModesH = paramAll(1, end-1);

   
    FFTCortex = [];

    for jk = 2:numModesC
        yy = modeCoeffPhiF(jk, :);
        
        fftValues = calculatePowerSpectrum(yy(end/4:end),Fs);  
        FFTCortex(jk, :) = fftValues;  
    end

    avgFFT = mean(FFTCortex, 1);
    stdFFT = std(FFTCortex, 0, 1);
   

    nPoints = size(FFTCortex, 2);
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
set(gca, 'YScale', 'log');

xlim([0.5 40])
ylim([0.01 80])
xticks([5:10:50]);
yticks([0.01, 10]);

set(gca, 'FontSize', 34);

function powerSpectrum = calculatePowerSpectrum(signal,Fs)
    n = length(signal);
    fftResult = fft(signal); % Normalize the FFT
    powerSpectrum = abs(fftResult).^2; % Square the magnitude for power
    powerSpectrum = powerSpectrum(1:floor(n/2)+1); % Keep only positive frequencies
    powerSpectrum(2:end-1) = 2 * powerSpectrum(2:end-1)/(n*Fs); % Adjust for single-sided spectrum
end
