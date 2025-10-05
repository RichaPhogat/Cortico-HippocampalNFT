clc;
clear all;
close all;

load('UncoupledCH.mat')
load('hippocampusMesh.mat');
load('hippoEigen200.mat');
fontsz = 27;
folderPath = 'D:\Two Sheets\NCIDataNFT\SeizureSimulations\dataCheckNCI';

fileNames = {'RampOldCheck1F.mat'};

colors = {[1.0, 0.45, 0.65], [0.1, 0.4, 0.4]};


    Fs = 1 / (dt * PointsToJump);

    numModesC = paramAll(1, end-1);
    numModesH = paramAll(1, end);

    numPoints = size(modeCoeffPhiF, 2);
    t = (0:numPoints-1) * dt * PointsToJump;

    mode2 = modeCoeffEtaF(2, :);
    mode3 = modeCoeffEtaF(3, :);
    mode25 = modeCoeffEtaF(25, :);

    windowSize = 10; 
    runningAverage2 = movmean(mode2, windowSize);
    runningAverage110 = movmean(mode25, windowSize);


    detrendedTimeSeries2 = mode2;
    detrendedTimeSeries25 = mode25;

    fig = figure('Position', [100, 100, 1200, 420]);
    hold on;
    plot(t, detrendedTimeSeries2, 'Color', colors{1}, 'LineWidth', 3,...
        'DisplayName', '$\zeta^h_{2}$');
    plot(t, detrendedTimeSeries25, 'Color', colors{2}, 'LineWidth', 3,...
        'DisplayName', '$\zeta^h_{25}$');
    hold off;

    xlim([53.5 57.5]);
    xticks([54, 55, 56, 57]);
    yticks([-10, 0, 10]);
    set(gca, 'FontSize', fontsz);
    % xlabel('Time (s)');
    % ylabel('Mode Coeff. (\zeta_{nc})');
    
    legend('show', 'Interpreter', 'latex', 'FontSize', fontsz, 'FontWeight',...
        'bold','Location', 'best', 'box','off');
