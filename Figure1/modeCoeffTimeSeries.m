clc;
clear all;
close all;

% Load data
load('UncoupledCH.mat');
load('cortexMesh.mat');
load('hippocampusMesh.mat');
load('corticalEigen200.mat');
load('hippoEigen200.mat');
fontSz = 27;


colors = {[0.45, 0.65, 1.0], [0.1, 0.4, 0.4]};



    Fs = 1 / (dt * PointsToJump);

    numModesC = paramAll(1, end-1);
    numModesH = paramAll(1, end);

    numPoints = size(modeCoeffPhiF, 2);
    t = (0:numPoints-1) * dt * PointsToJump;

    mode2 = modeCoeffPhiF(2, :);
    mode3 = modeCoeffPhiF(3, :);
    mode110 = modeCoeffPhiF(110, :);

    windowSize = 10; 
    runningAverage2 = movmean(mode2, windowSize);
    runningAverage110 = movmean(mode110, windowSize);


    detrendedTimeSeries2 = mode2;
    detrendedTimeSeries110 = mode110;

    fig = figure('Position', [100, 100, 1200, 420]);
    hold on;
    plot(t, detrendedTimeSeries2, 'Color', colors{1}, 'LineWidth', 3,...
        'DisplayName', '$\zeta^c_{2}$');
    plot(t, detrendedTimeSeries110, 'Color', colors{2}, 'LineWidth', 3,...
        'DisplayName', '$\zeta^c_{110}$');
    hold off;

    xlim([53.5 57.5]);
    xticks([54, 55, 56, 57]);
    yticks([-20, 0, 20]);
    set(gca, 'FontSize', fontSz);
    % xlabel('Time (s)');
    % ylabel('\zeta_{nc}');
    
    legend('show', 'Interpreter', 'latex', 'FontSize', fontSz,...
        'FontWeight', 'bold', 'Location', 'best', 'box','off');

