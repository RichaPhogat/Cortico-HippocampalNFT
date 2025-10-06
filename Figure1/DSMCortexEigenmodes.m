clc;
clear all;
close all;

load('UncoupledCH.mat');
load('cortexMesh.mat');
load('corticalEigen200.mat');

fontSz = 25;

Fs = 1 / (dt * PointsToJump);

% Select eigenmode time range
StartingPoint = floor(length(modeCoeffPhiF(1,:)) / 4);
EndPoint = floor(length(modeCoeffPhiF(1,:))/1);
modeData = modeCoeffPhiF(:, StartingPoint:EndPoint);
nmodes = length(modeCoeffPhiF(:,1));
% CCC settings
windowLengthSec = 6;
overlapSec = 3;
maxLagSec = 1;

windowLength = round(Fs * windowLengthSec);
overlap = round(Fs * overlapSec);
stepSize = windowLength - overlap;
maxLagSamples = round(Fs * maxLagSec);

nModes = size(modeData, 1);
nTimePoints = size(modeData, 2);
nWindows = floor((nTimePoints - windowLength) / stepSize) + 1;

CorrCoefficientMatrix = zeros(nModes, nModes);

for jk = 1:nModes
    for kk = jk+1:nModes  % only upper triangle, symmetric matrix
        CCCs = zeros(1, nWindows);

        for w = 1:nWindows
            idxStart = (w - 1) * stepSize + 1;
            idxEnd = idxStart + windowLength - 1;

            x = modeData(jk, idxStart:idxEnd);
            y = modeData(kk, idxStart:idxEnd);

            % Remove mean (detrend option)
            x = x - mean(x);
            y = y - mean(y);

            % Compute cross-correlation with max lag
            [c, lags] = xcorr(x, y, maxLagSamples, 'coeff');

            % Save max absolute CCC
            CCCs(w) = max(abs(c));
            
        end

        % Average over all windows
        CorrCoefficientMatrix(jk, kk) = mean(CCCs);
    end
end

% Mirror to lower triangle
CorrCoefficientMatrix = CorrCoefficientMatrix + CorrCoefficientMatrix';

% Optional: mask diagonal or keep it zero
CorrCoefficientMatrix(1:nModes+1:end) = 0;

    mask = tril(ones(size(CorrCoefficientMatrix)), -1);
    CorrCoefficientMatrix = CorrCoefficientMatrix .* mask;

    figure;

    imagesc(CorrCoefficientMatrix);
    hColorbar = colorbar; 
    hColorbar.Orientation = 'horizontal'; 
    hColorbar.FontSize = fontSz; 
    % hColorbar.FontWeight = 'bold'; 
    hColorbar.Ticks = [0.1 0.8];  % or use [0.35 0.5] if preferred
    hColorbar.TickLabels = {'0.1', '0.8'};  % match the tick values
    clim([0.1 0.8]);
    axis square;
    colormap(bluewhiteredEigen(256,1));

    hColorbar.Position(1) = 0.27; % Adjust x position
    hColorbar.Position(2) = 0.55; % Adjust y position
    hColorbar.Position(3) = 0.5;  % Adjust width
    hColorbar.Position(4) = 0.03; % Adjust height

    view([-45 90]); 

    box off  % Remove full box

    % Draw custom left and bottom lines manually
    hold on;
    plot([1.5 1.5], [1.5 nmodes-0.5], 'k', 'LineWidth', 1.5);  % Left edge
    plot([1.5 nmodes-0.5], [nmodes-0.5 nmodes-0.5], 'k', 'LineWidth', 1.5);  % Bottom edge
    hold off;
    
    % Customize ticks
    ax = gca;
    ax.XTick = [2 102];
    ax.YTick = [2 102];
    ax.TickDir = 'out';  % Ticks pointing outward
    ax.LineWidth = 2;
    ax.FontSize = fontSz;
    % ax.FontWeight = 'bold';
    
    % Hide the top and right ticks
    ax.XRuler.Axle.Visible = 'off';
    ax.YRuler.Axle.Visible = 'off';


    xlim([1.5 109.5]);
    ylim([1.5 109.5]);
    
