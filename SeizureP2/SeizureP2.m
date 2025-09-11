clc
clear all
tic

    
rng(73924)
contn = 0;

intPath = 'integratorsNFT';
addpath(genpath(intPath));

% =========================================================================
%                      Time (recorded and simulation) 
% =========================================================================

totalTime = 185;                       % Total Time
tSpan = 0.1;                           % each epoch lasts tspan seconds
numEpochs = ceil(totalTime / tSpan);
StartRecording = 2;                  %The time to start recording data
StartRecording = ceil(StartRecording/tSpan);

dt = 10^-4;                             % Time Step
tcount = round(tSpan/dt);               % Time series length
t = linspace(0,tSpan,tcount+1);         % Real time values
options.InitialStep = dt;              % options.x values go to the solver
options.tstep = dt;
nsteps = tcount+1;                      % Total steps simulated

% =========================================================================
%                      System Parameters (units)
% =========================================================================

params = defineSystemParameters(totalTime, dt);
fields = fieldnames(params);
for i = 1:numel(fields)
    assignin('base', fields{i}, params.(fields{i}));
end


%==========================================================================
%         Initial conditions (rewritten in solver for 1st epoch)
%==========================================================================

[Y0, outputE1] = defineInitialConditionsE1(params,nsteps);
fields = fieldnames(outputE1);
for i = 1:numel(fields)
    assignin('base', fields{i}, outputE1.(fields{i}));
end
 
%%%%%%%%%%%%% Run the First epoch %%%%%%%%%%%%%%%%

[sol, paramsS] = sdeEMEigenmodeODE(@odefun1, @sdefun1, t, Y0, options, params, contn);

Y = sol.y(:,2:end);

L1 = length(Y(1,:));

%%%%%%%%%%%%% Diving into epochs %%%%%%%%%%%%%%%%

for jcontn = 2:numEpochs

nsteps = tcount+1;
t = linspace(0,tSpan,nsteps);

[Y0n, outputEn] = parseEpochState(Y, params);
fields = fieldnames(outputEn);

for i = 1:numel(fields)
    assignin('base', fields{i}, outputEn.(fields{i}));
end

[sol, paramsS] = sdeEMEigenmodeODEEpoch(@odefun1, @sdefun1, t, Y0n, options,...
    params, contn, Y);

Y = sol.y(:,2:end);   % Coefficient 1 is only used for storing the last time step data

phiModeCoeff = Y(1:numModesC,:);
etaModeCoeff = Y(2*numModesC+1:2*numModesC + numModesH,:);


if jcontn>=StartRecording

    contn = contn+1;
    PointsToJump = ceil(1/(1000*dt));
    destIdxStart = ((jcontn - StartRecording)*L1)/PointsToJump + 1;
    destIdxEnd = destIdxStart + (L1/PointsToJump - 1);
    srcIdx = 1:PointsToJump:L1;
    rcIdx2 = 1:PointsToJump:L1;

    modeCoeffPhiF(1:numModesC,destIdxStart:destIdxEnd) = phiModeCoeff(1:numModesC,srcIdx);
    modeCoeffEtaF(1:numModesH,destIdxStart:destIdxEnd) = etaModeCoeff(1:numModesH,srcIdx);
    Error(:,destIdxStart:destIdxEnd) = sol.Loss(:,rcIdx2);
    nuSeC(destIdxStart:destIdxEnd) = paramsS.nu_seC(srcIdx);
    nuSeH(destIdxStart:destIdxEnd) = paramsS.nu_seH(srcIdx);
    cc(destIdxStart:destIdxEnd) = paramsS.cc(srcIdx);
    ch(destIdxStart:destIdxEnd) = paramsS.ch(srcIdx);


% DataPlotFileName = '/scratch/hl36/rr8777//HippSeizureSlowP2I16.mat';
DataPlotFileName = 'SeizureP2.mat';

save(DataPlotFileName, 'modeCoeffPhiF', 'modeCoeffEtaF', 'dt', 'PointsToJump',...
    'params', 'paramsS', 'Error','nuSeC','nuSeH','cc','ch','-v7.3');
end
end

toc