function params = defineSystemParameters(totalTime, dt)
    %===============================================================
    %            Cortico-Hippocampal System Parameters
    %===============================================================

    thisDir = fileparts(mfilename('fullpath'));
    projectRoot = fullfile(thisDir);
    dataDir = fullfile(projectRoot, 'integratorsNFT');

    params.numModesC = 110;
    params.numModesH = 25;

    params.noiseType = 'eigenmode';
    
    % Load cortical mesh
    meshData = load(fullfile(dataDir, 'cortexMesh.mat'));
    params.cortexFaces    = meshData.cortexFaces;
    params.cortexVertices = meshData.cortexVertices;
    
    % Load hippocampal mesh
    hippoData = load(fullfile(dataDir, 'hippocampusMesh.mat'));
    params.hippocampusFaces    = hippoData.hippocampusFaces;
    params.hippocampusVertices = hippoData.hippocampusVertices;
    
    % Load cortical eigenmodes
    cortexEig = load(fullfile(dataDir, 'corticalEigen200.mat'));
    params.corticalEigenvalues  = cortexEig.corticalEigenvalues(1:params.numModesC);
    params.corticalEigenvectors = cortexEig.corticalEigenvectors(:,1:params.numModesC);

    % Load hippocampal eigenmodes
    hippoEig = load(fullfile(dataDir, 'hippoEigen200.mat'));
    params.hippoEigenvalues  = hippoEig.hippoEigenvalues(1:params.numModesH);
    params.hippoEigenvectors = hippoEig.hippoEigenvectors(:,1:params.numModesH);
    
    % Load couplings
    couplingCtoH = load(fullfile(dataDir, 'couplingCortexToHippo.mat'));
    couplingHtoC = load(fullfile(dataDir, 'couplingHippoToCortex.mat'));
    params.couplingCortexToHippo = couplingCtoH.couplingCortexToHippo;
    params.couplingHippoToCortex = couplingHtoC.couplingHippoToCortex;

    % Node counts
    params.nodesC = length(params.cortexVertices);  % Number of cortical nodes
    params.nodesH = length(params.hippocampusVertices);   % Number of hippocampal nodes
    
    % Noise Sources
    params.NoiseSourcesC = params.nodesC;
    params.NoiseSourcesH = params.nodesH;
    
    % Axonal conduction velocities (ms^-1)
    params.vC = 11.6;
    params.vH = 5.8;
    
    % Characteristic spatial ranges (m)
    params.reC = 0.1;
    params.reH = 0.1;
    
    % Maximum firing rates (s^-1)
    params.QmaxC = 340;
    params.QmaxH = 340;
    
    % Firing thresholds (V)
    params.thetaC = 12.92 / 1000; % Cortex
    params.thetaH = 12.92 / 100;  % Hippocampus
    
    % Threshold standard deviations (V)
    params.sigmaC = 3.4 / 1000;
    params.sigmaH = 3.05 / 100;
    
    % Synaptodendritic rise rates (s^-1)
    params.alphaC = 83.33;
    params.alphaH = 83.33;
    
    % Decay to rise ratios
    params.betabyalphaC = 9.23;
    params.betabyalphaH = 9.23;
    
    % Time delays (s)
    params.tauC = 0.5 * (80 / 1000);
    params.tauH = 0.5 * (10 / 1000);
    
    % Synaptic gains (V·s)
    params.nu_eeC = 3.03 / 1000;
    params.nu_eeH = 3.03 / 100;
    
    params.nu_eiC = -6 / 1000;
    params.nu_eiH = -6 / 100;
    
    params.nu_esC = 2.06 / 1000;
    params.nu_esH = 2.06 / 100;
    
    params.nu_seC = 1.4 / 1000;
    params.nu_seH = 1.5 / 100;
    
    params.nu_srC = -0.83 / 1000;
    params.nu_srH = -0.83 / 100;
    
    params.nu_snC = 0.98 / 1000;
    params.nu_snH = 0.98 / 100;
    
    params.nu_reC = 0.33 / 1000;
    params.nu_reH = 0.33 / 100;
    
    params.nu_rsC = 0.03 / 1000;
    params.nu_rsH = 0.03 / 100;
    
    params.phi_n0C = 1;
    
    %===============================================================
    %                     Pulse Train Parameters
    %===============================================================
    params.pulseDuration = 0.02;       % s
    params.pulseInterval = 0.1;        % s
    params.pulseStartTime = 8;        % s
    params.pulse_start_time_count = round(params.pulseStartTime / dt);
    
    % Compute number of pulses
    params.numPulses = ceil((totalTime - params.pulseStartTime) / params.pulseInterval) + 10;
    
    % External input amplitude and placeholder
    params.extInputAmplitude = 0;
    params.extInput = zeros(length(params.cortexVertices), round(1 / dt) + 1); % One epoch's time series length
    
end
