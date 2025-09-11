function [Y0, outputE1] = defineInitialConditionsE1(params,nsteps)
    %===============================================================
    %       Generate Initial Conditions for First Epoch
    %===============================================================
    % params : Struct from defineSystemParameters
    % Y0     : Initial state vector for solver

    numModesC = size(params.corticalEigenvectors, 2);
    numModesH = size(params.hippoEigenvectors, 2);
    nodesC = params.nodesC;
    nodesH = params.nodesH;

    % Start from baseline (or small perturbation if needed)
    randScale = 0.00;  % Set >0 for random noise
    outputE1.phiModeCoeff  = randScale * randn(numModesC, nsteps)+1;
    outputE1.phiModeCoeff2 = zeros(numModesC, nsteps);
    outputE1.etaModeCoeff  = randScale * randn(numModesH, nsteps)+1;
    outputE1.etaModeCoeff2 = zeros(numModesH, nsteps);

    outputE1.VEC  = randScale * randn(nodesC, nsteps);
    outputE1.VEH  = randScale * randn(nodesH, nsteps);
    outputE1.VRC  = randScale * randn(nodesC, nsteps);
    outputE1.VRH  = randScale * randn(nodesH, nsteps);
    outputE1.VSC  = randScale * randn(nodesC, nsteps);
    outputE1.VSH  = randScale * randn(nodesH, nsteps);

    outputE1.VEC2 = zeros(nodesC, nsteps);
    outputE1.VEH2 = zeros(nodesH, nsteps);
    outputE1.VRC2 = zeros(nodesC, nsteps);
    outputE1.VRH2 = zeros(nodesH, nsteps);
    outputE1.VSC2 = zeros(nodesC, nsteps);
    outputE1.VSH2 = zeros(nodesH, nsteps);

    % Pack all initial states into a single column vector
    Y0 = [
        outputE1.phiModeCoeff(:,1);
        outputE1.phiModeCoeff2(:,1);
        outputE1.etaModeCoeff(:,1);
        outputE1.etaModeCoeff2(:,1);
        outputE1.VEH(:,1);
        outputE1.VEH2(:,1);
        outputE1.VRH(:,1);
        outputE1.VRH2(:,1);
        outputE1.VSH(:,1);
        outputE1.VSH2(:,1);
        outputE1.VEC(:,1);
        outputE1.VEC2(:,1);
        outputE1.VRC(:,1);
        outputE1.VRC2(:,1);
        outputE1.VSC(:,1);
        outputE1.VSC2(:,1)
    ];
end
