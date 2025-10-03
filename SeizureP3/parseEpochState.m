function [Y0, outputs] = parseEpochState(Y, params)

    numModesC = size(params.corticalEigenvectors, 2);
    numModesH = size(params.hippoEigenvectors, 2);
    nodesC = params.nodesC;
    nodesH = params.nodesH;

    outputs.phiModeCoeff = Y(1:numModesC,:);
    outputs.phiModeCoeff2 = Y(numModesC+1:2*numModesC,:);
    
    outputs.etaModeCoeff = Y(2*numModesC + 1 : 2*numModesC + numModesH,:);
    outputs.etaModeCoeff2 = Y(2*numModesC + numModesH + 1 : 2*numModesC + 2*numModesH,:);

    outputs.VEC = Y(2*numModesC + 2*numModesH + 6*nodesH + 1 : 2*numModesC + 2*numModesH + 6*nodesH + nodesC,:);
    outputs.VEC2 = Y(2*numModesC + 2*numModesH + 6*nodesH + nodesC + 1 : 2*numModesC + 2*numModesH + 6*nodesH + 2*nodesC,:);
    
    outputs.VRC = Y(2*numModesC + 2*numModesH + 6*nodesH + 2*nodesC + 1 : 2*numModesC + 2*numModesH + 6*nodesH + 3*nodesC,:);
    outputs.VRC2 = Y(2*numModesC + 2*numModesH + 6*nodesH + 3*nodesC + 1 : 2*numModesC + 2*numModesH + 6*nodesH + 4*nodesC,:);
    
    outputs.VSC = Y(2*numModesC + 2*numModesH + 6*nodesH + 4*nodesC + 1 : 2*numModesC + 2*numModesH + 6*nodesH + 5*nodesC,:);
    outputs.VSC2 = Y(2*numModesC + 2*numModesH + 6*nodesH + 5*nodesC + 1 : 2*numModesC + 2*numModesH + 6*nodesH + 6*nodesC,:);
    
    outputs.VEH = Y(2*numModesC + 2*numModesH + 1:2*numModesC + 2*numModesH + nodesH,:);
    outputs.VEH2 = Y(2*numModesC + 2*numModesH + nodesH + 1 : 2*numModesC + 2*numModesH + 2*nodesH,:);
    
    outputs.VRH = Y(2*numModesC + 2*numModesH + 2*nodesH + 1 : 2*numModesC + 2*numModesH + 3*nodesH,:);
    outputs.VRH2 = Y(2*numModesC + 2*numModesH + 3*nodesH + 1 : 2*numModesC + 2*numModesH + 4*nodesH,:);
    
    outputs.VSH = Y(2*numModesC + 2*numModesH + 4*nodesH + 1 : 2*numModesC + 2*numModesH + 5*nodesH,:);
    outputs.VSH2 = Y(2*numModesC + 2*numModesH + 5*nodesH + 1 : 2*numModesC + 2*numModesH + 6*nodesH,:);

    Y0 = [
        outputs.phiModeCoeff(:, end);
        outputs.phiModeCoeff2(:, end);
        outputs.etaModeCoeff(:, end);
        outputs.etaModeCoeff2(:, end);
        outputs.VEH(:, end);
        outputs.VEH2(:, end);
        outputs.VRH(:, end);
        outputs.VRH2(:, end);
        outputs.VSH(:, end);
        outputs.VSH2(:, end);
        outputs.VEC(:, end);
        outputs.VEC2(:, end);
        outputs.VRC(:, end);
        outputs.VRC2(:, end);
        outputs.VSC(:, end);
        outputs.VSC2(:, end)
    ];
end
