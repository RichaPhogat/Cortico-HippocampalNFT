function  dY = odefun1(~, Y, YDelayed, params, nu_seC, nu_seH, cc, ch)
    
    betaC = params.alphaC*params.betabyalphaC;
    betaH = params.alphaH*params.betabyalphaH;
    gammaC = (params.vC)/(params.reC);
    gammaH = (params.vH)/(params.reH);

    YCortexMode = Y(1 : 2*params.numModesC);
    YHippoMode = Y(2*params.numModesC + 1 : 2*params.numModesC + 2*params.numModesH);
    YCortex = Y(2*params.numModesC + 2*params.numModesH + 6*params.nodesH ...
        + 1 : 2*params.numModesC + 2*params.numModesH + 6*params.nodesH + 6*params.nodesC);   
    YHippo = Y(2*params.numModesC + 2*params.numModesH + 1: 2*params.numModesC...
        + 2*params.numModesH + 6*params.nodesH);

    YCortexModeDelayed = YDelayed(1 : 2*params.numModesC);
    YHippoModeDelayed = YDelayed(2*params.numModesC + 1 : 2*params.numModesC...
        + 2*params.numModesH);
    YCortexDelayed = YDelayed(2*params.numModesC + 2*params.numModesH + 6*params.nodesH...
        + 1 : 2*params.numModesC + 2*params.numModesH + 6*params.nodesH + 6*params.nodesC);   
    YHippoDelayed = YDelayed(2*params.numModesC + 2*params.numModesH + 1: 2*params.numModesC...
        + 2*params.numModesH + 6*params.nodesH);
    YHippoModeDelayedC = YDelayed(2*params.numModesC + 2*params.numModesH + 6*params.nodesH...
        + 6*params.nodesC +1: 2*params.numModesC + 2*params.numModesH + 6*params.nodesH + 6*params.nodesC...
        + 2*params.numModesH);

    % restore 2D data format

    VSCortexDelayed = YCortexDelayed(4*params.nodesC+1:5*params.nodesC);
    VSHippoDelayed = YHippoDelayed(4*params.nodesH+1:5*params.nodesH);

    phiMode = YCortexMode(1:params.numModesC);
    phiMode2 = YCortexMode(params.numModesC+1:2*params.numModesC);
    etaMode = YHippoMode(1:params.numModesH);
    etaMode2 = YHippoMode(params.numModesH+1:2*params.numModesH);

    phiE = params.corticalEigenvectors*phiMode;
    phiEdelayed = params.corticalEigenvectors*YCortexModeDelayed(1:params.numModesC);
    eta = params.hippoEigenvectors*etaMode;                %<============================
    etadelayed = params.hippoEigenvectors*YHippoModeDelayed(1:params.numModesH); %<============================
    etadelayedC = params.hippoEigenvectors*YHippoModeDelayedC(1:params.numModesH);

    VEC = YCortex(1:params.nodesC);
    VEC2 = YCortex(params.nodesC+1:2*params.nodesC);
    VEH = YHippo(1:params.nodesH);
    VEH2 = YHippo(params.nodesH+1:2*params.nodesH);

    VRC = YCortex(2*params.nodesC+1:3*params.nodesC);
    VRC2 = YCortex(3*params.nodesC+1:4*params.nodesC);
    VRH = YHippo(2*params.nodesH+1:3*params.nodesH);
    VRH2 = YHippo(3*params.nodesH+1:4*params.nodesH);

    VSC = YCortex(4*params.nodesC+1:5*params.nodesC);
    VSC2 = YCortex(5*params.nodesC+1:6*params.nodesC);
    VSH = YHippo(4*params.nodesH+1:5*params.nodesH);
    VSH2 = YHippo(5*params.nodesH+1:6*params.nodesH);

    QEC = params.QmaxC./(1 + exp(-(VEC-params.thetaC)/params.sigmaC));
    QIC = QEC;
    QRC = params.QmaxC./(1 + exp(-(VRC-params.thetaC)/params.sigmaC));
    QSC = params.QmaxC./(1 + exp(-(VSC-params.thetaC)/params.sigmaC));

    QEH = params.QmaxH./(1 + exp(-(VEH-params.thetaH)/params.sigmaH));
    QIH = QEH;
    QRH = params.QmaxH./(1 + exp(-(VRH-params.thetaH)/params.sigmaH));
    QSH = params.QmaxH./(1 + exp(-(VSH-params.thetaH)/params.sigmaH));

    QSCdelayed= params.QmaxC./(1 + exp(-(VSCortexDelayed-params.thetaC)/params.sigmaC));
    QSHdelayed= params.QmaxH./(1 + exp(-(VSHippoDelayed-params.thetaH)/params.sigmaH));

    % if contn == 3
    %   keyboard
    % end

    QphiEModes = calc_eigendecomposition(QEC + cc*params.couplingCortexToHippo*etadelayedC,...
        params.corticalEigenvectors, 'matrix');
    QetaModes = calc_eigendecomposition(QEH + ch*params.couplingHippoToCortex*phiEdelayed,...
        params.hippoEigenvectors, 'matrix');

    QPhiError = QEC - params.corticalEigenvectors*QphiEModes;
    QEtaError = QEH - params.hippoEigenvectors*QetaModes;
    QEtaRelative = abs(QEtaError./QEH);
    QPhiRelative = abs(QPhiError./QEC);
    QEtaRelativePercentage = max(QEtaRelative*100);
    QPhiRelativePercentage = max(QPhiRelative*100);

    % QPhi1 = corticalEigenvectors*QphiEModes;
    % QEta1 = hippoEigenvectors*QetaModes;
    % cc1 = corrcoef(QEC,QPhi1);
    % cc2 = corrcoef(QEH,QEta1);
    % QEtaRelativePercentage = cc1(1,2);
    % QPhiRelativePercentage = cc2(1,2);


    dphiMode = phiMode2;
    dphiMode2 = gammaC^2 * (QphiEModes - (2 / gammaC) * phiMode2...
        - phiMode .* (1 + params.reC^2 * params.corticalEigenvalues));

    detaMode = etaMode2;
    detaMode2 = gammaH^2 * (QetaModes - (2 / gammaH) * etaMode2...
        - etaMode .* (1 + params.reH^2 * params.hippoEigenvalues));

    PeC = params.nu_eeC*phiE + params.nu_eiC*QIC + params.nu_esC*QSCdelayed;
    PrC = params.nu_rsC*QSC + params.nu_reC*phiEdelayed;
    PsC = nu_seC*phiEdelayed + params.nu_srC*QRC ...
        + params.nu_snC*(params.phi_n0C);   

    PeH = params.nu_eeH*eta + params.nu_eiH*QIH + params.nu_esH*QSHdelayed;
    PrH = params.nu_rsH*QSH + params.nu_reH*etadelayed;
    PsH = nu_seH*etadelayed + params.nu_srH*QRH ...
        + params.nu_snH*(params.phi_n0C); 

    dVEC = VEC2;
    dVEC2 = -(betaC + params.alphaC)*VEC2 ...
        - params.alphaC*betaC*VEC + params.alphaC*betaC*PeC;

    dVRC = VRC2;
    dVRC2 = -(betaC + params.alphaC)*VRC2 ...
        - params.alphaC*betaC*VRC + params.alphaC*betaC*PrC;

    dVSC = VSC2;
    dVSC2 = -(betaC + params.alphaC)*VSC2 ...
        - params.alphaC*betaC*VSC + params.alphaC*betaC*PsC; % + stochastic input

    dVEH = VEH2;
    dVEH2 = -(betaH + params.alphaH)*VEH2 ...
        - params.alphaH*betaH*VEH + params.alphaH*betaH*PeH;

    dVRH = VRH2;
    dVRH2 = -(betaH + params.alphaH)*VRH2 ...
        - params.alphaH*betaH*VRH + params.alphaH*betaH*PrH;

    dVSH = VSH2;
    dVSH2 = -(betaH + params.alphaH)*VSH2 ...
        - params.alphaH*betaH*VSH + params.alphaH*betaH*PsH; % + stochastic input

    dY=[dphiMode;dphiMode2;detaMode;detaMode2;dVEH;dVEH2;dVRH;dVRH2;dVSH;dVSH2;...
        dVEC;dVEC2;dVRC;dVRC2;dVSC;dVSC2;QPhiRelativePercentage;QEtaRelativePercentage];
end