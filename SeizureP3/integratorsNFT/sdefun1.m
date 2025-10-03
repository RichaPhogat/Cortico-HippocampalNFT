function dg = sdefun1(params, mC, mH, numModesC, numModesH)

    betaC = params.alphaC*params.betabyalphaC;
    betaH = params.alphaH*params.betabyalphaH;
    %dg=sparse(6*n+2*num_modes,6*n+2*num_modes);
    %dg(2*num_modes+5*n+1:6*n+2*num_modes,2*num_modes+5*n+1:6*n+2*num_modes)=alpha*beta*nu_sn.*[speye(n,n)];
    
    dg = sparse(2*numModesC + 2*numModesH + 6*mH + 6*mC, 2*numModesC + ...
        2*numModesH + 6*mH + 6*mC);
    dg(2*numModesC + 2*numModesH + 5*mH + 1 : 2*numModesC + 2*numModesH + ...
        6*mH, 2*numModesC + 2*numModesH + 5*mH + 1 : 2*numModesC + 2*numModesH + 6*mH) = ...
        20*params.alphaH*betaH*params.nu_snH.*speye(mH,mH);
    dg(2*numModesC + 2*numModesH + 6*mH + 5*mC + 1 : 2*numModesC + 2*numModesH + 6*mH + 6*mC, 2*numModesC ...
        + 2*numModesH + 6*mH + 5*mC + 1 : 2*numModesC + 2*numModesH + 6*mH + 6*mC) = ...
        150*params.alphaC*betaC*params.nu_snC.*speye(mC,mC);

end