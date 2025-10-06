function [nu_seC, nu_seH, cc, ch] = updateCouplingAndGains(params, contn, tcount, tSpan)
    % Persistent storage
    persistent nuseCL nuseHL ccL chL

    if isempty(nuseCL) % Initialize persistent variables
        nuseCL = params.nu_seC;
        nuseHL = params.nu_seH;
        ccL = 0.0;
        chL = 0.0;
    end


        nu_seC = nuseCL * ones(1, tcount);
        nu_seH = nuseHL * ones(1, tcount);
        cc = ccL*ones(1, tcount);
        ch = chL*ones(1, tcount);
    
    
    nuseCL = nu_seC(end);
    nuseHL = nu_seH(end);
    ccL = cc(end);
    chL = ch(end);

end
