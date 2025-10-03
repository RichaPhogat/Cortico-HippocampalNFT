function [nu_seC, nu_seH, cc, ch] = updateCouplingAndGains(params, contn, tcount, tSpan)
    % Persistent storage
    persistent nuseCL nuseHL ccL chL

    if isempty(nuseCL) % Initialize persistent variables
        nuseCL = params.nu_seC;
        nuseHL = params.nu_seH;
        ccL = 2.8*0.001;
        chL = 0.001;
    end

    % Define time vector
    tvec = linspace(0, tSpan, tcount);

    % Threshold checks
    countCheck0 = 47/tSpan;
    contnCheck1 = 56/tSpan;
    contnCheck2 = 68/tSpan;
    countCheck3 = 128/tSpan;
    slopeIncr = (0.01)/1000;
    slopeReduc = slopeIncr;

    if contn > countCheck0
        
        if contn <= contnCheck1
            nu_seH = nuseHL + 116 * slopeIncr * tvec;
            nu_seC = nuseCL + 0 * slopeIncr * tvec;
            cc = ccL + 0 * slopeIncr * tvec;
            ch = chL + 0 * slopeIncr * tvec;
        else
            nu_seH = nuseHL - 55 * slopeReduc * tvec;
            nu_seC = nuseCL + 1.584*slopeIncr * tvec;
            cc = ccL + 0 * slopeIncr * tvec;
            ch = chL + 0 * slopeIncr * tvec;
        end
    
        if contn > contnCheck2
            nu_seH = nuseHL + 25.5 * slopeReduc * tvec;
            cc = ccL + 2.8*120 * slopeIncr * tvec;
            ch = chL + 120 * slopeIncr * tvec;
            nu_seC = nuseCL - 0.0*slopeIncr * tvec;
        end
    
        if contn > (contnCheck2) && contn <= (contnCheck2+(7/tSpan))
            nu_seC = nuseCL - 1.58571*slopeIncr * tvec;
        end
    
        if contn > (contnCheck2 + (18/tSpan)) && contn <= (contnCheck2+(27/tSpan))
            cc = ccL - 2.8*100 * slopeReduc * tvec;
            ch = chL - 100 * slopeReduc * tvec;
        end
    
        if contn > (contnCheck2 + (27/tSpan)) && contn <= (contnCheck2+(33/tSpan))
            cc = ccL + 2.8*150 * slopeReduc * tvec;
            ch = chL + 150 * slopeReduc * tvec;
        end
    
        if contn > (contnCheck2 + (33/tSpan)) && contn <= (contnCheck2+(39/tSpan))
            nu_seH = nuseHL - 2 * slopeReduc * tvec;
            cc = ccL - 2.8*150 * slopeReduc * tvec;
            ch = chL - 150 * slopeReduc * tvec;
        end
    
        if contn > (contnCheck2 + (39/tSpan)) && contn <= (contnCheck2+(45/tSpan))
            nu_seH = nuseHL - 2 * slopeReduc * tvec;
            cc = ccL + 2.8*150 * slopeReduc * tvec;
            ch = chL + 150 * slopeReduc * tvec;
        end
    
        if contn > (contnCheck2 + (45/tSpan)) && contn <= (contnCheck2+(60/tSpan))
            cc = ccL - 2.8*35 * slopeReduc * tvec;
            ch = chL - 35 * slopeReduc * tvec;
            nu_seH = nuseHL - 60 * slopeReduc * tvec;
        end
    
        if contn > (countCheck3)
            cc = ccL - 2.8*20 * slopeReduc * tvec;
            ch = chL - 20 * slopeReduc * tvec;
            nu_seH = nuseHL - 3.5 * slopeReduc * tvec;
        end
    
        if contn > (countCheck3 + (33/tSpan)) && contn <= (countCheck3 + (43/tSpan))
            nu_seH = nuseHL - 42 * slopeReduc * tvec;
            cc = ccL - 2.8*106 * slopeReduc * tvec;
            ch = chL - 106 * slopeReduc * tvec;
        end
    
    
    else
        nu_seC = nuseCL * ones(1, tcount);
        nu_seH = nuseHL * ones(1, tcount);
        cc = ccL*ones(1, tcount);
        ch = chL*ones(1, tcount);
    end
    
    nuseCL = nu_seC(end);
    nuseHL = nu_seH(end);
    ccL = cc(end);
    chL = ch(end);

end
