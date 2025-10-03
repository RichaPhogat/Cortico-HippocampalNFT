function [nu_seC, nu_seH, cc, ch] = updateCouplingAndGains(params, contn, tcount, tSpan)
    % Persistent storage
    persistent nuseCL nuseHL ccL chL

    if isempty(nuseCL) % Initialize persistent variables
        nuseCL = params.nu_seC;
        nuseHL = params.nu_seH;
        ccL = 2.8*0.008;
        chL = 0.008;
    end

    % Define time vector
    tvec = linspace(0, tSpan, tcount);

    % Threshold checks
   countCheck0 = 110/tSpan;%110%5;
    contnCheck1 = 120/tSpan;
    contnCheck2 = 135/tSpan;
    contnCheck3 = 160/tSpan;
    contnCheck4 = 173/tSpan;
    contnCheck5 = 183/tSpan;
    slopeIncr = (0.01)/1000;
    slopeReduc = slopeIncr;
    
    if contn > countCheck0
        
        if contn <= contnCheck1
            nu_seC = (1.4/1000) * ones(1, tcount);
            nu_seH = nuseHL + 100 * slopeIncr * tvec;
            cc = ccL + 0 * slopeReduc * tvec;
            ch = chL + 0 * slopeReduc * tvec;
        else
            nu_seC = (1.4/1000) * ones(1, tcount);
            nu_seH = nuseHL - 38 * slopeReduc * tvec;
            cc = ccL + 0 * slopeReduc * tvec;
            ch = chL + 0 * slopeReduc * tvec;
        end
    
        if contn > (contnCheck2) && contn <= (contnCheck3-(7/tSpan))
            nu_seH = nuseHL + 0.0 * slopeReduc * tvec;
            nu_seC = nuseCL + 0.00 * slopeReduc * tvec;
            cc = ccL + 2.8*75 * slopeReduc * tvec;
            ch = chL + 75 * slopeReduc * tvec;
        end
    
        if contn > (contnCheck3-(7/tSpan)) && contn <= (contnCheck3)
            nu_seH = nuseHL + 0.0 * slopeReduc * tvec;
            nu_seC = nuseCL + 0.0 * slopeReduc * tvec;
            cc = ccL + 0 * slopeReduc * tvec;
            ch = chL + 0 * slopeReduc * tvec;
        end
    
        if contn > (contnCheck3) && contn <= (contnCheck3+(3/tSpan)) 
            nu_seH = nuseHL + 20 * slopeReduc * tvec;
            nu_seC = nuseCL + 0.00 * slopeReduc * tvec;
            cc = ccL + 0 * slopeReduc * tvec;
            ch = chL + 0 * slopeReduc * tvec;
        end
    
        if contn > (contnCheck3+(3/tSpan)) && contn <= (contnCheck4) 
            nu_seH = nuseHL + 0 * slopeReduc * tvec;
            nu_seC = nuseCL + 0.00 * slopeReduc * tvec;
            cc = ccL + 0 * slopeReduc * tvec;
            ch = chL + 0 * slopeReduc * tvec;
        end
    
        if contn > (contnCheck4) && contn <= (contnCheck5) 
            nu_seH = nuseHL - 36 * slopeReduc * tvec;
            nu_seC = nuseCL + 4.2 * slopeReduc * tvec;
            cc = ccL - 2.8*75 * slopeReduc * tvec;
            ch = chL - 75 * slopeReduc * tvec;
        end
    
        if contn > (contnCheck5) && contn <= (contnCheck5+(6/tSpan)) 
            nu_seH = nuseHL - 0 * slopeReduc * tvec;
            nu_seC = nuseCL - 7 * slopeReduc * tvec;
            cc = ccL - 2.8*75 * slopeReduc * tvec;
            ch = chL - 75 * slopeReduc * tvec;
        end
    
        if contn > (contnCheck5+(6/tSpan)) && contn <= (contnCheck5+(9/tSpan)) 
            nu_seH = nuseHL - 0 * slopeReduc * tvec;
            nu_seC = nuseCL - 0.0 * slopeReduc * tvec;
            cc = ccL - 0 * slopeReduc * tvec;
            ch = chL - 0 * slopeReduc * tvec;
        end
    
        if contn > (contnCheck5+(9/tSpan))
            nu_seH = nuseHL - 10 * slopeReduc * tvec;
        end
    
    %         cc = ccL + 10 * slopeReduc * tvec;
    %         ch = chL + 10 * slopeReduc * tvec;    
    
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
