function [nu_seC, nu_seH, cc, ch] = updateCouplingAndGains(params, contn, tcount, tSpan)
    % Persistent storage
    persistent nuseCL nuseHL ccL chL

    if isempty(nuseCL) % Initialize persistent variables
        nuseCL = params.nu_seC;
        nuseHL = params.nu_seH;
        ccL = 2.8*0.007;
        chL = 0.007;
    end

    % Define time vector
    tvec = linspace(0, tSpan, tcount);
    

    countCheck0 = 35/tSpan;
    contnCheck1 = 41/tSpan;
    contnCheck2 = 59/tSpan;
    contnCheck3 = 110/tSpan;
    contnCheck4 = 146/tSpan;
    slopeIncr = (0.01)/1000;
    slopeReduc = slopeIncr;
    
    if contn > countCheck0
    
    
        if contn <= contnCheck1
            nu_seH = nuseHL - 50.2 * slopeIncr * tvec;  %37 init change 13.2 or 1.32
            nu_seC = nuseCL - 0.1 * slopeIncr * tvec;
            cc = ccL - 2.8*60* slopeReduc * tvec;
            ch = chL - 60* slopeReduc * tvec;    
    
        else
            nu_seH = nuseHL + 165 * slopeIncr * tvec;
            nu_seC = nuseCL + 0.1 * slopeReduc * tvec;
            cc = ccL + 2.8*135* slopeReduc * tvec;
            ch = chL + 135* slopeReduc * tvec;    
        end
    
        if contn > (contnCheck2) && contn <= (contnCheck3)
            nu_seH = nuseHL - 1.5 * slopeIncr * tvec;
            nu_seC = nuseCL - 0.03 * slopeIncr * tvec;
            cc = ccL - 2.8*2.25* slopeReduc * tvec;
            ch = chL - 2.25* slopeReduc * tvec;    
        end
    
        if contn > (contnCheck1+(3/tSpan)) && contn <= (contnCheck2)
            nu_seH = nuseHL + 0 * slopeIncr * tvec;   
        end
    
        if contn > (contnCheck2-(6/tSpan)) && contn <= (contnCheck2+(6/tSpan))
            nu_seH = nuseHL + 0 * slopeIncr * tvec;
            cc = ccL - 2.8*8.5* slopeReduc * tvec;
            ch = chL - 8.5 * slopeReduc * tvec;    
        end
    
        if contn > (contnCheck2+(6/tSpan)) && contn <= (contnCheck2+(12/tSpan))
            nu_seH = nuseHL + 100 * slopeIncr * tvec;
            cc = ccL - 2.8*168.75 * slopeReduc * tvec;
            ch = chL - 168.75 * slopeReduc * tvec;    
        end
    
        if contn > (contnCheck2+(12/tSpan)) && contn <= (contnCheck2+(18/tSpan))
            nu_seH = nuseHL - 100 * slopeIncr * tvec;
            cc = ccL + 2.8*160* slopeReduc * tvec;
            ch = chL + 160 * slopeReduc * tvec;    
        end
    
        if contn > (contnCheck2+(18/tSpan)) && contn <= (contnCheck2+(24/tSpan))
            nu_seH = nuseHL + 76 * slopeIncr * tvec;
            cc = ccL - 2.8*67.5* slopeReduc * tvec;
            ch = chL - 67.5 * slopeReduc * tvec;
        end
    
        if contn > (contnCheck2+(24/tSpan)) && contn <= (contnCheck3)
            nu_seH = nuseHL - 30 * slopeIncr * tvec;
            cc = ccL - 2.8*15.5* slopeReduc * tvec;
            ch = chL - 15.5* slopeReduc * tvec;
        end
    
        if contn > (contnCheck3) && contn <= (contnCheck4)
            nu_seC = nuseCL + 0.02 * slopeReduc * tvec;
            nu_seH = nuseHL + 8 * slopeReduc * tvec;
            cc = ccL + 2.8 * 4 * slopeReduc * tvec;
            ch = chL + 4 * slopeReduc * tvec;    
        end
    
        if contn > (contnCheck4) && contn <= (contnCheck4+(15/tSpan))
            nu_seC = nuseCL - 0.1 * slopeReduc * tvec;
            nu_seH = nuseHL - 70 * slopeReduc * tvec;
            cc = ccL - 2.8*45 * slopeReduc * tvec;
            ch = chL - 45 * slopeReduc * tvec;    
        end
    else
        nu_seC = nuseCL * ones(1, tcount);
        nu_seH = nuseHL * ones(1, tcount);
        cc = ccL*ones(1, tcount);
        ch = chL*ones(1, tcount);
    
        if contn > (countCheck0-(12/tSpan)) && contn <= (countCheck0-(9/tSpan))
            nu_seH = nuseHL + 248 * slopeReduc * tvec;
        end
        if contn > (countCheck0-(9/tSpan)) && contn < (countCheck0)
            nu_seH = nuseHL - 88 * slopeReduc * tvec;
        end
    
    
    end
    
    nuseCL = nu_seC(end);
    nuseHL = nu_seH(end);
    ccL = cc(end);
    chL = ch(end);

end
