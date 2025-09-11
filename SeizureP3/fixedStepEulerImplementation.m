function [sol, params] = fixedStepEulerImplementation(sdefun, sdegun, sol, tcount, delaytimeC,...
    delaytimeH, numModesC, numModesH, nodesH, nodesC, dt, params, contn)
   
   [nu_seC, nu_seH, cc, ch] = updateCouplingAndGains(params, contn, tcount, sol.x(end));

    %%%%%%%%%%%%%%%%% Delay Implementation %%%%%%%%%%%%%%%%
    for indx=1:tcount-1

            params.nu_seC(indx) = nu_seC(indx);
            params.nu_seH(indx) = nu_seH(indx);
            params.cc(indx) = cc(indx);
            params.ch(indx) = ch(indx);

        if delaytimeC ~= 0
            if indx - delaytimeC > 0
                yModeCortexDelayed = sol.y(1 : 2*numModesC, indx - delaytimeC);
                yCortexDelayed = sol.y(2*numModesC + 2*numModesH + 6*nodesH + 1 : 2*numModesC...
                    + 2*numModesH + 6*nodesH + 6*nodesC, indx - delaytimeC);
                yModeHippoDelayedC = sol.y(2*numModesC + 1 : 2*numModesC...
                    + 2*numModesH, indx - delaytimeC);
    
            else
                yModeCortexDelayed = sol.y(1 : 2*numModesC, 1);
                yCortexDelayed = sol.y(2*numModesC + 2*numModesH + 6*nodesH + 1 : 2*numModesC...
                    + 2*numModesH + 6*nodesH + 6*nodesC, 1);
                yModeHippoDelayedC = sol.y(2*numModesC + 1 : 2*numModesC + 2*numModesH, 1);
            end
        else
            yCortexDelayed = sol.y*0;
        end
    
        if delaytimeH ~= 0
            if indx - delaytimeH > 0
                yModeHippoDelayed = sol.y(2*numModesC + 1 : 2*numModesC + 2*numModesH,...
                    indx - delaytimeH);
                yHippoDelayed = sol.y(2*numModesC + 2*numModesH + 1 : 2*numModesC +...
                    2*numModesH + 6*nodesH, indx - delaytimeH);
            else
                yModeHippoDelayed = sol.y(2*numModesC + 1 : 2*numModesC + 2*numModesH, 1);
                yHippoDelayed = sol.y(2*numModesC + 2*numModesH + 1 : 2*numModesC...
                    + 2*numModesH + 6*nodesH, 1);
            end
        else
            yCortexDelayed = sol.y*0;
        end
        
        yDelayed = [yModeCortexDelayed; yModeHippoDelayed; yHippoDelayed; yCortexDelayed;...
            yModeHippoDelayedC];
        
            Fe = sdefun(sol.x(indx), sol.y(:,indx), yDelayed, params, nu_seC(1,indx),...
                nu_seH(1,indx), cc(1,indx), ch(1,indx));
            G = sdegun(params, nodesC, nodesH, numModesC, numModesH);
         
        
        Loss = Fe(end-1:end);
        F = Fe(1:end-2);
    
        % Euler step
        sol.yp(:,indx) = F*dt + G*sol.dW(:,indx);            % dy(t) = F(t,y(t))*dt + G(t,y(t))*dW(t)
        sol.y(:,indx+1) = sol.y(:,indx) + sol.yp(:,indx);    % y(t+1) = y(t) + dy(t)
        sol.Loss(1:2,indx+1) = Loss;
    end

    %%%%%%%%%%%%%%%%%%%%%% Last Step %%%%%%%%%%%%%%%%%%%%%%%%%
    if delaytimeC ~= 0

        if tcount - delaytimeC > 0
            yModeCortexDelayed = sol.y(1:2*numModesC, tcount-delaytimeC);
            yCortexDelayed = sol.y(2*numModesC + 2*numModesH + 6*nodesH + 1 : 2*numModesC...
                + 2*numModesH + 6*nodesH + 6*nodesC, tcount-delaytimeC);
            yModeHippoDelayedC = sol.y(2*numModesC + 1 : 2*numModesC + 2*numModesH,...
                tcount-delaytimeC);
        else
            yModeCortexDelayed = sol.y(1:2*numModesC,1);
            yCortexDelayed = sol.y(2*numModesC + 2*numModesH + 6*nodesH + 1 : 2*numModesC...
                + 2*numModesH + 6*nodesH + 6*nodesC, 1);
            yModeHippoDelayedC = sol.y(2*numModesC + 1 : 2*numModesC + 2*numModesH, 1);
        end
    else
        yCortexDelayed = yn*0;
    end

    if delaytimeH ~= 0
        if tcount - delaytimeH > 0
            yModeHippoDelayed = sol.y(2*numModesC + 1 : 2*numModesC + 2*numModesH, ...
                tcount-delaytimeH);
            yHippoDelayed = sol.y(2*numModesC + 2*numModesH + 1 : 2*numModesC +...
                2*numModesH + 6*nodesH, tcount-delaytimeH);
        else
            yModeHippoDelayed = sol.y(2*numModesC + 1 : 2*numModesC + 2*numModesH,1);
            yHippoDelayed = sol.y(2*numModesC + 2*numModesH + 1 : 2*numModesC...
                + 2*numModesH + 6*nodesH, 1);
        end
    else
        yHippoDelayed = yn*0;
    end

yDelayed = [yModeCortexDelayed; yModeHippoDelayed; yHippoDelayed; yCortexDelayed;...
    yModeHippoDelayedC];

% Complete the final Euler step
Fe = sdefun(sol.x(end), sol.y(:,end), yDelayed, params, nu_seC(1,indx), nu_seH(1,indx),...
    cc(1,end), ch(1,end));
G = sdegun(params, nodesC, nodesH, numModesC, numModesH);
Loss = Fe(end-1:end);
F = Fe(1:end-2);
sol.yp(:,end) = F*dt + G*sol.dW(:,end);
sol.Loss(1:2,end) = Loss;

end
