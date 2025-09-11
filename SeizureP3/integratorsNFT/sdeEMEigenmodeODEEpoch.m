function [sol, paramsS] = sdeEMEigenmodeODEEpoch(sdefun, sdegun, tspan, y0, options,...
    params, contn, Y)
if isfield(options,'InitialStep')
    dt = options.InitialStep;       % may return dt=[]
else
    dt = [];                        % use an empty time step
end

% If the given time step is empty then ....
if isempty(dt)
    dt = 0.001;     % default time step
    warning('sdeEM:InitialStep','Step size is undefined. Using InitialStep=%g',dt);
end

if isfield(options,'randn') && ~isempty(options.randn)
    % The rows of dW determine the number of noise sources (m).
    % The cols of dW determine the number of time steps (tcount).
    [mm,tcount] = size(options.randn);

    % The number of time steps in tspan determine the step size (dt)
    dt = (tspan(end) - tspan(1))/(tcount-1);

    % Assert that m matches the NoiseSources option
    assert(mm==nmodesC,'The number of rows in options.randn must equal options.NoiseSources');
end

% span the time domain in fixed steps
sol.x = tspan(1):dt:tspan(end);
tcount = numel(sol.x);

tauC = params.tauC;
tauH = params.tauH;
numModesC = params.numModesC;
numModesH = params.numModesH;

hippoEigenvectors = params.hippoEigenvectors;
corticalEigenvectors = params.corticalEigenvectors;

delaytimeC = round((tauC/dt));
delaytimeH = round((tauH/dt));
Yi = Y;
% delaytime in samples;
% allocate space for the results
sol.y = Yi;                           % values of y(t)
sol.yp = sol.y;                       % values of dy(t)/dt 

% compute the Wiener increments
if isfield(options,'randn') && ~isempty(options.randn)
    sol.dW = sqrt(dt) .* options.randn;
else
    sol.dW = generateNoiseCH(tcount, dt, numModesC, numModesH, ...
                params.nodesC, params.nodesH, corticalEigenvectors,...
                hippoEigenvectors, 5e-4, params.noiseType);
end
% miscellaneous output
sol.solver = mfilename;
sol.extdata.odefun = sdefun;
sol.extdata.sdefun = sdegun;
sol.extdata.options = options;
sol.ex21tdata.params = params;

% Get the OutputFcn callback.
OutputFcn = odeget(options,'OutputFcn');

% Implementing the Fixed-step Euler method
sol.y(:,1) = y0;
[sol, paramsS] = fixedStepEulerImplementation(sdefun, sdegun, sol, tcount, delaytimeC, delaytimeH,...
    numModesC, numModesH, params.nodesH, params.nodesC, dt, params, contn);

end
