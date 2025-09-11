function dW = generateNoiseCH(tcount, dt, numModesC, numModesH, ...
    mC, mH, corticalEigenvectors, hippoEigenvectors, scale, noiseType)

    %%% Spatiotemporal Noise Generation %%%
    switch noiseType
    case 'eigenmode'
        %%% Eigenmode-based noise generation
        noiseHippo = zeros(mH, tcount);
        for eigVal = 1:numModesH
            noiseHippo = noiseHippo + hippoEigenvectors(:,eigVal)*randn(1,tcount);
        end

        noiseCortex = zeros(mC, tcount);
        for eigVal = 1:numModesC
            noiseCortex = noiseCortex + corticalEigenvectors(:,eigVal)*randn(1,tcount);
        end

    case 'gaussian'
        %%% Simple Gaussian noise per node
        noiseHippo = randn(mH, tcount);
        noiseCortex = randn(mC, tcount);

    otherwise
        error('Unknown noiseType: %s. Use ''eigenmode'' or ''gaussian''.', params.noiseType);
    end


    % Build sparse noise matrix
    nn = sparse(2*numModesC + 2*numModesH + 6*mH + 6*mC , tcount);
    nn(2*numModesC + 2*numModesH + 5*mH+1 : 2*numModesC + 2*numModesH + 6*mH , :) = noiseHippo;
    nn(2*numModesC + 2*numModesH + 6*mH + 5*mC+1 : 2*numModesC + 2*numModesH + 6*mH + 6*mC , :) = noiseCortex;

    % Scale Wiener increments
    dW = sqrt(scale * dt) .* nn;
end
