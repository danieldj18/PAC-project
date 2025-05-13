function [MI, meanHFampPerLFphaseBin, nSamplesPerLFphaseBin] = calcModulationIndex_Tort2010(LFphase, HFamp, nPhaseBins)
    
    LFphase_binIds = discretize(LFphase, linspace(-pi, pi, nPhaseBins+1));

    meanHFampPerLFphaseBin = nan(1, nPhaseBins);  
    nSamplesPerLFphaseBin = nan(1, nPhaseBins);   

    for binIdx = 1:nPhaseBins

        meanHFampPerLFphaseBin(binIdx) = mean(HFamp(LFphase_binIds == binIdx));

        nSamplesPerLFphaseBin(binIdx) = length(HFamp(LFphase_binIds == binIdx));
    end

    normMeanHFampPerLFphaseBin = meanHFampPerLFphaseBin / sum(meanHFampPerLFphaseBin);

    if any(normMeanHFampPerLFphaseBin == 0)
        normMeanHFampPerLFphaseBin(normMeanHFampPerLFphaseBin == 0) = eps;
    end

    P = normMeanHFampPerLFphaseBin;

    N = nPhaseBins;

    each_bin_entropy = nan(1, N);  
    
    for j = 1:N
        each_bin_entropy(j) = P(j) * log(P(j));
    end

    entropy = -sum(each_bin_entropy);

    D_KL = log(N) - entropy;

    MI = D_KL / log(N);

end






