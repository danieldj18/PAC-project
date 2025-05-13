function [MI, meanHFampPerLFphaseBin, nSamplesPerLFphaseBin] = calcModulationIndex_Tort2010(LFphase, HFamp, nPhaseBins)

    %{
    
    Inputs:
    - LFphase: Vector of low-frequency phase values (radians).
    - HFamp: Vector of high-frequency amplitude values.
    - nPhaseBins: The number of phase bins to categorize LFphase (e.g., 18 phase bins, 
                  corresponding to 18 equally spaced phase intervals between -π and π).
    
    Outputs:
    - MI: The Modulation Index, a value between 0 and 1. A value close to 1 indicates strong phase-amplitude coupling.
        - meanHFampPerLFphaseBin: The mean high-frequency amplitude for each phase bin.
        - nSamplesPerLFphaseBin: The number of samples in each phase bin.

    %}

    % Step 1: Bin the low-frequency (LF) phase values into discrete phase bins
    %{
    Bin the low-frequency phase values (LFphase) into nPhaseBins equally spaced bins between -π and π.
    Ex: If nPhaseBins = 4, the bins are [-π, -π/2, 0, π/2, π] <= the edges
    Each value in LFphase is assigned to one of these bins, and the bin index is stored in LFphase_binIds.
    Example: If LFphase = [-2.5, -1.5, 0.3, 1.5, 2.8] (in radians) and nPhaseBins = 4,
    the output LFphase_binIds = [1, 1, 2, 3, 4] indicates the bin index for each phase value.
    %}
    
    LFphase_binIds = discretize(LFphase, linspace(-pi, pi, nPhaseBins+1));

    % Initialize arrays to store the mean HF amplitude and the number of samples in each phase bin
    meanHFampPerLFphaseBin = nan(1, nPhaseBins);   % Mean high-frequency amplitude for each phase bin
    nSamplesPerLFphaseBin = nan(1, nPhaseBins);    % Number of samples in each phase bin

    % Step 2: Loop through each phase bin to compute the mean HF amplitude and the number of samples
    for binIdx = 1:nPhaseBins

        % Calculate the mean HF amplitude for the current bin (HFamp values where LFphase belongs to the current bin)
        meanHFampPerLFphaseBin(binIdx) = mean(HFamp(LFphase_binIds == binIdx));

        % Count the number of samples that fall into this bin
        nSamplesPerLFphaseBin(binIdx) = length(HFamp(LFphase_binIds == binIdx));
    end

    % Step 3: Normalize the mean HF amplitude per phase bin
    % This is done to compute the probability distribution of HF amplitude across the phase bins.
    % According to Tort et al., 2010, we normalize by dividing each bin's mean HF amplitude by the total sum.
    normMeanHFampPerLFphaseBin = meanHFampPerLFphaseBin / sum(meanHFampPerLFphaseBin);

    % Step 4: Handle cases where a phase bin might have a mean HF amplitude of zero
    % To avoid issues in the MI calculation (like log(0)), we replace zeros with a very small value (epsilon).
    if any(normMeanHFampPerLFphaseBin == 0)
        normMeanHFampPerLFphaseBin(normMeanHFampPerLFphaseBin == 0) = eps;
    end

    % Step 5: Calculate the Shannon Entropy (H) and Kullback-Leibler Divergence (D_KL)
    % Based on Tort et al., 2010, we first calculate the Shannon Entropy (H), 
    % which measures the "disorder" in the distribution of HF amplitude across phase bins.

    % P is the normalized mean HF amplitude, which can be treated as a probability distribution.
    P = normMeanHFampPerLFphaseBin;

    % N is the number of phase bins
    N = nPhaseBins;

    % Initialize entropy summation term (x) as a vector to store the P(j) * log(P(j)) for each bin
    each_bin_entropy = nan(1, N);  % Initialize to avoid errors for summation
    
    for j = 1:N
        % Compute P(j) * log(P(j)) for each phase bin
        each_bin_entropy(j) = P(j) * log(P(j));
    end

    % Shannon Entropy (H) is the negative sum of P(j) * log(P(j)) across all bins
    entropy = -sum(each_bin_entropy);

    % Step 6: Compute the Kullback-Leibler Divergence (D_KL)
    % D_KL is a measure of how "far" the distribution of HF amplitude per phase bin is from a uniform distribution.
    % A uniform distribution indicates no phase-amplitude coupling, whereas a non-uniform distribution suggests coupling.
    D_KL = log(N) - entropy;

    % Step 7: Compute the Modulation Index (MI)
    % The MI is the Kullback-Leibler Divergence normalized by log(N), which scales the MI between 0 and 1.
    MI = D_KL / log(N);

end






