function PAC = calculatePAC_mst_trial(fullSig, fullMask, segIdx, fs, ...
    LFrange, LFwidth, LFstep, HFrange, HFwidth, HFstep, ...
    nPhaseBins, nShuffles)

    % frequency‐band centers
    LF_centers = LFrange(1):LFstep:LFrange(2);
    HF_centers = HFrange(1):HFstep:HFrange(2);
    nLF = length(LF_centers);
    nHF = length(HF_centers);

    PAC.LF_frequencyCenters = LF_centers;
    PAC.HF_frequencyCenters = HF_centers;

    % take out the trial mask
    trialMask = fullMask(segIdx);

    % allocate prams
    observedMI = nan(nLF, nHF);
    MI_surr = nan(nLF, nHF, nShuffles);

    % loop over each LF‐HF pair
    for iLF = 1:nLF
        % filter + Hilbert for this LF on full signal
        lfWin = [LF_centers(iLF)-LFwidth/2, LF_centers(iLF)+LFwidth/2];
        lfFilt = zscore(eegfilt(fullSig, fs, lfWin(1), lfWin(2)));
        lfPhase = angle(hilbert(lfFilt));

        % select condition trial segs
        LF_phase = lfPhase(segIdx);
        LF_phaseTrial = LF_phase(trialMask);

        % brak theta phase into chunks for surrogates
        segLen = fs;
        total = length(LF_phaseTrial);
        nSegs = ceil(total/segLen);
        pad = nSegs * segLen - total;
        phiPad = [LF_phaseTrial nan(1, pad)];
        phiSeg = reshape(phiPad, segLen, nSegs);

        for iHF = 1:nHF

            % filter + Hilbert for this HF on full signal
            hfWin = [HF_centers(iHF) - HFwidth/2, HF_centers(iHF) + HFwidth/2];
            hfFilt = zscore(eegfilt(fullSig, fs, hfWin(1), hfWin(2)));
            hfAmpVec = abs(hilbert(hfFilt));

            % extract & mask trial-ampitude
            HF_amplitude = hfAmpVec(segIdx);
            HF_amplitudeTrial = HF_amplitude(trialMask);

            % observed MI
            observedMI(iLF,iHF) = calcModulationIndex_Tort2010(LF_phaseTrial, HF_amplitudeTrial, nPhaseBins);

            % surrogate MI
            for s = 1:nShuffles
                order = randperm(nSegs);
                phiShufM = phiSeg(:,order);
                phaseShuf = phiShufM(:);
                phaseShuf = phaseShuf(~isnan(phaseShuf));
                MI_surr(iLF,iHF,s) = calcModulationIndex_Tort2010(phaseShuf, HF_amplitudeTrial, nPhaseBins);
            end
        end
    end

    % z-score against surrogates
    surMean = mean(MI_surr, 3, 'omitnan');
    surStd = std( MI_surr, 0, 3, 'omitnan');
    PAC.MI = (observedMI - surMean) ./ surStd;
end
