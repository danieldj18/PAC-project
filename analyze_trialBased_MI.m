function results = analyze_trialBased_MI(data_name, excluded_channels, ...
                                         patient_number, excludeIED, fs)

    Data = load(data_name);
    blocks = Data.data.taskBlock;
    numCh  = size(blocks(1).LFP_portion(1).rawSignal,1);

    % concatenate raw EEG & phase labels
    raw_signal = cell(numCh,1);
    taskPhase_all = [];
    for b = 1:numel(blocks)
        if mod(b,2)==0, continue; end % skip observation blokcs
        for p = 1:numel(blocks(b).LFP_portion)
            pct = blocks(b).LFP_portion(p);
            idx1 = pct.markIdsInLFPdata(1);
            idx2 = pct.markIdsInLFPdata(2);
            seg = pct.rawSignal(:, idx1:idx2);
            for ch = 1:numCh
                raw_signal{ch} = [raw_signal{ch}, seg(ch,:)];
            end
            taskPhase_all = [taskPhase_all, pct.taskPhasePerLFPtimepoint];
        end
    end

    % find memory & movement trial‐index 
    mem_segs = find_segments(taskPhase_all, 2);
    mov_segs = find_segments(taskPhase_all, 1);

    % IED mask
    IEDmask_all = false(numCh, numel(taskPhase_all));
    if excludeIED
        for ch = 1:numCh
            IEDmask_all(ch,:) = removeIED(raw_signal{ch}, fs);
        end
    end

    % PAC parms
    LFrange = [1 3];    LFwidth = 2;  LFstep = 0.5;
    HFrange = [35 45];  HFwidth = 4;  HFstep = 2;
    nPhaseBins = 18;
    nShuffles = 500;

    edgeSamples = 50;
    minTrialLength = 3 * 1500; % for eggfilt

    % init output
    results.navigation = struct();

    % loop channels
    for ch = 1:numCh
        subjKey = ['subj',num2str(patient_number)];
        if isfield(excluded_channels,subjKey) && ...
           ismember(ch, excluded_channels.(subjKey))
            continue;
        end

        chanKey = ['channel_',num2str(ch)];

        F = {'memory_MI_trials','memory_thetaPower_trials','memory_trial_lengths', ...
             'movement_MI_trials','movement_thetaPower_trials','movement_trial_lengths'};
        for f = F, results.navigation.(chanKey).(f{1}) = []; end

        fullSig = raw_signal{ch};
        fullMask = ~IEDmask_all(ch,:);

        % memory trials
        for i = 1:numel(mem_segs)

            segIdx = mem_segs{i};

            if numel(segIdx) <= 2 * edgeSamples, continue; end
            segIdx = segIdx(edgeSamples+1:end-edgeSamples);

            if numel(segIdx) < minTrialLength, continue; end

            % compute MI
            zMImap = calculatePAC_mst_trial( ...
               fullSig, fullMask, segIdx, fs, ...
               LFrange, LFwidth, LFstep, ...
               HFrange, HFwidth, HFstep, ...
               nPhaseBins, nShuffles);

            MI_val = mean(zMImap.MI(:),'omitnan');

            % theta‐power
            thetaFilt = eegfilt(fullSig(segIdx), fs, 1, 3);
            thetaPower = abs(hilbert(thetaFilt)).^2;
            thetaMask = fullMask(segIdx);
            avg_theta_power = mean(thetaPower(thetaMask),'omitnan');

            % store
            results.navigation.(chanKey).memory_MI_trials(end+1) = MI_val;
            results.navigation.(chanKey).memory_thetaPower_trials(end+1) = avg_theta_power;
            results.navigation.(chanKey).memory_trial_lengths(end+1) = numel(segIdx);
        end

        % movement trials (optional)
        movementControl = false;
        if movementControl
            for i = 1:numel(mov_segs)

                segIdx = mov_segs{i};

                if numel(segIdx)<=2 * edgeSamples, continue; end
                segIdx = segIdx(edgeSamples + 1 : end - edgeSamples);
                if numel(segIdx)<minTrialLength, continue; end

                zMImap = calculatePAC_mst_trial( ...
                   fullSig, fullMask, segIdx, fs, ...
                   LFrange, LFwidth, LFstep, ...
                   HFrange, HFwidth, HFstep, ...
                   nPhaseBins, nShuffles);
                
                MI_val = mean(zMImap.MI(:),'omitnan');

                thetaFilt = eegfilt(fullSig(segIdx), fs, 1, 3);
                thetaPower = abs(hilbert(thetaFilt)).^2;
                thetaMask = fullMask(segIdx);
                avg_theta_power = mean(thetaPower(thetaMask),'omitnan');

                results.navigation.(chanKey).movement_MI_trials(end+1) = MI_val;
                results.navigation.(chanKey).movement_thetaPower_trials(end+1) = avg_theta_power;
                results.navigation.(chanKey).movement_trial_lengths(end+1) = numel(segIdx);
            end
        end
    end
end




%% 

function IED_mask = removeIED(raw_signal, fs)

    [nChannels, nTimepoints] = size(raw_signal);
    IED_mask = false(nChannels, nTimepoints);

    smoothingWindowSec = 0.05;
    detectionThreshold = 0.01; 
    exclusionWindowSec = 0.5;
    exclusionSamples = round(exclusionWindowSec * fs);

    for ch = 1:nChannels
        channelData = raw_signal(ch, :);

        rawIEDmask = detect_IED(fs, channelData);

        smoothedMask = smooth_gaussian_1d(double(rawIEDmask), fs, smoothingWindowSec);
        thresholdedMask= smoothedMask > detectionThreshold;

        flaggedIndices = find(thresholdedMask);

        for i = flaggedIndices
            lowerBound = max(1, i - exclusionSamples);
            upperBound = min(nTimepoints, i + exclusionSamples);
            thresholdedMask(lowerBound:upperBound) = true;
        end

        IED_mask(ch, :) = thresholdedMask;
    end
end

