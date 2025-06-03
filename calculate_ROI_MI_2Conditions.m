function PAC = calculate_ROI_MI_2Conditions(signal, ...
    samplingFrequency, LF_selectRange, LFwidth, LF_stepsize, ...
    HF_selectRange, HFwidth, HF_stepsize, ...
    samplesToUse_condition1, samplesToUse_condition2, nDownsamples, ROI_mask, doSurrogate)

    % Set default for optional inputs
    if nargin < 13
        doSurrogate = false;
    end

    if nargin < 12
        ROI_mask = [];
    end

    % Convert to logical indices (if not already)
    samplesToUse_condition1 = logical(samplesToUse_condition1);
    samplesToUse_condition2 = logical(samplesToUse_condition2);

    % Settings
    nPhaseBins = 18; % For calcModulationIndex_Tort2010
    segmentDuration = 1;
    samplesPerSegment = samplingFrequency * segmentDuration;

    % Compute LF & HF frequency centers
    LF_frequencyCenters = LF_selectRange(1):LF_stepsize:LF_selectRange(2);
    HF_frequencyCenters = HF_selectRange(1):HF_stepsize:HF_selectRange(2);

    % Adjust ROI mask size if provided
    if ~isempty(ROI_mask)
        ROI_mask = ROI_mask(1:length(LF_frequencyCenters), 1:length(HF_frequencyCenters));
    end

    % Initialize output struct
    PAC.dimensionInfo = '2D: [LF x HF], with real_MI per condition';
    PAC.LF_frequencyCenters = LF_frequencyCenters;
    PAC.HF_frequencyCenters = HF_frequencyCenters;
    PAC.real_MI_cond1 = nan(length(LF_frequencyCenters), length(HF_frequencyCenters));
    PAC.real_MI_cond2 = nan(length(LF_frequencyCenters), length(HF_frequencyCenters));
    PAC.downsampled_MI_cond1 = [];
    PAC.downsampled_MI_cond2 = [];

    % Check for trial length mismatch
    lengthCond1 = sum(samplesToUse_condition1);
    lengthCond2 = sum(samplesToUse_condition2);
    
    threshold = 0.05; % 5% threshold for trial length difference
    meanLength = (lengthCond1 + lengthCond2) / 2;
    downsampleNeeded = abs(lengthCond1 - lengthCond2) > (threshold * meanLength);

    % Determine which condition needs to be downsampled
    if downsampleNeeded
        if lengthCond1 > lengthCond2
            downsampleTarget = lengthCond2;
            downsampleCondition = 1; 
        else
            downsampleTarget = lengthCond1;
            downsampleCondition = 2; 
        end
    else
        downsampleCondition = 0;
    end

    totalSteps = length(LF_frequencyCenters) * length(HF_frequencyCenters);
    

    % Pre-allocate temporary results for PAC values for each LF-HF pair
    real_MI_cond1_temp = nan(length(LF_frequencyCenters), length(HF_frequencyCenters));
    real_MI_cond2_temp = nan(length(LF_frequencyCenters), length(HF_frequencyCenters));

    % Loop over LF_frequencyCenters
    parfor iLF = 1:length(LF_frequencyCenters)
        currentLFwindow = [LF_frequencyCenters(iLF) - LFwidth/2, LF_frequencyCenters(iLF) + LFwidth/2];

        % Preallocate results for this LF center
        temp_MI_cond1 = nan(1, length(HF_frequencyCenters));
        temp_MI_cond2 = nan(1, length(HF_frequencyCenters));
        
        for iHF = 1:length(HF_frequencyCenters)
            currentHFwindow = [HF_frequencyCenters(iHF) - HFwidth/2, HF_frequencyCenters(iHF) + HFwidth/2];
            
            disp(['PAC: ' num2str(LF_frequencyCenters(iLF)) ' Hz with ' num2str(HF_frequencyCenters(iHF)) ...
                  ' Hz (step ' num2str((iLF-1)*length(HF_frequencyCenters)+iHF) ' of ' num2str(totalSteps) ')']);

            % Filter signal at LF for current window
            LFsignal = eegfilt(signal, samplingFrequency, currentLFwindow(1), currentLFwindow(2));
            LF_phase = angle(hilbert(zscore(LFsignal)));

            % Filter signal at HF for current window
            HFsignal = eegfilt(signal, samplingFrequency, currentHFwindow(1), currentHFwindow(2));
            HF_amplitude = abs(hilbert(zscore(HFsignal)));

            % Compute MI for both conditions
            if downsampleNeeded
                if downsampleCondition == 1
                    % For condition 1, use downsampled surrogate MI
                    MI_stack = perform_downsampling(LF_phase, HF_amplitude, samplesToUse_condition1, ...
                        downsampleTarget, nDownsamples, nPhaseBins, doSurrogate, samplingFrequency);
                    temp_MI_cond1(iHF) = mean(MI_stack);
                    % For condition 2, use full-data surrogate MI
                    temp_MI_cond2(iHF) = perform_surrogate_calculation(LF_phase, HF_amplitude, samplesToUse_condition2, nPhaseBins, doSurrogate, samplingFrequency);
                else
                    MI_stack = perform_downsampling(LF_phase, HF_amplitude, samplesToUse_condition2, ...
                        downsampleTarget, nDownsamples, nPhaseBins, doSurrogate, samplingFrequency);
                    temp_MI_cond2(iHF) = mean(MI_stack);
                    temp_MI_cond1(iHF) = perform_surrogate_calculation(LF_phase, HF_amplitude, samplesToUse_condition1, nPhaseBins, doSurrogate, samplingFrequency);
                end
            else
                % No downsampling: use surrogate procedure for both conditions
                temp_MI_cond1(iHF) = perform_surrogate_calculation(LF_phase, HF_amplitude, samplesToUse_condition1, nPhaseBins, doSurrogate, samplingFrequency);
                temp_MI_cond2(iHF) = perform_surrogate_calculation(LF_phase, HF_amplitude, samplesToUse_condition2, nPhaseBins, doSurrogate, samplingFrequency);
            end
        end  % end for iHF
        
        % Store results for this LF center
        real_MI_cond1_temp(iLF, :) = temp_MI_cond1;
        real_MI_cond2_temp(iLF, :) = temp_MI_cond2;
    end  % end parfor iLF

    % Assign the temporary arrays to PAC
    PAC.real_MI_cond1 = real_MI_cond1_temp;
    PAC.real_MI_cond2 = real_MI_cond2_temp;

    % ROI MI Extraction if mask provided
    if ~isempty(ROI_mask)
        PAC.ROI_MI_cond1 = mean(PAC.real_MI_cond1(ROI_mask), 'all');
        PAC.ROI_MI_cond2 = mean(PAC.real_MI_cond2(ROI_mask), 'all');
    end
end

%% Function to Perform Downsampling + Surrogate MI Calculation

function MI_stack = perform_downsampling(LF_phase, HF_amplitude, samplesToUse, ...
    downsampleTarget, nDownsamples, nPhaseBins, doSurrogate, samplingFrequency)

    samplesPerSegment = samplingFrequency * 1;

    validIndices = find(samplesToUse);

    % Identify segment breaks
    segmentBreaks = find(diff(validIndices) > 1);
    segmentStartIndices = [validIndices(1), validIndices(segmentBreaks + 1)];
    segmentEndIndices   = [validIndices(segmentBreaks), validIndices(end)];
    segmentLengths      = segmentEndIndices - segmentStartIndices + 1;

    % Precompute all valid non-overlapping 1-second segments
    possibleSegments = {};
    for i = 1:length(segmentStartIndices)
        startIdx = segmentStartIndices(i);
        numFullSegments = floor(segmentLengths(i) / samplesPerSegment);
        for j = 0:numFullSegments-1
            possibleSegments{end+1} = startIdx + (j * samplesPerSegment) : (startIdx + (j+1) * samplesPerSegment - 1);
        end
    end

    totalAvailableSegments = length(possibleSegments);
    numSegmentsToSample = min(totalAvailableSegments, round(downsampleTarget / samplesPerSegment));

    MI_stack = zeros(nDownsamples, 1);
    nSurrogates = 200;

    for dIdx = 1:nDownsamples

        % Randomly select segments (without replacement)
        selectedSegmentIndices = randsample(1:totalAvailableSegments, numSegmentsToSample, false);
        selectedSegments = possibleSegments(selectedSegmentIndices);
        finalSelectedSamples = horzcat(selectedSegments{:});

        % Compute observed MI
        MI_obs = calcModulationIndex_Tort2010( ...
            LF_phase(finalSelectedSamples), HF_amplitude(finalSelectedSamples), nPhaseBins);

        if doSurrogate
            phase_vec = LF_phase(finalSelectedSamples);
            amp_vec   = HF_amplitude(finalSelectedSamples);

            nSeg = ceil(length(phase_vec) / samplesPerSegment);
            nFill = nSeg*samplesPerSegment - length(phase_vec);
            phase_padded = [phase_vec; nan(nFill,1)];

            phase_segments = reshape(phase_padded, samplesPerSegment, nSeg);

            surrogateMI = zeros(nSurrogates, 1);
            parfor s = 1:nSurrogates
                permIdx = randperm(nSeg);
                phase_shuffled = phase_segments(:,permIdx);
                phase_shuffled = phase_shuffled(:);
                phase_shuffled = phase_shuffled(~isnan(phase_shuffled));
                surrogateMI(s) = calcModulationIndex_Tort2010(phase_shuffled, amp_vec, nPhaseBins);
            end

            surrogateMean = mean(surrogateMI);
            surrogateStd  = std(surrogateMI);

            if surrogateStd == 0
                zMI = 0;
            else
                zMI = (MI_obs - surrogateMean) / surrogateStd;
            end
            MI_stack(dIdx) = zMI;

        else
            MI_stack(dIdx) = MI_obs;
        end
    end
end


%% Function to Perform Surrogate MI Calculation on Full (Non-Downsampled) Data

function zMI = perform_surrogate_calculation( ...
        LF_phase, HF_amplitude, samplesToUse, nPhaseBins, doSurrogate, samplingFrequency)

    phase_vec = LF_phase(samplesToUse).'; % row
    amp_vec = HF_amplitude(samplesToUse).'; % row
    N = numel(phase_vec);

    MI_obs = calcModulationIndex_Tort2010(phase_vec, amp_vec, nPhaseBins);

    if ~doSurrogate
        zMI = MI_obs;
        return
    end

    nSamplesPerBlock = samplingFrequency * 1;      
    nBlocks = ceil(N / nSamplesPerBlock);
    nPad = nBlocks * nSamplesPerBlock - N;
    phase_vec(end+1:end+nPad) = NaN;    

    phase_mat = reshape(phase_vec, nSamplesPerBlock, nBlocks);

    nSur = 200;
    MI_surr = zeros(1, nSur);

    parfor s = 1:nSur
        p_shuf = phase_mat(:, randperm(nBlocks)); % permute
        p_shuf = p_shuf(:); % flatten
        p_shuf = p_shuf(~isnan(p_shuf));  % drop padded NaNs
        MI_surr(s) = calcModulationIndex_Tort2010(p_shuf, amp_vec, nPhaseBins);
    end

    muS = mean(MI_surr);
    sdS = std(MI_surr);
    zMI = (MI_obs - muS) / sdS;
end

