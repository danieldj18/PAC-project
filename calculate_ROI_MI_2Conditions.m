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

%% Function for Downsampling + Surrogate MI Calculation
function MI_stack = perform_downsampling(...
    LF_phase, HF_amplitude, samplesToUse, ...
    downsampleTarget, nDownsamples, nPhaseBins, doSurrogate, samplingFrequency)

    samplesPerSegment = samplingFrequency;

    % Find all valid data chunks
    validIndices      = find(samplesToUse);
    segmentBreaks     = find(diff(validIndices) > 1);
    segmentStartIdx   = [validIndices(1), validIndices(segmentBreaks+1)];
    segmentEndIdx     = [validIndices(segmentBreaks), validIndices(end)];
    possibleSegments  = {};
    for i = 1:length(segmentStartIdx)
        segLen = segmentEndIdx(i) - segmentStartIdx(i) + 1;
        nBlocks = floor(segLen / samplesPerSegment);
        for b = 0:nBlocks-1
            idx = segmentStartIdx(i) + b*samplesPerSegment : ...
                  segmentStartIdx(i)+(b+1)*samplesPerSegment-1;
            possibleSegments{end+1} = idx;
        end
    end

    totalBlocks = length(possibleSegments);
    nBlocksToSample = min(totalBlocks, round(downsampleTarget/samplesPerSegment));

    MI_stack = zeros(nDownsamples,1);
    nSurrogates = 500; % Make smaller for testing

    for dIdx = 1:nDownsamples
        % Pick random chunks to match target data length
        pickIdx      = randsample(totalBlocks, nBlocksToSample, false);
        selectedSegs = possibleSegments(pickIdx);
        finalSamples = horzcat(selectedSegs{:});

        % Observed MI
        MI_obs = calcModulationIndex_Tort2010( ...
            LF_phase(finalSamples), HF_amplitude(finalSamples), nPhaseBins);

        if doSurrogate
            surrogateMI = zeros(nSurrogates,1);
            L = length(finalSamples);

            for s = 1:nSurrogates
                if doSurrogate == 2
                    % Circshift surrogate (cross-check)
                    minOff = round(0.2*samplingFrequency);
                    offset = randi([minOff+1, L-1]);
                    amp_surr = circshift(HF_amplitude(finalSamples), offset);
                    surrogateMI(s) = calcModulationIndex_Tort2010( ...
                        LF_phase(finalSamples), amp_surr, nPhaseBins);
                else % Block shuffle surrogate (default)
                    permOrder  = randperm(numel(selectedSegs));
                    permSegs   = selectedSegs(permOrder);
                    permSamples= horzcat(permSegs{:});
                    surrogateMI(s) = calcModulationIndex_Tort2010( ...
                        LF_phase(permSamples), HF_amplitude(permSamples), nPhaseBins);
                end
            end

            mu = mean(surrogateMI);
            sigma = std(surrogateMI);
            if sigma == 0
                MI_stack(dIdx) = 0;
            else
                MI_stack(dIdx) = (MI_obs - mu)/sigma;
            end
        else
            % If no surrogate => just raw MI
            MI_stack(dIdx) = MI_obs;
        end
    end
end


%% Function to Perform Surrogate MI Calculation on Full (Non-Downsampled) Data
function zMI = perform_surrogate_calculation(LF_phase, HF_amplitude, samplesToUse, nPhaseBins, doSurrogate, samplingFrequency)
    
    samplesPerSegment = samplingFrequency; 

    validIndices = find(samplesToUse);
    segmentBreaks = find(diff(validIndices) > 1);
    segmentStartIndices = [validIndices(1), validIndices(segmentBreaks + 1)];
    segmentEndIndices   = [validIndices(segmentBreaks), validIndices(end)];
    segmentLengths      = segmentEndIndices - segmentStartIndices + 1;

    possibleSegments = {}; % Possible number of 1s segments
    for i = 1:length(segmentStartIndices)
        startIdx = segmentStartIndices(i);
        numFullSegments = floor(segmentLengths(i) / samplesPerSegment);
        for j = 0:numFullSegments-1
            possibleSegments{end+1} = startIdx + (j * samplesPerSegment) : (startIdx + (j+1) * samplesPerSegment - 1);
        end
    end

    totalAvailableSegments = length(possibleSegments);

    % For full-data surrogate, use all available segments
    numSegmentsToSample = totalAvailableSegments;
    finalSelectedSamples = horzcat(possibleSegments{:});
    
    % Compute observed MI
    MI_obs = calcModulationIndex_Tort2010( ...
        LF_phase(finalSelectedSamples), HF_amplitude(finalSelectedSamples), nPhaseBins);

    if doSurrogate
        nSurrogates = 500; % Make smaller for testing
        surrogateMI = zeros(nSurrogates, 1);
        minOffset = round(0.2 * samplingFrequency);
        L = length(finalSelectedSamples);

        for s = 1:nSurrogates
            if doSurrogate == 2 % Circshift surrogate (cross-check)
                offset = randi([minOffset+1, L-1]);
                HF_amplitude_sur = circshift(HF_amplitude(finalSelectedSamples), offset);
                
                surrogateMI(s) = calcModulationIndex_Tort2010( ...
                    LF_phase(finalSelectedSamples), HF_amplitude_sur, nPhaseBins);
            else % Block shuffle surrogate (default)
                permOrder = randperm(length(possibleSegments));
                permSegs  = possibleSegments(permOrder);
                permSamples = horzcat(permSegs{:});
                
                surrogateMI(s) = calcModulationIndex_Tort2010( ...
                    LF_phase(permSamples), HF_amplitude(permSamples), ...
                    nPhaseBins);
            end
        end
        surrogateMean = mean(surrogateMI);
        surrogateStd = std(surrogateMI);
        if surrogateStd == 0
            zMI = 0;
        else
            zMI = (MI_obs - surrogateMean) / surrogateStd;
        end
    else
        zMI = MI_obs;
    end
end

