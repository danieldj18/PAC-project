clear all
close all

% Save result to this name 
filename = 'PAC.mat';

% Define the number of subjects to process
subject_start = 1;
patient_number = 5;

% Excluded channels
excluded_channels = struct();
excluded_channels.subj1 = [4]; % Exclude channel 4 for subject 1
excluded_channels.subj2 = [];
excluded_channels.subj3 = [2,4]; % Exclude channels 2 and 4 for subject 3
excluded_channels.subj4 = [];
excluded_channels.subj5 = [4]; % Exclude channel 4 for subject 5

Fsample = 250;
excludeIED = true;
calculate_MI = true;

%{

State Definitions:
-------------------
1  - Memory            | 2   - Movement
3  - Waiting           | 4   - Boundary
5  - Inner             | 6   - Good Memory
7  - Bad Memory        | 8   - Into Boundary
9  - Out of Boundary   | 19  - Memory_NoStanding


Memory + Confound Splits:
-------------------------
Memory Movement Split:   11 - Low Movement   | 12 - High Movement
Memory Turning Split:    13 - Low Turn       | 14 - High Turn
Memory Pupil Split:      15 - Slow Pupil     | 16 - Fast Pupil

Visually Guided Movement Split:   
                    
                         20 - Low Movement   | 21 - High Movement
%}

state = 1;  % Adjust this value to select the analysis condition

printsteps = true; % Print the PAC processing steps statement or not

numShuffles = 500; % Number of shuffles for surrogate analysis
moving_threshold = 0.2; % threshold to separate standing vs moving (m/s)
detour_threshold = 100; % threshold to exclude detour (ratio: actual/perfect) traveling to hidden target 

% Define PAC calculation parameters
LFrange = [1 12];    % Low-frequency range (theta band), default [1 12]
LFwidth = 2;         % Frequency window for low frequency, default 2
LF_stepsize = 0.5;   % Low-frequency window step size (frequency centers), default 0.5

HFrange = [30 90];   % High-frequency range (gamma band), default [30 90]
HFwidth = 4;         % Frequency window for high frequency, default 4
HF_stepsize = 2;     % High-frequency window step size (frequency centers), default 2

% Initialize struct for storing PAC results
PAC_results = struct();

% Loop over each subject to process their data
for subj = subject_start:patient_number
    
    data_name = sprintf('data_preprocessed_subj%d.mat', subj);
    
    % Call the process_eeg_signal function to process the subject's EEG signals
    [eegSignal, samplesToUse, movement_data_amount, waiting_data_amount] = ...
        process_eeg_signal(data_name, excluded_channels, subj, Fsample, excludeIED, moving_threshold, detour_threshold);

    % ===== State Selection for Analysis =====
    % Memory(1) | Movement(2) | Waiting(3) | Boundary(4) | Inner(5) | Good Mem(6) | Bad Mem(7)
    % IntoBoundary(8) | OutofBoundary(9) | Mem-LowMove(11) | Mem-HighMove(12)
    % Low Turn(13) | High Turn(14) | Slow Pupil(15) | Fast Pupil(16)
    % Mem-Standing(17) | Mem-Moving(18) | Mem-NoStanding(19)  
    
    switch state
        case 1, samplesToUse = samplesToUse.memory;
        case 2, samplesToUse = samplesToUse.movement;
        case 3, samplesToUse = samplesToUse.waiting;

        case 4, samplesToUse = samplesToUse.boundary;
        case 5, samplesToUse = samplesToUse.inner;

        case 6, samplesToUse = samplesToUse.goodmemory;
        case 7, samplesToUse = samplesToUse.badmemory;

        case 8, samplesToUse = samplesToUse.IntoBoundary;
        case 9, samplesToUse = samplesToUse.OutofBoundary;

        case 11, samplesToUse = samplesToUse.memory_lowMovement;
        case 12, samplesToUse = samplesToUse.memory_highMovement;

        case 13, samplesToUse = samplesToUse.lowTurn;
        case 14, samplesToUse = samplesToUse.highTurn;

        case 15, samplesToUse = samplesToUse.slowPupil;
        case 16, samplesToUse = samplesToUse.fastPupil;

        case 17, samplesToUse = samplesToUse.memory_standing;
        case 18, samplesToUse = samplesToUse.memory_moving;

        case 19, samplesToUse = samplesToUse.memory_NoStanding;

        case 20, samplesToUse = samplesToUse.visual_lowMovement;
        case 21, samplesToUse = samplesToUse.visual_highMovement;
            
        otherwise, error('Invalid state selection. Use a valid state number.');
    end

    
    % Compute PAC for each channel independently
    num_channels = size(eegSignal, 1);  % Number of channels for the subject

    if calculate_MI

        for channel_id = 1:num_channels
        
            % Skip if the signal contains NaN values for the current channel
            if sum(isnan(eegSignal(channel_id, :))) ~= 0
               continue;  % Skip this channel
            end
    
            % Extract the current channel's EEG data and samplesToUse mask
            eegSignal_currentChannel = eegSignal(channel_id, :);
            samplesToUse_currentChannel = samplesToUse(channel_id, :);
            
            % Calculate PAC for the current channel
            PAC = calculatePAC_mst_fast(eegSignal_currentChannel, eegSignal_currentChannel, Fsample, ...
                LFrange, LFwidth, LF_stepsize, HFrange, HFwidth, HF_stepsize, numShuffles,...
                samplesToUse_currentChannel, printsteps);
            
            % Save the PAC results for the current subject and channel
            PAC_results.(['subj', num2str(subj)]).(['channel', num2str(channel_id)]) = PAC;
            
        end
    end
    
    disp(['Processed PAC for subject ', num2str(subj)]);

end

save(filename, 'PAC_results');

%% Compute p-values separately before adding to the structure (optional)

%filename = 'PAC_movement';


allCh_pval = compute_allCh_pvalues(PAC_results); % Standard Tort 
%allCh_pval = compute_allCh_pvalues_bootstrap(PAC_results, 10000); % Bootstrap Method
%allCh_pval = compute_allCh_pvalues_block_bootstrap(PAC_results, 10000);

% Add the computed p-values to the structure after calculation
PAC_results.allCh_pval = allCh_pval;

% Save the PAC results
save(filename, 'PAC_results');

%{

'PAC_results_allSubjects.mat' structure:

PAC_results (1x1 struct)
├── subj1
│   ├── channel1
│   │   ├── MI          % Modulation Index for this channel
│   │   ├── MI_pVal     % p-value for the MI
│   │   ├── MIsurrMean  % Mean MI from surrogate data
│   │   ├── MIsurrStd   % Standard deviation of MI from surrogate data
│   │   ├── meanHFampPerLFphaseBin  % Mean high-frequency amplitude per phase bin
│   │   └── nSamplesPerLFphaseBin   % Number of samples per phase bin
│   ├── channel2
│   │   ├── MI          % Same structure for channel 2
│   │   └── ...         % Other PAC metrics for channel 2
│   └── ...
├── subj2
│   └── ...
├── subj3
│   └── ...
└── allCh_pval

%}


%% Functions

function [eegSignal_chXt, samplesToUse_chXt_combined, movement_data_amount, waiting_data_amount] = ...
   process_eeg_signal(data_name, excluded_channels, patient_number, fs, excludeIED, moving_threshold, detour_threshold)
    
    % Load the data within the function
    data = load(data_name);

    % Determine the number of channels from the data
    num_channels = size(data.data.taskBlock(1).LFP_portion(1).rawSignal, 1);

    % Preallocate large matrices for eegSignal and the corresponding mask
    eegSignal_chXt = nan(num_channels, 10000000);  % Large size, later trimmed
    samplesToUse_chXt = nan(num_channels, 10000000);  % Same size for the mask

    % Preallocate variables to concatenate taskPhase, xPosition, yPosition (for each subject)
    taskPhase_all = [];
    xPosition_all = [];
    yPosition_all = [];
    speed_all = [];
    xPupil_all = [];
    yPupil_all = [];
    movementDir_all = [];
    
    totalSigLen = 0;  % Keep track of how much of the preallocated matrix is used

    % Loop over all task blocks
    for block_number = 1:length(data.data.taskBlock)
    
        % Skip observation task blocks (assuming these are even-numbered)
        if mod(block_number, 2) == 0
            continue;
        end
  
    
        % Loop over each portion within the block (once per block)
        total_portions = length(data.data.taskBlock(block_number).LFP_portion);
    
        for LFP_portion_number = 1:total_portions
    
            % Access the current portion
            portion = data.data.taskBlock(block_number).LFP_portion(LFP_portion_number);
    
            % Concatenate taskPhase, xPosition, and yPosition for this portion
            taskPhase_all = [taskPhase_all, portion.taskPhasePerLFPtimepoint];
            xPosition_all = [xPosition_all, portion.xPositionPerLFPtimepoint];
            yPosition_all = [yPosition_all, portion.yPositionPerLFPtimepoint];
            speed_all = [speed_all, portion.movSpeedPerLFPtimepoint];
            
            xPupil_all = [xPupil_all, portion.pupilXperLFPtimepoint];
            yPupil_all = [yPupil_all, portion.pupilYperLFPtimepoint];
            movementDir_all = [movementDir_all, portion.movDirPerLFPtimepoint];

    
            % Get the signal length for this portion
            currentSigLen = length(portion.rawSignal(1, portion.markIdsInLFPdata(1):portion.markIdsInLFPdata(2)));
    
            % Loop over all channels
            for channel_id = 1:num_channels
    
                % Skip excluded channels
                subject_key = ['subj', num2str(patient_number)];

                if isfield(excluded_channels, subject_key) && ismember(channel_id, excluded_channels.(subject_key))
                    continue;
                end
                
                % Get the raw signal for the current channel
                raw_signal = portion.rawSignal(channel_id, portion.markIdsInLFPdata(1):portion.markIdsInLFPdata(2));
                
                % Initialize the mask as all ones (valid samples)
                currentSamplesToUse = ones(1, currentSigLen);
                
                %% Exclude IEDs if necessary
                if excludeIED
                    IED_samples = detect_IED(fs, raw_signal); % Detect IEDs in the raw signal
                    
                    % Smoothing parameters
                    smoothing_win = 0.05; % seconds
                    threshold = 0.01; % for detecting IED
                    
                    % Smooth the IED_samples logical array
                    smoothed_bad_signal = smooth_gaussian_1d(IED_samples, fs, smoothing_win);
                    
                    % Threshold the smoothed signal
                    IED_samples = smoothed_bad_signal > threshold;

                    % Optional further IED window exclusion
                    excludeIEDWin = 1;
                    if excludeIEDWin
                        % Exclude a ±1.5 s window around each IED sample
                        excludeWin = round(0.5 * fs);  % number of samples in 1.5 s
                        nSamples   = length(IED_samples);
                    
                        % Loop over all flagged samples
                        flaggedIndices = find(IED_samples);
                        for idx = 1:numel(flaggedIndices)
                            iEvt       = flaggedIndices(idx); 
                            lowerBound = max(iEvt - excludeWin, 1);
                            upperBound = min(iEvt + excludeWin, nSamples);
                    
                            % Mark those samples as IEDs as well
                            IED_samples(lowerBound:upperBound) = true;
                        end
                    end
                    
                    % Mark IED samples as invalid (0)
                    currentSamplesToUse(IED_samples) = 0;
                end
    
                %% Exclude first and last second, because of mark artifact
                currentSamplesToUse(1:fs) = 0;
                currentSamplesToUse(end-fs+1:end) = 0;
    
                % Concatenate this portion of the signal and the mask into the main matrices
                eegSignal_chXt(channel_id, totalSigLen + 1 : totalSigLen + currentSigLen) = raw_signal;
                samplesToUse_chXt(channel_id, totalSigLen + 1 : totalSigLen + currentSigLen) = currentSamplesToUse;
            end
    
            % Update the total length of the signal processed so far
            totalSigLen = totalSigLen + currentSigLen;
        end
    end
    
    % Remove any unused preallocated space in the matrices
    eegSignal_chXt(:, totalSigLen + 1 : end) = [];
    samplesToUse_chXt(:, totalSigLen + 1 : end) = [];
    
    % Now calculate the memory, movement, and waiting logical indices using concatenated taskPhase, xPosition, and yPosition
    [logical_indices_memory, memory_good_logical, memory_bad_logical, logical_indices_navigation, logical_indices_waiting, ~] = ...
        find_memoryPerform_movement_indices(taskPhase_all, xPosition_all, yPosition_all, detour_threshold);
    

    % Further seperate taskphase 2,3 into moving standing logical indices +
    [logical_indices_movement, logical_indices_standing] = ...
        find_moving_standing_indices(logical_indices_navigation, logical_indices_waiting, speed_all, moving_threshold);

    % Check if moving/standing indices from above have approx same leng
    % data
%     movement_data_amount = sum(logical_indices_movement);
%     standing_data_amount = sum(logical_indices_standing);   
    
    % Find boundary vs. inner logical indices
    stratified = false;
    [boundary_logical, inner_logical] = find_boundary_inner_indices(xPosition_all, yPosition_all, stratified);

    % Find Into & Outof Boundary logical inidces
    [IntoBoundary_logical, OutofBoundary_logical] = find_InOut_boundary_indices(logical_indices_memory, logical_indices_navigation, logical_indices_waiting, boundary_logical);


    % ----------------------------------------------------------------------
    % Indices For Confounding Variables 
    % ----------------------------------------------------------------------

    % Pupil movement
    [fast_pupil, slow_pupil, ~] = find_pupil_speed_indices(xPupil_all, yPupil_all, ...
                                                     logical_indices_memory);

    % Movement speed + turning 
    [logical_indices_memory_lowMovement, logical_indices_memory_highMovement,...
        logical_indices_lowTurn, logical_indices_highTurn, ...
        logical_indices_memory_standing, logical_indices_memory_moving, ...
        ~, ~] = ...
        find_lowHigh_movement_indices(logical_indices_memory, speed_all, movementDir_all);
    

    [logical_indices_movement_lowMovement, logical_indices_movement_highMovement,...
        ~, ~, ...
        ~, ~, ...
        ~, ~] = find_lowHigh_movement_indices(logical_indices_navigation, speed_all, movementDir_all);


    % ======================================= Matching Shape ======================================
    
    % Replicate the memory, movement, and waiting masks to match the shape of samplesToUse_chXt
    memory_mask_replicated = repmat(logical_indices_memory, size(samplesToUse_chXt, 1), 1);  
    memory_good_replicated = repmat(memory_good_logical, size(samplesToUse_chXt, 1), 1);  
    memory_bad_replicated  = repmat(memory_bad_logical, size(samplesToUse_chXt, 1), 1);  
    
    movement_mask_replicated = repmat(logical_indices_movement, size(samplesToUse_chXt, 1), 1);  
    waiting_mask_replicated  = repmat(logical_indices_standing, size(samplesToUse_chXt, 1), 1);  
    
    % Replicate the boundary and inner masks
    boundary_mask_replicated = repmat(boundary_logical, size(samplesToUse_chXt, 1), 1);
    inner_mask_replicated    = repmat(inner_logical, size(samplesToUse_chXt, 1), 1);
    
    IntoBoundary_mask_replicated = repmat(IntoBoundary_logical, size(samplesToUse_chXt, 1), 1);
    OutofBoundary_mask_replicated = repmat(OutofBoundary_logical, size(samplesToUse_chXt, 1), 1);
    
    % Replicate movement, turn, pupil splits within memory
    memory_lowMovement_mask_replicated  = repmat(logical_indices_memory_lowMovement, size(samplesToUse_chXt, 1), 1);
    memory_highMovement_mask_replicated = repmat(logical_indices_memory_highMovement, size(samplesToUse_chXt, 1), 1);
    
    memory_standing_mask_replicated = repmat(logical_indices_memory_standing, size(samplesToUse_chXt, 1), 1); % NEW
    memory_moving_mask_replicated   = repmat(logical_indices_memory_moving, size(samplesToUse_chXt, 1), 1);   % NEW
    
    lowTurn_mask_replicated  = repmat(logical_indices_lowTurn, size(samplesToUse_chXt, 1), 1);
    highTurn_mask_replicated = repmat(logical_indices_highTurn, size(samplesToUse_chXt, 1), 1);

    fastPupil_mask_replicated = repmat(fast_pupil, size(samplesToUse_chXt, 1), 1);
    slowPupil_mask_replicated = repmat(slow_pupil, size(samplesToUse_chXt, 1), 1);

    % Replicate movement speed median split for visually guided navigation
    visual_lowMovement_mask_replicated  = repmat(logical_indices_movement_lowMovement, size(samplesToUse_chXt, 1), 1);
    visual_highMovement_mask_replicated = repmat(logical_indices_movement_highMovement, size(samplesToUse_chXt, 1), 1);
    
    
    % =======================================================================================
    % =============== Apply logical AND for all masks on non-NaN values =====================
    % =======================================================================================
    
    non_nan_mask = ~isnan(samplesToUse_chXt);
    
    % Memory-related masks
    combined_memory_mask = nan(size(samplesToUse_chXt));
    combined_memory_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & memory_mask_replicated(non_nan_mask);
    
    combined_good_memory_mask = nan(size(samplesToUse_chXt));
    combined_good_memory_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & memory_good_replicated(non_nan_mask);
    
    combined_bad_memory_mask = nan(size(samplesToUse_chXt));
    combined_bad_memory_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & memory_bad_replicated(non_nan_mask);
    
    % Movement and waiting
    combined_movement_mask = nan(size(samplesToUse_chXt));
    combined_movement_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & movement_mask_replicated(non_nan_mask);
    
    combined_waiting_mask = nan(size(samplesToUse_chXt));
    combined_waiting_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & waiting_mask_replicated(non_nan_mask);
    
    % Boundary and inner
    combined_boundary_mask = nan(size(samplesToUse_chXt));
    combined_boundary_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & boundary_mask_replicated(non_nan_mask);
    
    combined_inner_mask = nan(size(samplesToUse_chXt));
    combined_inner_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & inner_mask_replicated(non_nan_mask);
    
    combined_IntoBoundary_mask = nan(size(samplesToUse_chXt));
    combined_IntoBoundary_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & IntoBoundary_mask_replicated(non_nan_mask);
    
    combined_OutofBoundary_mask = nan(size(samplesToUse_chXt));
    combined_OutofBoundary_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & OutofBoundary_mask_replicated(non_nan_mask);
    
    % Memory-based standing and moving splits
    combined_memory_standing_mask = nan(size(samplesToUse_chXt)); % NEW
    combined_memory_standing_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & memory_standing_mask_replicated(non_nan_mask);
    
    combined_memory_moving_mask = nan(size(samplesToUse_chXt));   % NEW
    combined_memory_moving_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & memory_moving_mask_replicated(non_nan_mask);

    % Memory without standing (memory_NoStanding)
    combined_memory_NoStanding_mask = nan(size(samplesToUse_chXt));
    combined_memory_NoStanding_mask(non_nan_mask) = ...
        samplesToUse_chXt(non_nan_mask) & memory_mask_replicated(non_nan_mask) & ...
        ~memory_standing_mask_replicated(non_nan_mask);  % Exclude memory_standing
    
    % Existing memory movement splits
    combined_memory_lowMovement_mask = nan(size(samplesToUse_chXt));
    combined_memory_lowMovement_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & memory_lowMovement_mask_replicated(non_nan_mask);
    
    combined_memory_highMovement_mask = nan(size(samplesToUse_chXt));
    combined_memory_highMovement_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & memory_highMovement_mask_replicated(non_nan_mask);
    
    % Turning behavior
    combined_lowTurn_mask = nan(size(samplesToUse_chXt));
    combined_lowTurn_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & lowTurn_mask_replicated(non_nan_mask);
    
    combined_highTurn_mask = nan(size(samplesToUse_chXt));
    combined_highTurn_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & highTurn_mask_replicated(non_nan_mask);
    
    % Pupil movement
    combined_fastPupil_mask = nan(size(samplesToUse_chXt));
    combined_fastPupil_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & fastPupil_mask_replicated(non_nan_mask);
    
    combined_slowPupil_mask = nan(size(samplesToUse_chXt));
    combined_slowPupil_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & slowPupil_mask_replicated(non_nan_mask);

    % Visually Guided Navigation Movement Splits
    combined_visual_lowMovement_mask = nan(size(samplesToUse_chXt)); 
    combined_visual_lowMovement_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & visual_lowMovement_mask_replicated(non_nan_mask);

    combined_visual_highMovement_mask = nan(size(samplesToUse_chXt)); 
    combined_visual_highMovement_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & visual_highMovement_mask_replicated(non_nan_mask);
    
    
    % ======================== Store All Masks into Struct ========================
    samplesToUse_chXt_combined = struct(...
        'memory', combined_memory_mask, ...
        'goodmemory', combined_good_memory_mask, ...
        'badmemory', combined_bad_memory_mask, ...
        'movement', combined_movement_mask, ...
        'waiting', combined_waiting_mask, ...
        'boundary', combined_boundary_mask, ...
        'inner', combined_inner_mask, ...
        'IntoBoundary', combined_IntoBoundary_mask, ...
        'OutofBoundary', combined_OutofBoundary_mask, ...
        'memory_standing', combined_memory_standing_mask, ...   
        'memory_NoStanding', combined_memory_NoStanding_mask, ...   
        'memory_moving', combined_memory_moving_mask, ...      
        'memory_lowMovement', combined_memory_lowMovement_mask, ...
        'memory_highMovement', combined_memory_highMovement_mask, ...
        'lowTurn', combined_lowTurn_mask, ...
        'highTurn', combined_highTurn_mask, ...
        'fastPupil', combined_fastPupil_mask, ...
        'slowPupil', combined_slowPupil_mask, ...
        'visual_lowMovement', combined_visual_lowMovement_mask, ...
        'visual_highMovement', combined_visual_highMovement_mask ...
    );
    
    movement_data_amount = sum(combined_movement_mask(1,:));
    waiting_data_amount  = sum(combined_waiting_mask(1,:));
    memory_standing_data_amount  = sum(combined_memory_standing_mask(1,:));
    memory_moving_data_amount  = sum(combined_memory_moving_mask(1,:));

end





