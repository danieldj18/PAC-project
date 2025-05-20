clear all
close all


%% Cluster-Based Permutation Test (single condition)

load('PAC_Memory.mat', 'PAC_results')

nPermutations = 5000;
alpha = 0.05;

zscore_map = 1; % raw or z-scored MI

if zscore_map
    
    [ROI_map, t_map] = CBPT(PAC_results, nPermutations, alpha);
else
    ROI_map = ClusterBased_PermutationTest_single_ChLevel(PAC_results, nPermutations, alpha);
    %[ROI_map, t_map] = ClusterBased_PermutationTest_single_SubjLevel(PAC_results, nPermutations, alpha);
end


% T-Statistic Plot
LF_range = [1 3];  % theta range

LF_frequencyCenters = PAC_results.subj1.channel1.LF_frequencyCenters;
HF_frequencyCenters = PAC_results.subj1.channel1.HF_frequencyCenters;

plot_TStatistic_vs_Gamma(t_map, LF_frequencyCenters, HF_frequencyCenters, LF_range);


%% Two-conditions comparison stats

ROI_mask = [];

filename = 'test123.mat';
%maxNumCompThreads(feature('numcores')); %if use all CPU cores
nDownsamples = 500;

% Define the number of subjects to process
patient_number = 5;

% Excluded channels
excluded_channels = struct();
excluded_channels.subj1 = [4]; 
excluded_channels.subj2 = [];
excluded_channels.subj3 = [2,4]; 
excluded_channels.subj4 = [];
excluded_channels.subj5 = [4]; 

Fsample = 250;
excludeIED = true;
calculate_MI = true;


%{

State Definitions:
-------------------
1  - Memory            | 2  - Movement
3  - Waiting           | 4  - Boundary
5  - Inner             | 6  - Good Memory
7  - Bad Memory        | 8  - Into Boundary
9  - Out of Boundary   | 21 - Memory NoStanding

Memory + Confound Splits:
-------------------------
Movement Split:  11 - Low Movement   | 12 - High Movement
Turning Split:   13 - Low Turn       | 14 - High Turn
Pupil Split:     15 - Slow Pupil     | 16 - Fast Pupil

Visually Guided Movement Split:   
 ------------------------------                   
Movement Split:  17 - Low Movement   | 18 - High Movement
Turning Split:   22 - Low Turn       | 23 - High Turn
Pupil Split:     24 - Slow Pupil     | 25 - High Pupil

Memory Learning Effects Split:   
 ------------------------------                   
                 19 - Memory Before  | 20 - Memory After

%}

condition_states = [1 2];
doSurrogate = true;

moving_threshold = 0.2; % standing vs moving
detour_threshold = 10000; % exclude detour (ratio: actual/perfect) traveling to hidden target 

LFrange = [1 3];   
LFwidth = 2;      
LF_stepsize = 0.5; 

HFrange = [35 45];
HFwidth = 4;     
HF_stepsize = 2; 


% Run analysis
process_PAC_analysis(patient_number, condition_states, excluded_channels, ...
    Fsample, excludeIED, moving_threshold, detour_threshold, LFrange, LFwidth, ...
    LF_stepsize, HFrange, HFwidth, HF_stepsize, nDownsamples, calculate_MI, filename, ROI_mask, doSurrogate)


%% Optional (use a new ROI mask)

% Define frequency parameters
% LFrange = [1 12];        
% LF_stepsize = 0.5;   
% 
% HFrange = [30 90];          
% HF_stepsize = 2;    

% LFrange = [1 6];        
% LF_stepsize = 0.5;   
% 
% HFrange = [30 80];          
% HF_stepsize = 2;   
% 
% % Define new ROI region
% LF_ROI = [1 3];   % ROI range for LF
% HF_ROI = [35 45]; % ROI range for HF

% Generate ROI mask
% ROI_mask = build_ROI_mask(LFrange, LF_stepsize, ...
%                           HFrange, HF_stepsize, ...
%                           LF_ROI, HF_ROI);


ROI_mask = ones(5,6);

newfile = 'newPAC.mat';
update_PAC_ROI(filename, ROI_mask, newfile);



%% Run statistics

close all

load(filename, 'PAC_results')

condition_names = {'Condition1', 'Condition2'};

% Run analysis
analyze_MI_statistics(PAC_results, condition_names);


%% analyze_MI_statistics()

function analyze_MI_statistics(PAC_results, condition_names)

MI_cond1_all = [];
MI_cond2_all = [];

% Subject names
subjectNames = fieldnames(PAC_results);
numSubjects = length(subjectNames);

% Figure for per-subject subplots
figure('Name', 'Per-Subject MI Comparisons', 'Color', 'w');
tiledlayout(1, numSubjects, 'Padding', 'compact', 'TileSpacing', 'compact');

for subjIdx = 1:numSubjects
    subjName = subjectNames{subjIdx};
    channelNames = fieldnames(PAC_results.(subjName));

    MI_cond1_subj = [];
    MI_cond2_subj = [];

    for chIdx = 1:length(channelNames)
        channelName = channelNames{chIdx};

        % Extract MI 
        if isfield(PAC_results.(subjName).(channelName), 'ROI_MI_cond1')
            MI1 = PAC_results.(subjName).(channelName).ROI_MI_cond1;
            MI2 = PAC_results.(subjName).(channelName).ROI_MI_cond2;
        else
            MI1 = PAC_results.(subjName).(channelName).real_MI_cond1;
            MI2 = PAC_results.(subjName).(channelName).real_MI_cond2;
        end

        MI_cond1_all = [MI_cond1_all; MI1];
        MI_cond2_all = [MI_cond2_all; MI2];

        MI_cond1_subj = [MI_cond1_subj; MI1];
        MI_cond2_subj = [MI_cond2_subj; MI2];
    end

    % Per subj plot
    nexttile;
    plot_MI_statistics(MI_cond1_subj, MI_cond2_subj, condition_names, subjName, false, true, true);
end

% All subj plot
figure('Name', 'Overall MI Comparison', 'Color', 'w');
plot_MI_statistics(MI_cond1_all, MI_cond2_all, condition_names, 'All Subjects & Channels', true, true, true);

end


%% Helper Function: plot_MI_statistics
function plot_MI_statistics(MI_cond1, MI_cond2, condition_names, group_label, do_stats, plot_lines, do_permutation, num_permutations)

    if nargin < 7
        do_permutation = false; 
    end
    if nargin < 8
        num_permutations = 10000; 
    end

    all_MI = [MI_cond1; MI_cond2];
    mean_all = mean(all_MI);
    std_all  = std(all_MI);
    MI_cond1_z = (MI_cond1 - mean_all) / std_all;
    MI_cond2_z = (MI_cond2 - mean_all) / std_all;

    mean_MI_z = [mean(MI_cond1_z), mean(MI_cond2_z)];

    hold on;
    bar(1:2, mean_MI_z, 0.6, 'FaceColor', [0.7, 0.7, 0.7]);

    % Scatter plot
    jitter_amount = 0.35;
    x1 = 1 + (rand(size(MI_cond1_z))-0.5)*jitter_amount;
    x2 = 2 + (rand(size(MI_cond2_z))-0.5)*jitter_amount;
    scatter(x1, MI_cond1_z, 50, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
    scatter(x2, MI_cond2_z, 60, 'r', 'filled', 'MarkerFaceAlpha', 0.6);

    % Connecting lines for paired data
    if plot_lines
        for i = 1:length(MI_cond1_z)
            plot([x1(i), x2(i)], [MI_cond1_z(i), MI_cond2_z(i)], 'k-', ...
                 'LineWidth', 0.5, 'Color', [0.5, 0.5, 0.5]);
        end
    end

    set(gca, 'XTick', 1:2, 'XTickLabel', condition_names, 'FontSize', 12);
    ylabel('Z-Scored MI');
    title(sprintf('%s', group_label), 'Interpreter', 'none');

    % Statistical testing
    if do_stats
        alpha = 0.05;
        diff_vals = MI_cond1 - MI_cond2;

        if do_permutation

            % Paired permutation test
            observed_diff = mean(MI_cond1) - mean(MI_cond2);
            
            perm_diffs = zeros(num_permutations,1);
            n = length(diff_vals);
            
            for p = 1:num_permutations

                flip_mask = (rand(n,1) > 0.5);
                
                perm_cond1 = MI_cond1;
                perm_cond2 = MI_cond2;
                
                perm_cond1(flip_mask) = MI_cond2(flip_mask);
                perm_cond2(flip_mask) = MI_cond1(flip_mask);
                
                perm_diffs(p) = mean(perm_cond1) - mean(perm_cond2);
            end
            
            % Two-sided p-value
            p_value = mean(abs(perm_diffs) >= abs(observed_diff));
            
            % Effect size (Cohen's d):
            effect_size = observed_diff / std(diff_vals);

            test_type = sprintf('Paired Permutation (%d perms)', num_permutations);

        else % (optional, not used in reality)
            % t-test (plus Shapiro-Wilk)
            [hDiff, ~] = swtest(diff_vals, alpha);  % from Stats Toolbox
            if hDiff == 0
                % Paired t-test
                test_type = 'Paired t-test';
                [~, p_value] = ttest(MI_cond1, MI_cond2);
                effect_size = mean(diff_vals) / std(diff_vals);
            else
                % Wilcoxon sign-rank
                test_type = 'Wilcoxon sign-rank';
                [p_value, ~, stats] = signrank(MI_cond1, MI_cond2);
                effect_size = abs(stats.zval) / sqrt(length(diff_vals));
            end
        end

        % Error bars (SEM)
        sem_MI_z = [std(MI_cond1_z)/sqrt(length(MI_cond1_z)), ...
                    std(MI_cond2_z)/sqrt(length(MI_cond2_z))];
        errorbar(1:2, mean_MI_z, sem_MI_z, 'k', 'LineStyle', 'none', 'LineWidth', 2);

        % Display significance if p < alpha
        if p_value < alpha
            text(1.5, max(mean_MI_z)+0.2, '*', 'FontSize', 18, ...
                 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end

        annotation('textbox', [0.65, 0.75, 0.3, 0.15], ...
                   'String', sprintf('%s\np = %.4f\nEffect size = %.3f', ...
                                     test_type, p_value, effect_size), ...
                   'EdgeColor', 'none');
    end

    hold off;
end



%% process_PAC_analysis()

function PAC_results = process_PAC_analysis(patient_number, condition_states, excluded_channels, ...
    Fsample, excludeIED, moving_threshold, detour_threshold, LFrange, LFwidth, ...
    LF_stepsize, HFrange, HFwidth, HF_stepsize, nDownsamples, calculate_MI, filename, ROI_mask, doSurrogate)

    % Initialize result structure
    PAC_results = struct();

    % Loop over each subject to process their data
    for subj = 1:patient_number
        
        data_name = sprintf('data_preprocessed_subj%d.mat', subj);
        
        % Process the subject's EEG signals
        [eegSignal, samplesToUse, ~, ~] = process_eeg_signal(data_name, excluded_channels, ...
            subj, Fsample, excludeIED, moving_threshold, detour_threshold);

        % Initialize structure for samplesToUse conditions
        samplesToUse_conditions = struct();

        % Loop through conditions
        for condIdx = 1:length(condition_states)
            state = condition_states(condIdx);
            samplesToUse_conditions.(['condition', num2str(condIdx)]) = select_samples(state, samplesToUse);
        end

        samplesToUse1 = samplesToUse_conditions.condition1;
        samplesToUse2 = samplesToUse_conditions.condition2;

        % Compute PAC if enabled
        if calculate_MI
            num_channels = size(eegSignal, 1);

            for channel_id = 1:num_channels
                if any(isnan(eegSignal(channel_id, :)))
                    continue; % Skip if signal has NaN values
                end

                % Extract current channel data
                eegSignal_currentChannel = eegSignal(channel_id, :);
                samplesToUse_currentChannel_condition1 = samplesToUse1(channel_id, :);
                samplesToUse_currentChannel_condition2 = samplesToUse2(channel_id, :);

                % Calculate PAC for both conditions
                PAC = calculate_ROI_MI_2Conditions(eegSignal_currentChannel, Fsample, ...
                    LFrange, LFwidth, LF_stepsize, HFrange, HFwidth, HF_stepsize, ...
                    samplesToUse_currentChannel_condition1, samplesToUse_currentChannel_condition2, nDownsamples, ROI_mask, doSurrogate);

                % Save PAC results
                PAC_results.(['subj', num2str(subj)]).(['channel', num2str(channel_id)]) = PAC;
            end
        end
        
        disp(['Processed PAC for subject ', num2str(subj)]);
    end

    save(filename, 'PAC_results');
end




function samplesToUse_condition = select_samples(state, samplesToUse)

    switch state
        case 1, samplesToUse_condition = samplesToUse.memory;
        case 2, samplesToUse_condition = samplesToUse.movement;
        case 3, samplesToUse_condition = samplesToUse.waiting;
        case 4, samplesToUse_condition = samplesToUse.boundary;
        case 5, samplesToUse_condition = samplesToUse.inner;
        case 6, samplesToUse_condition = samplesToUse.goodmemory;
        case 7, samplesToUse_condition = samplesToUse.badmemory;
        case 8, samplesToUse_condition = samplesToUse.IntoBoundary;
        case 9, samplesToUse_condition = samplesToUse.OutofBoundary;
        case 11, samplesToUse_condition = samplesToUse.memory_lowMovement;
        case 12, samplesToUse_condition = samplesToUse.memory_highMovement;
        case 13, samplesToUse_condition = samplesToUse.lowTurn;
        case 14, samplesToUse_condition = samplesToUse.highTurn;
        case 15, samplesToUse_condition = samplesToUse.slowPupil;
        case 16, samplesToUse_condition = samplesToUse.fastPupil;
        case 17, samplesToUse_condition = samplesToUse.visual_lowMovement;
        case 18, samplesToUse_condition = samplesToUse.visual_highMovement;
        case 19, samplesToUse_condition = samplesToUse.memory_before;
        case 20, samplesToUse_condition = samplesToUse.memory_after;
        case 21, samplesToUse_condition = samplesToUse.memory_NoStanding;
        case 22, samplesToUse_condition = samplesToUse.visual_lowTurn;
        case 23, samplesToUse_condition = samplesToUse.visual_highTurn;
        case 24, samplesToUse_condition = samplesToUse.visual_slow_pupil;
        case 25, samplesToUse_condition = samplesToUse.visual_fast_pupil;
        otherwise, error('Invalid state selection. Please use a valid state number.');
    end
end





%% process_eeg_signal(): a different version comparing to the main code

function [eegSignal_chXt, samplesToUse_chXt_combined, movement_data_amount, waiting_data_amount] = ...
   process_eeg_signal(data_name, excluded_channels, patient_number, fs, excludeIED, moving_threshold, detour_threshold)
    
    % Load the data within the function
    data = load(data_name);

    % Determine the number of channels from the data
    num_channels = size(data.data.taskBlock(1).LFP_portion(1).rawSignal, 1);

    % Preallocate large matrices for eegSignal and the corresponding mask
    eegSignal_chXt = nan(num_channels, 10000000);  
    samplesToUse_chXt = nan(num_channels, 10000000); 

    % Preallocate variables to concatenate taskPhase, xPosition, yPosition (for each subject)
    taskPhase_all = [];
    xPosition_all = [];
    yPosition_all = [];
    speed_all = [];
    xPupil_all = [];
    yPupil_all = [];
    movementDir_all = [];
    
    totalSigLen = 0;  

    % Loop over all task blocks
    for block_number = 1:length(data.data.taskBlock)
    
        skipLastBlock = false;
        if skipLastBlock
            if ismember(patient_number, [2, 3, 4, 5])
                even_blocks = 2:2:length(data.data.taskBlock);
                if block_number == even_blocks(end)
                    continue; 
                end
            end     
        else
            if mod(block_number, 2) == 0
                continue;
            end
        end

        % Loop over each portion within the block
        total_portions = length(data.data.taskBlock(block_number).LFP_portion);
    
        for LFP_portion_number = 1:total_portions
    
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
                
                % Initialize the mask (valid samples)
                currentSamplesToUse = ones(1, currentSigLen);
                
                %% Exclude IEDs
                if excludeIED
                    IED_samples = detect_IED(fs, raw_signal); % detect IEDs
                    
                    % Smoothing parameters
                    smoothing_win = 0.05; % seconds
                    threshold = 0.01; % for detecting IED
                    
                    % Smooth the IED_samples logical array
                    smoothed_bad_signal = smooth_gaussian_1d(IED_samples, fs, smoothing_win);
                    
                    % Threshold the smoothed signal
                    IED_samples = smoothed_bad_signal > threshold;
                    
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
    
    [logical_indices_memory, memory_good_logical, memory_bad_logical, logical_indices_navigation, logical_indices_waiting, ~] = ...
        find_memoryPerform_movement_indices(taskPhase_all, xPosition_all, yPosition_all, detour_threshold);
    

    % Memory learning effects split
    nbins = 3; % tercile
    memory_bins_logical = find_memory_learning_indices(taskPhase_all, xPosition_all, yPosition_all, detour_threshold, nbins);

    memory_before_logical = memory_bins_logical{1};
    memory_after_logical = memory_bins_logical{3};

    % Further seperate taskphase 2,3 into moving standing logical indices
    [logical_indices_movement, logical_indices_standing] = ...
        find_moving_standing_indices(logical_indices_navigation, logical_indices_waiting, speed_all, moving_threshold);

    % Check if moving/standing indices from above have approx same leng
    % data
%     movement_data_amount = sum(logical_indices_movement);
%     standing_data_amount = sum(logical_indices_standing);   
    
    % Find boundary vs. inner logical indices
    stratified = true;
    [boundary_logical, inner_logical] = find_boundary_inner_indices(xPosition_all, yPosition_all, stratified);

    % Find Into & Outof Boundary logical inidces
    [IntoBoundary_logical, OutofBoundary_logical] = find_InOut_boundary_indices(logical_indices_memory, logical_indices_navigation, logical_indices_waiting, boundary_logical);


    % ----------------------------------------------------------------------
    % Indices For Confounding Variables 
    % ----------------------------------------------------------------------

    % Pupil movement
    [fast_pupil, slow_pupil, ~] = find_pupil_speed_indices(xPupil_all, yPupil_all, ...
                                                     logical_indices_memory);

    [fast_pupil_movement, slow_pupil_movement, ~] = find_pupil_speed_indices(xPupil_all, yPupil_all, ...
                                                     logical_indices_movement);

    % Movement speed + turning 
    [logical_indices_memory_lowMovement, logical_indices_memory_highMovement,...
        logical_indices_lowTurn, logical_indices_highTurn, ...
        logical_indices_memory_standing, logical_indices_memory_moving, ...
        ~, ~] = ...
        find_lowHigh_movement_indices(logical_indices_memory, speed_all, movementDir_all);


    [logical_indices_movement_lowMovement, logical_indices_movement_highMovement,...
        logical_indices_movement_lowTurn, logical_indices_movement_highTurn, ...
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
    
    memory_standing_mask_replicated = repmat(logical_indices_memory_standing, size(samplesToUse_chXt, 1), 1); 
    memory_moving_mask_replicated   = repmat(logical_indices_memory_moving, size(samplesToUse_chXt, 1), 1);  
    
    lowTurn_mask_replicated  = repmat(logical_indices_lowTurn, size(samplesToUse_chXt, 1), 1);
    highTurn_mask_replicated = repmat(logical_indices_highTurn, size(samplesToUse_chXt, 1), 1);

    fastPupil_mask_replicated = repmat(fast_pupil, size(samplesToUse_chXt, 1), 1);
    slowPupil_mask_replicated = repmat(slow_pupil, size(samplesToUse_chXt, 1), 1);

    % Replicate movement speed median split for visually guided navigation
    visual_lowMovement_mask_replicated  = repmat(logical_indices_movement_lowMovement, size(samplesToUse_chXt, 1), 1);
    visual_highMovement_mask_replicated = repmat(logical_indices_movement_highMovement, size(samplesToUse_chXt, 1), 1);

    visual_lowTurn_mask_replicated  = repmat(logical_indices_movement_lowTurn, size(samplesToUse_chXt, 1), 1);
    visual_highTurn_mask_replicated = repmat(logical_indices_movement_highTurn, size(samplesToUse_chXt, 1), 1);

    visual_slow_pupil_mask_replicated = repmat(slow_pupil_movement, size(samplesToUse_chXt, 1), 1);
    visual_fast_pupil_mask_replicated = repmat(fast_pupil_movement, size(samplesToUse_chXt, 1), 1);

    
    % Replicate Memory Learning Effects 
    memory_before_replicated = repmat(memory_before_logical, size(samplesToUse_chXt, 1), 1);  
    memory_after_replicated  = repmat(memory_after_logical, size(samplesToUse_chXt, 1), 1);   

    
    % =============== Apply logical AND for all masks on non-NaN values =====================
    
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
        ~memory_standing_mask_replicated(non_nan_mask);  
    
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
 
    combined_visual_lowTurn_mask = nan(size(samplesToUse_chXt)); 
    combined_visual_lowTurn_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & visual_lowTurn_mask_replicated(non_nan_mask);
    
    combined_visual_highTurn_mask = nan(size(samplesToUse_chXt)); 
    combined_visual_highTurn_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & visual_highTurn_mask_replicated(non_nan_mask);

    combined_visual_slow_pupil_mask = nan(size(samplesToUse_chXt)); 
    combined_visual_slow_pupil_mask(non_nan_mask) = ...
        samplesToUse_chXt(non_nan_mask) & visual_slow_pupil_mask_replicated(non_nan_mask);

    combined_visual_fast_pupil_mask = nan(size(samplesToUse_chXt)); 
    combined_visual_fast_pupil_mask(non_nan_mask) = ...
        samplesToUse_chXt(non_nan_mask) & visual_fast_pupil_mask_replicated(non_nan_mask);
    % Memory Learning Effects
    combined_memory_before_mask = nan(size(samplesToUse_chXt));  
    combined_memory_before_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & memory_before_replicated(non_nan_mask);
    
    combined_memory_after_mask = nan(size(samplesToUse_chXt));   
    combined_memory_after_mask(non_nan_mask) = samplesToUse_chXt(non_nan_mask) & memory_after_replicated(non_nan_mask);

    
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
        'visual_highMovement', combined_visual_highMovement_mask, ...
        'visual_lowTurn', combined_visual_lowTurn_mask, ...
        'visual_highTurn', combined_visual_highTurn_mask, ...
        'visual_slow_pupil', combined_visual_slow_pupil_mask, ...
        'visual_fast_pupil', combined_visual_fast_pupil_mask, ...
        'memory_before', combined_memory_before_mask, ...   
        'memory_after', combined_memory_after_mask ...       
    );
    
    movement_data_amount = sum(combined_movement_mask(1,:));
    waiting_data_amount  = sum(combined_waiting_mask(1,:));
    memory_standing_data_amount  = sum(combined_memory_standing_mask(1,:));
    memory_moving_data_amount  = sum(combined_memory_moving_mask(1,:));

end












