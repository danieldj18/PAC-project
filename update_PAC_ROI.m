function update_PAC_ROI(PAC_results_file, new_ROI_mask, output_filename)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function: update_PAC_ROI
%
% Purpose:
% - Updates `PAC.ROI_MI_cond1` and `PAC.ROI_MI_cond2` in an **already computed** PAC results file.
% - Uses a **new ROI mask** without recomputing PAC from raw EEG signals.
%
% Inputs:
% - PAC_results_file: Filename of the stored PAC results (MAT file).
% - new_ROI_mask: Binary mask (logical matrix) specifying the new ROI.
% - output_filename: Filename to save the updated PAC results.
%
% Key Steps:
% 1. Load the existing PAC results.
% 2. Loop through all subjects and channels.
% 3. Apply the new ROI mask to `PAC.real_MI_cond1` and `PAC.real_MI_cond2`.
% 4. Compute new mean MI values for the new ROI.
% 5. Replace old `PAC.ROI_MI_cond1` & `PAC.ROI_MI_cond2` with the new values.
% 6. Save the updated PAC results.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Load the existing PAC results
    load(PAC_results_file, 'PAC_results');

    % Extract all subject names
    subjectNames = fieldnames(PAC_results);

    % Loop through all subjects
    for subjIdx = 1:length(subjectNames)
        subjName = subjectNames{subjIdx};
        channelNames = fieldnames(PAC_results.(subjName));

        % Loop through all channels for this subject
        for chIdx = 1:length(channelNames)
            channelName = channelNames{chIdx};

            % Extract PAC data
            PAC = PAC_results.(subjName).(channelName);

            % Ensure ROI mask size matches MI matrix dimensions
            if size(new_ROI_mask) ~= size(PAC.real_MI_cond1)
                error('ROI mask size mismatch. Ensure it matches MI matrix dimensions.');
            end

            % Compute new ROI-based MI (mean over new ROI)
            PAC.ROI_MI_cond1 = mean(PAC.real_MI_cond1(new_ROI_mask), 'all');
            PAC.ROI_MI_cond2 = mean(PAC.real_MI_cond2(new_ROI_mask), 'all');

            % Store the updated PAC back in the results
            PAC_results.(subjName).(channelName) = PAC;
        end
    end

    % Save the updated results
    save(output_filename, 'PAC_results');
    
    disp(['Updated PAC results saved to ', output_filename]);

end
