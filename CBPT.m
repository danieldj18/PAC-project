function [ROI_map, t_map] = CBPT(PAC_results, nPermutations, alpha)

% Cluster based permutation test (channel level + z-scored)

    %% Extract Real and Surrogate MI Per Channel
    allFields = fieldnames(PAC_results);
    subjectNames = allFields(startsWith(allFields, 'subj'));
    channelCount = 0;

    for subj = 1:length(subjectNames)
        channels = fieldnames(PAC_results.(subjectNames{subj}));
        for ch = 1:length(channels)
            PAC_struct = PAC_results.(subjectNames{subj}).(channels{ch});
            if channelCount == 0
                LF_centers = PAC_struct.LF_frequencyCenters;
                HF_centers = PAC_struct.HF_frequencyCenters;
                [nLF, nHF] = size(PAC_struct.MI);
                nShuffles = size(PAC_struct.MIsurr, 3);
                real_MI_all = nan(nLF, nHF, 0);
                surrogate_MI_all = nan(nLF, nHF, nShuffles, 0);
            end

            channelCount = channelCount + 1;
            
            % Extract MI and Surrogate MI for each channel
            real_MI_all(:, :, channelCount) = PAC_struct.MI;
            surrogate_MI_all(:, :, :, channelCount) = PAC_struct.MIsurr;
        end
    end

    [nLF, nHF, nChannels] = size(real_MI_all);

    %% Z-Score Based On Surrogates For Each Channel Separately
    z_scored_MI_all = nan(nLF, nHF, nChannels);

    for ch = 1:nChannels
        % Compute mean and std for each channel independently
        surrogate_mean = mean(squeeze(surrogate_MI_all(:, :, :, ch)), 3);
        surrogate_std = std(squeeze(surrogate_MI_all(:, :, :, ch)), 0, 3);
        
        % Compute Z-score for each channel
        z_scored_MI_all(:, :, ch) = (real_MI_all(:, :, ch) - surrogate_mean) ./ surrogate_std;
    end

    %% Compute Mean & Std of Z-Scored MI Across Channels
    real_MI_mean = mean(z_scored_MI_all, 3);
    real_MI_std = std(z_scored_MI_all, 0, 3);

    %% Compute T-Statistic Across Channels
    t_map = real_MI_mean ./ (real_MI_std ./ sqrt(nChannels));

    %% Apply T-Threshold
    df = nChannels - 1;
    t_thresh = tinv(1 - alpha, df);
    binary_map = t_map > t_thresh;

    %% Identify Clusters
    CC = bwconncomp(binary_map, 4);
    cluster_sums_real = cellfun(@(idx) sum(t_map(idx)), CC.PixelIdxList);

    %% Permutation Testing (Cluster-Based Correction)
    nPermutations = max([1, nPermutations]);
    null_distribution = zeros(nPermutations, 1);
    
    for permIdx = 1:nPermutations
        permuted_real = zeros(nLF, nHF, nChannels);
        for ch = 1:nChannels
            sign_flip = randi([0 1]) * 2 - 1;
            permuted_real(:, :, ch) = sign_flip * z_scored_MI_all(:, :, ch);
        end
        perm_mean = mean(permuted_real, 3);
        perm_std  = std(permuted_real, 0, 3);
        perm_t_map = perm_mean ./ (perm_std ./ sqrt(nChannels));
        perm_binary_map = perm_t_map > t_thresh;
        CC_perm = bwconncomp(perm_binary_map, 4);
        perm_cluster_sums = cellfun(@(idx) sum(perm_t_map(idx)), CC_perm.PixelIdxList);
        null_distribution(permIdx) = max([0, perm_cluster_sums]);
    end

    %% Apply Cluster Threshold Based on Null Distribution
    critical_cluster_sum = prctile(null_distribution, 100 * (1 - alpha));
    final_map = zeros(size(binary_map));
    for i = 1:length(cluster_sums_real)
        if cluster_sums_real(i) >= critical_cluster_sum
            final_map(CC.PixelIdxList{i}) = 1;
        end
    end

    %% Dynamic MI Strength Filtering (optional; alternative filter method)
    percentile_threshold = 75; 

    % Collect MI values for all significant clusters
    cluster_MI_values = [];
    for i = 1:length(CC.PixelIdxList)
        if sum(final_map(CC.PixelIdxList{i})) > 0 
            cluster_MI_values = [cluster_MI_values; mean(real_MI_mean(CC.PixelIdxList{i}))];
        end
    end

    % Compute percentile-based threshold
    if ~isempty(cluster_MI_values)
        min_MI_threshold = prctile(cluster_MI_values, percentile_threshold);
    else
        min_MI_threshold = 0;
    end

    % Apply the MI strength filter
    final_map_filtered = zeros(size(final_map));
    for i = 1:length(CC.PixelIdxList)
        if sum(final_map(CC.PixelIdxList{i})) > 0
            mean_MI_in_cluster = mean(real_MI_mean(CC.PixelIdxList{i}));
            if mean_MI_in_cluster >= min_MI_threshold
                final_map_filtered(CC.PixelIdxList{i}) = 1;
            end
        end
    end

    % Clean up
    final_map_filtered = bwmorph(final_map_filtered, 'spur', Inf);
    final_map_filtered = imopen(final_map_filtered, strel('line', 4, 90));
    final_map_filtered = imopen(final_map_filtered, strel('line', 4, 0));
    final_map_filtered = imfill(final_map_filtered, 'holes');
    final_map_filtered = bwmorph(final_map_filtered, 'clean');

    %% Smoothing (optional) + Plot Results

    smoothMap = true;
    smoothFactor = 7;

    % For plotting, real_MI_mean is size [nLF x nHF].
    % Need to transpose for imagesc so that X=LF, Y=HF, and Z=real_MI_mean'.
    % If smoothing, also need to interpolate both the frequencies and the MI matrix.
    
    if smoothMap

        HF_centers_smooth = interpn(HF_centers, smoothFactor);
        LF_centers_smooth = interpn(LF_centers, smoothFactor);

        [LFgrid, HFgrid] = meshgrid(LF_centers, HF_centers);
        [LFgrid_smooth, HFgrid_smooth] = meshgrid(LF_centers_smooth, HF_centers_smooth);

        real_MI_mean_smooth = interp2(LFgrid, HFgrid, real_MI_mean', LFgrid_smooth, HFgrid_smooth, 'linear');

        figure;
        imagesc(LF_centers_smooth, HF_centers_smooth, real_MI_mean_smooth);
        set(gca, 'YDir', 'normal');
        hold on;

    else

        % Plot the unsmoothed version
        figure;
        imagesc(LF_centers, HF_centers, real_MI_mean');
        set(gca, 'YDir', 'normal');
        hold on;
    end

    % Plotting parameters
    xticks(LF_centers(1:2:end));
    xticklabels(string(LF_centers(1:2:end)));
    yticks(HF_centers(1:2:end));
    yticklabels(string(HF_centers(1:2:end)));

    xlabel('LF Frequency [Hz]');
    ylabel('HF Frequency [Hz]');

    % Need to upscale the final_map with the up-sampled grid, otherwise the 
    % boundaries won't perfectly align
    [B, ~] = bwboundaries(final_map_filtered', 'noholes');
    for k = 1:length(B)
        boundary = B{k};

        plot(LF_centers(boundary(:, 2)), HF_centers(boundary(:, 1)), 'k', 'LineWidth', 2);
    end

    c = colorbar;
    caxis([-0.5 3]);
    title(sprintf('Memory-Guided Navigation'));
    set(gcf, 'color', 'w');
    colormap(jet);
    hold on;

    % Plot dotted black lines at 35 and 45 Hz starting from 0 (prior
    % literature's gamma range)
    x_limits = xlim;
    y_limits = ylim;
    plot(x_limits, [35, 35], '--r', 'LineWidth', 1);
    plot(x_limits, [45, 45], '--r', 'LineWidth', 1);
    hold off;


    ROI_map = final_map_filtered;

end

