function [balanced_boundary_logical, balanced_inner_logical] = find_boundary_inner_indices(x_position_data, y_position_data,...
    use_stratified_sampling, window_size, boundary_dist)
  
    % Identifies boundary and inner indices with optional balanced stratified sampling
    % 
    % Parameters:
    % x_position_data - x-coordinates of positions
    % y_position_data - y-coordinates of positions
    % boundary_dist - threshold distance to classify boundary vs. inner
    % window_size - Number of consecutive samples per segment (for stratified sampling)
    % use_stratified_sampling - Boolean flag to enable/disable sampling
    %
    % Returns:
    % balanced_boundary_logical - Logical array for balanced boundary samples
    % balanced_inner_logical - Logical array for balanced inner samples

    if nargin < 3 || isempty(use_stratified_sampling)
        use_stratified_sampling = false;
    end
    if nargin < 4 || isempty(window_size)
        window_size = 250; 
    end
    if nargin < 5 || isempty(boundary_dist)
        boundary_dist = 1.2;
    end

    % Define boundary region limits
    x_max = max(x_position_data);
    x_min = min(x_position_data);
    y_max = max(y_position_data);
    y_min = min(y_position_data);

    y1 = y_max - boundary_dist;
    y0 = y_min + boundary_dist;
    x1 = x_max - boundary_dist;
    x0 = x_min + boundary_dist;

    % Classify positions
    inner_logical = (y_position_data > y0 & y_position_data < y1 & ...
                     x_position_data > x0 & x_position_data < x1);
    boundary_logical = ~inner_logical;

    if use_stratified_sampling
        boundary_count = sum(boundary_logical);
        inner_count = sum(inner_logical);

        % Find the smaller group size
        min_count = min(boundary_count, inner_count);

        % Keep all samples from the smaller group, only sample from the larger group
        if boundary_count > inner_count
            % Downsample boundary indices
            boundary_indices = find(boundary_logical);
            selected_boundary = select_segments(boundary_indices, min_count, window_size);
            balanced_boundary_logical = false(size(boundary_logical));
            balanced_boundary_logical(selected_boundary) = true;
            balanced_inner_logical = inner_logical; % Keep all inner samples
        elseif inner_count > boundary_count
            % Downsample inner indices
            inner_indices = find(inner_logical);
            selected_inner = select_segments(inner_indices, min_count, window_size);
            balanced_inner_logical = false(size(inner_logical));
            balanced_inner_logical(selected_inner) = true;
            balanced_boundary_logical = boundary_logical; % Keep all boundary samples
        else
            % Already balanced
            balanced_boundary_logical = boundary_logical;
            balanced_inner_logical = inner_logical;
        end

    else
        % Use original classifications without balancing
        balanced_boundary_logical = boundary_logical;
        balanced_inner_logical = inner_logical;
    end
end

%% func

function selected_indices = select_segments(indices, n_samples, window_size)

    selected_indices = [];   

    % Loop til the required number of samples is selected
    while length(selected_indices) < n_samples
        % Randomly sample the start index from the remaining indices
        start_idx = randsample(indices, 1);
        
        % Define the segment range
        end_idx = min(start_idx + window_size - 1, max(indices));
        new_segment = indices(indices >= start_idx & indices <= end_idx);

        % Ensure consistent dimensions for concatenation
        new_segment = new_segment(:); % need column vector
        selected_indices = selected_indices(:); % san check

        % Add the new segment to the selected indices
        selected_indices = [selected_indices; new_segment];

        % Remove the selected segment from the pool of indices
        indices(indices >= start_idx & indices <= end_idx) = [];

        % Stop if enough samples are collected
        if length(selected_indices) >= n_samples
            selected_indices = selected_indices(1:n_samples); % Trim excess samples
            break;
        end
    end
end
