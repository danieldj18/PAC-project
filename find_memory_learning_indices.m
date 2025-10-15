function memory_bins_logical = ...
    find_memory_learning_indices(taskPhase, xPosition, yPosition, detour_threshold, nBins)

% FIND_MEMORY_LEARNING_INDICES
%
% Splits valid memory segments into N bins for investigating potential learning effects,
% based on a global data-driven split (chronological).
%
% INPUTS:
%   taskPhase, xPosition, yPosition, detour_threshold : same usage as before
%   nBins (int) - how many bins to split into:
%        2 = median split (before/after)
%        3 = tercile, 4 = quartile, etc.
%
% OUTPUTS:
%   memory_bins_logical : cell array of length=nBins; each cell is a logical vector 
%                         marking where the data belongs to that bin
%   movement_logical    : [1 x nSamples] indicates movement (phase=1)
%   waiting_logical     : [1 x nSamples] indicates waiting  (phase=3)
%
% EXAMPLE:
%   % For a tercile split:
%   [mem_log, mov_log, wait_log] = find_memory_learning_indices(..., 1000, 3);
%   % Then mem_log{1}, mem_log{2}, mem_log{3} are each [1 x nSamples].


    % Handle input / setup
    if nargin < 4
        detour_threshold = 1000; % large => effectively no exclusion
    end
    
    if nargin < 5
        nBins = 2; % default = 2 (median split)
    end

    % Replace NaN in positions with linear interpolation
    if any(isnan(xPosition))
        xPosition = fillmissing(xPosition, 'linear');
    end
    if any(isnan(yPosition))
        yPosition = fillmissing(yPosition, 'linear');
    end

    nSamples = length(taskPhase);

    % Allocate
    movement_logical = false(1, nSamples);
    waiting_logical  = false(1, nSamples);

    % Instead of "before" and "after," create a cell array
    % of logical vectors, one per bin. 
    memory_bins_logical = cell(1, nBins);
    for b = 1:nBins
        memory_bins_logical{b} = false(1, nSamples);
    end


    % 1) Identify phase segments
    segments_movement = find_segments(taskPhase, 1);  % Movement
    segments_memory   = find_segments(taskPhase, 2);  % Memory
    segments_waiting  = find_segments(taskPhase, 3);  % Waiting

    % Mark movement
    for i = 1:length(segments_movement)
        seg_idx = segments_movement{i};
        movement_logical(seg_idx) = true;
    end

    % Mark waiting
    for i = 1:length(segments_waiting)
        seg_idx = segments_waiting{i};
        waiting_logical(seg_idx) = true;
    end


    % 2) Collect memory segments
    memory_segments = {};
    for i = 1:length(segments_memory)
        seg_idx = segments_memory{i};
        if i == 1
            % If you want to skip the first memory segment:
            continue;
        end
        if isempty(seg_idx), continue; end

        % Enforce max length of 10,000 from the END
        seg_end     = seg_idx(end);
        seg_len     = length(seg_idx);
        window_size = min(10000, seg_len);
        seg_start   = seg_end - window_size + 1;

        memory_segments{end+1} = seg_start:seg_end;
    end


    % 3) Filter memory segments by detour ratio
    final_memory_indices = [];
    for i = 1:length(memory_segments)
        seg_idx = memory_segments{i};
        if length(seg_idx) < 2, continue; end

        dx = diff(xPosition(seg_idx));
        dy = diff(yPosition(seg_idx));
        actual_route_distance = nansum(sqrt(dx.^2 + dy.^2));
        perfect_route_distance = sqrt( ...
            (xPosition(seg_idx(end)) - xPosition(seg_idx(1)))^2 + ...
            (yPosition(seg_idx(end)) - yPosition(seg_idx(1)))^2 );

        if actual_route_distance == 0
            detour_ratio = 1;
        else
            detour_ratio = actual_route_distance / perfect_route_distance;
        end

        if detour_ratio <= detour_threshold
            final_memory_indices = [final_memory_indices, seg_idx];
        end
    end

    final_memory_indices = sort(unique(final_memory_indices));
    total_mem_samples = length(final_memory_indices);
    if total_mem_samples < 1
        % No valid memory data
        return;
    end


    % 4) Split final_memory_indices into N bins
    % e.g. if nBins=2 -> median split
    % if nBins=3 -> terciles, etc.
    binEdges = round(linspace(1, total_mem_samples+1, nBins+1));

    for b = 1:nBins
        bin_startIdx = binEdges(b);
        bin_endIdx   = binEdges(b+1) - 1; 
        if bin_endIdx < bin_startIdx
            continue;
        end

        % Indices from final_memory_indices
        these_indices = final_memory_indices(bin_startIdx : bin_endIdx);

        % Mark them in the corresponding bin
        memory_bins_logical{b}(these_indices) = true;
    end

end

%% Sub function (find the condition indices)

function segments = find_segments(taskPhase, phase_of_interest)
    segments = {};
    in_segment = false;
    segment_start = 0;

    for i = 1:length(taskPhase)
        if taskPhase(i) == phase_of_interest && ~in_segment
            in_segment   = true;
            segment_start = i;
        elseif (taskPhase(i) ~= phase_of_interest || isnan(taskPhase(i))) && in_segment
            in_segment = false;
            segments{end+1} = segment_start:(i-1);
        end
    end
    if in_segment
        segments{end+1} = segment_start:length(taskPhase);
    end
end
