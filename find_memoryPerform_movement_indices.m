function [memory_logical, memory_good_logical, memory_bad_logical, ...
          movement_logical, waiting_logical, detour_percentage] = ...
    find_memoryPerform_movement_indices(taskPhase, xPosition, yPosition, detour_threshold)

    % Parameters for local (per-segment) filter (optional)
    use_memory_segment_filter           = false;
    memory_segment_filter_type          = 'middle';    % choose bet/ 'first','last','middle'
    memory_segment_filter_duration_sec  = 2;          
    sample_rate                         = 250;         

    % Parameters for global split (optional)
    use_memory_global_split      = false;
    memory_global_split_type     = 'after'; % 'before' or 'after'
    memory_global_split_fraction = 0.5;     % e.g., 0.5 => half

    % Downsampling
    use_memory_downsampling = false; % toggle on/off

    % Memory performance input (optional)
    if nargin < 4
        detour_threshold = 10000; % very large => no exclusion
    end

    % Interpolate x/y position data
    if any(isnan(xPosition))
        xPosition = fillmissing(xPosition, 'linear');
    end

    if any(isnan(yPosition))
        yPosition = fillmissing(yPosition, 'linear');
    end

    nSamples = length(taskPhase);

    % Initialize outputs
    movement_logical     = zeros(1, nSamples);
    memory_logical       = zeros(1, nSamples);
    memory_good_logical  = zeros(1, nSamples);
    memory_bad_logical   = zeros(1, nSamples);
    waiting_logical      = zeros(1, nSamples);
    detour_percentage    = [];

    % Identify contiguous segments for each phase
    segments_phase_1 = find_segments(taskPhase, 1);  % Movement
    segments_phase_2 = find_segments(taskPhase, 2);  % Memory
    segments_phase_3 = find_segments(taskPhase, 3);  % Waiting

    % 1) Mark movment segments (phase 1)
    for i = 1:length(segments_phase_1)
        seg_idx = segments_phase_1{i};
        if ~isempty(seg_idx) && taskPhase(seg_idx(1)) == 1
            movement_logical(seg_idx) = 1;
        end
    end

    % 2) Get memory segments (phase 2) 
    raw_memory_segments = {};
    for k = 1:length(segments_phase_2)
        seg_idx = segments_phase_2{k};
        if isempty(seg_idx), continue; end

        if k == 1
            continue;
        end

        if taskPhase(seg_idx(1)) ~= 2
            continue;
        end

        seg_end     = seg_idx(end);
        max_window  = 100000; % large number to turn off
        window_size = min(max_window, length(seg_idx));
        seg_start   = seg_end - window_size + 1;
        seg_crop    = seg_start:seg_end;

        raw_memory_segments{end+1} = seg_crop;
    end

    % 3) Filter 'first','middle','last' (optional)
    filtered_memory_segments = {};
    for iSeg = 1:length(raw_memory_segments)
        seg_range = raw_memory_segments{iSeg};

        if use_memory_segment_filter
            remove_or_keepN = round(memory_segment_filter_duration_sec * sample_rate);

            switch lower(memory_segment_filter_type)
                case 'first'
                    keepN  = min(remove_or_keepN, length(seg_range));
                    seg_range = seg_range(1:keepN);

                case 'last'
                    keepN  = min(remove_or_keepN, length(seg_range));
                    seg_range = seg_range(end-keepN+1 : end);

                case 'middle'
                    removeN = remove_or_keepN; 
                    if length(seg_range) <= 2*removeN
                        seg_range = [];
                    else
                        seg_range = seg_range(removeN+1 : end-removeN);
                    end

                otherwise
                    warning('Unknown memory_segment_filter_type.');
            end
        end

        if ~isempty(seg_range)
            filtered_memory_segments{end+1} = seg_range; 
        end
    end

    % 4) Split into 'before' or 'after' portion (optional)
    final_memory_segments = {};
    if ~use_memory_global_split
        final_memory_segments = filtered_memory_segments;
    else
        % Gather all memory indices, then keep only the 'before' or 'after' 
        all_mem_idxs = [];
        for iSeg = 1:length(filtered_memory_segments)
            all_mem_idxs = [all_mem_idxs, filtered_memory_segments{iSeg}];
        end
        if isempty(all_mem_idxs)
            final_memory_segments = {};
        else
            all_mem_idxs = sort(unique(all_mem_idxs), 'ascend');
            total_mem    = length(all_mem_idxs);
            cutoff       = round(memory_global_split_fraction * total_mem);
            cutoff       = max(1, min(cutoff, total_mem)); 

            before_idxs  = all_mem_idxs(1 : cutoff);
            after_idxs   = all_mem_idxs(cutoff+1 : end);

            switch lower(memory_global_split_type)
                case 'before'
                    keep_idxs = before_idxs;
                case 'after'
                    keep_idxs = after_idxs;
                otherwise
                    warning('Unknown memory_global_split_type. Keeping all memory data.');
                    keep_idxs = all_mem_idxs;
            end

            keep_idxs = sort(keep_idxs, 'ascend');

            % Re-split into contiguous subranges
            subSegs = find_contiguous_subranges(keep_idxs);
            final_memory_segments = [final_memory_segments, subSegs]; 
        end
    end

    % 5) Detour calc + exclusion
    accepted_mem_indices = [];  
    seg_detr_vals        = []; 

    for iSeg = 1:length(final_memory_segments)
        seg_idx = final_memory_segments{iSeg};
        if isempty(seg_idx), continue; end

        seg_start = seg_idx(1);
        seg_end   = seg_idx(end);

        dx = diff(xPosition(seg_idx));
        dy = diff(yPosition(seg_idx));
        actual_route_distance = nansum(sqrt(dx.^2 + dy.^2));

        perfect_route_distance = sqrt( ...
            (xPosition(seg_end) - xPosition(seg_start))^2 + ...
            (yPosition(seg_end) - yPosition(seg_start))^2 );

        if actual_route_distance == 0
            this_detour = 1;
        else
            this_detour = actual_route_distance / perfect_route_distance;
        end

        if this_detour <= detour_threshold

            % Accept these indices
            accepted_mem_indices = [accepted_mem_indices, seg_idx]; 
            seg_detr_vals        = [seg_detr_vals, this_detour];   
        end
    end

    accepted_mem_indices = sort(unique(accepted_mem_indices));  

    % 5.1) Downsample memory to match movement data 
    if use_memory_downsampling && ~isempty(accepted_mem_indices)
        total_movement = sum(movement_logical);
        total_memory   = length(accepted_mem_indices);

        if total_memory > total_movement

            chunk_size = sample_rate; % 1 sec => 250 samples

            % Break accepted_mem_indices into contiguous subranges
            all_subranges = find_contiguous_subranges(accepted_mem_indices);

            % Break each contiguous subrange into 1s chunks
            possible_chunks = {}; % each entry is a 1s index range
            for r = 1:length(all_subranges)
                r_idx = all_subranges{r};
                len_r = length(r_idx);

                start_i = 1;
                while (start_i + chunk_size - 1) <= len_r
                    chunk_indices = r_idx(start_i : (start_i + chunk_size - 1));
                    possible_chunks{end+1} = chunk_indices;
                    start_i = start_i + chunk_size;
                end
            end

            % Shuffle the possible chunks randomly
            rand_order = randperm(length(possible_chunks));
            possible_chunks = possible_chunks(rand_order);

            % Pick chunks until match ~ total_movement
            target_needed = total_movement; 
            chosen_indices = [];
            n_collected    = 0;
            iC = 1;

            while iC <= length(possible_chunks) && n_collected < target_needed
                this_chunk = possible_chunks{iC};
                chosen_indices = [chosen_indices, this_chunk];
                n_collected   = n_collected + length(this_chunk);
                iC = iC + 1;
            end

            % Sort and keep unique
            chosen_indices = sort(unique(chosen_indices));

            % That is our new set of accepted memory indices
            accepted_mem_indices = chosen_indices;
        end
    end

    % 5.2) Mark memory_logical
    memory_logical(accepted_mem_indices) = 1;

    % For the final detour percentage array
    detour_percentage = seg_detr_vals;

    % Good vs Bad memory (tersile split)
    if ~isempty(detour_percentage)
        % Compute the lower and upper cutoffs based on a tercile split
        lower_cutoff = prctile(detour_percentage, 33);
        upper_cutoff = prctile(detour_percentage, 66);
        
        % Reinitialize memory
        memory_logical(:)      = 0;
        memory_good_logical(:) = 0;
        memory_bad_logical(:)  = 0;
        
        acceptedSegments_final = {}; 
        detourVals_final       = [];
        
        for iSeg = 1:length(final_memory_segments)
            seg_idx = final_memory_segments{iSeg};
            if isempty(seg_idx), continue; end

            dx = diff(xPosition(seg_idx));
            dy = diff(yPosition(seg_idx));
            actual_route_distance = nansum(sqrt(dx.^2 + dy.^2));

            seg_start = seg_idx(1);
            seg_end   = seg_idx(end);
            perfect_route_distance = sqrt( (xPosition(seg_end) - xPosition(seg_start))^2 + ...
                                           (yPosition(seg_end) - yPosition(seg_start))^2 );

            if actual_route_distance == 0
                this_detour = 1;
            else
                this_detour = actual_route_distance / perfect_route_distance;
            end

            if this_detour <= detour_threshold
  
                seg_intersect = intersect(seg_idx, accepted_mem_indices);
                if ~isempty(seg_intersect)
                    acceptedSegments_final{end+1} = seg_intersect; 
                    detourVals_final(end+1) = this_detour;      
                end
            end
        end

        % Perform the tercile split: 
        %   classify as GOOD if <= lower_cutoff
        %   classify as BAD if  >= upper_cutoff 
        if ~isempty(detourVals_final)
            for k = 1:length(acceptedSegments_final)
                seg_int = acceptedSegments_final{k};
                this_val = detourVals_final(k);

                if this_val <= lower_cutoff
                    memory_good_logical(seg_int) = 1;
                    memory_logical(seg_int) = 1;
                elseif this_val >= upper_cutoff
                    memory_bad_logical(seg_int) = 1;
                    memory_logical(seg_int) = 1;
                end
            end
        end

        % Update detour_percentage
        detour_percentage = detourVals_final;
    end

    % 7) Mark waiting segments (phase 3)
    for i = 1:length(segments_phase_3)
        seg_idx = segments_phase_3{i};
        if ~isempty(seg_idx) && taskPhase(seg_idx(1)) == 3
            waiting_logical(seg_idx) = 1;
        end
    end
end

%% SUB-FUNCTION: find_segments
% Identifies contiguous trial indices for each condition

function segments = find_segments(taskPhase, phase_of_interest)
    segments = {};
    in_segment = false;
    segment_start = 0;

    for i = 1:length(taskPhase)
        if taskPhase(i) == phase_of_interest && ~in_segment
            in_segment = true;
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

%% SUB-FUNCTION: find_contiguous_subranges
% Takes a sorted vector of indices (e.g. [10 11 12 20 21 22 23])
% and breaks it into contiguous subranges => {[10 11 12], [20 21 22 23]}

function subranges = find_contiguous_subranges(indexVec)
    subranges = {};
    if isempty(indexVec), return; end

    start_idx = indexVec(1);
    prev_val  = indexVec(1);

    for i = 2:length(indexVec)
        if indexVec(i) ~= prev_val + 1
            % Found a gap => close off previous range
            subranges{end+1} = start_idx:prev_val;
            start_idx = indexVec(i);
        end
        prev_val = indexVec(i);
    end

    % Add the final subrange
    subranges{end+1} = start_idx:prev_val;
end
