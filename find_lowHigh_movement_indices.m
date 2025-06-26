function [logical_indices_lowMovement, logical_indices_highMovement, ...
          logical_indices_lowTurn, logical_indices_highTurn, ...
          logical_indices_veryLowMovement, logical_indices_restMovement, ...
          turning_speed_all, angular_acceleration_all] = ...
          find_lowHigh_movement_indices(memory_mask, speed_all,movementDir_all, blockSize, smoothingWindowSize)

    if nargin < 4 || isempty(blockSize)
        blockSize = 250;
    end

    if nargin < 5 || isempty(smoothingWindowSize)
        smoothingWindowSize = 0;
    end

    % san check dimen
    memory_mask = logical(memory_mask(:).');
    speed_all = speed_all(:).';
    movementDir_all = movementDir_all(:).';


    N = numel(speed_all);
    
    % optional smoothing
    if smoothingWindowSize > 1
        movementDir_smooth = movmean(movementDir_all, smoothingWindowSize); % smoothed
    else
        movementDir_smooth = movementDir_all; %
    end
    
    % turning-speed trace
    turning_speed_all = [0 abs(diff(movementDir_smooth))]; % deg/sample
    turning_speed_all(N) = turning_speed_all(end); % pad last sample
    angular_acceleration_all = [turning_speed_all(2) - turning_speed_all(1) diff(turning_speed_all)]; % angular velo
    angular_acceleration_all(N) = angular_acceleration_all(end); % pad last sample
    
    % block split
    totalBlocks = floor(N / blockSize);
    block_speed = nan(1, totalBlocks); % mean speed per block
    block_turn = nan(1, totalBlocks); % sum turning per block
    maskCount_blk = zeros(1, totalBlocks); % masked samples per block
    isMemBlock = false(1, totalBlocks); % flag memory blocks
    
    for b = 1:totalBlocks
        idx = (b - 1) * blockSize + (1:blockSize);
        block_speed(b) = mean(speed_all(idx), 'omitnan'); % mean speed
        block_turn(b) = sum(turning_speed_all(idx)); % total turning
        maskCount_blk(b) = sum(memory_mask(idx)); % masked count
        isMemBlock(b) = maskCount_blk(b) > 0; % memory flag
    end
    
    % movement-speed split (sample balanced)
    logical_indices_lowMovement = false(1, N);
    logical_indices_highMovement = false(1, N);
    memBlocks_speed = find(isMemBlock);
    
    if numel(memBlocks_speed) >= 2

        spdVals = block_speed(memBlocks_speed); % speed of mem blocks
        mCounts = maskCount_blk(memBlocks_speed); % samples per block

        [~, sortOrd] = sort(spdVals); % ascending speed
        sortedBlocks = memBlocks_speed(sortOrd); % ordered blocks
        sortedCounts = mCounts(sortOrd); % ordered counts

        cumCounts = cumsum(sortedCounts); % cumulative samples
        halfSamp = cumCounts(end) / 2; % 50 % point

        splitIdx = find(cumCounts >= halfSamp, 1, 'first'); % split index
        lowBlocks_speed = sortedBlocks(1:splitIdx); % slow blocks
        highBlocks_speed = sortedBlocks(splitIdx + 1:end); % fast blocks

        % build splited blocks back to original length mask
        for b = lowBlocks_speed
            idx = (b - 1) * blockSize + (1:blockSize); % block indices
            memIdx = idx(memory_mask(idx)); % masked samples
            logical_indices_lowMovement(memIdx) = true; % mark slow
        end

        for b = highBlocks_speed
            idx = (b - 1) * blockSize + (1:blockSize); % block indices
            memIdx = idx(memory_mask(idx)); % masked samples
            logical_indices_highMovement(memIdx) = true; % mark fast
        end
    end
    
    % very-low vs rest movement
    logical_indices_veryLowMovement = memory_mask & (speed_all < 0.2); % near zero speed
    logical_indices_restMovement = memory_mask & (speed_all >= 0.2); % all others
    
    % turning-speed split (sample balanced)
    logical_indices_lowTurn = false(1, N);
    logical_indices_highTurn = false(1, N);
    memBlocks_turn = find(isMemBlock);
    
    if numel(memBlocks_turn) >= 2
        
        turnVals = block_turn(memBlocks_turn); % turning values
        tCounts = maskCount_blk(memBlocks_turn); % samples per block

        [~, sortOrd] = sort(turnVals); % ascending turning
        sortedBlocks = memBlocks_turn(sortOrd); % ordered blocks
        sortedCounts = tCounts(sortOrd); % ordered counts

        cumCounts = cumsum(sortedCounts); % cumulative samples
        halfSamp = cumCounts(end) / 2; % 50% point

        splitIdx = find(cumCounts >= halfSamp, 1, 'first'); % split index
        lowBlocks_turn = sortedBlocks(1:splitIdx); % low turn
        highBlocks_turn = sortedBlocks(splitIdx + 1:end); % high turn

        % build splited blocks back to original length mask
        for b = lowBlocks_turn
            idx = (b - 1) * blockSize + (1:blockSize); % block indices
            memIdx = idx(memory_mask(idx)); % masked samples
            logical_indices_lowTurn(memIdx) = true; % mark low turn
        end

        for b = highBlocks_turn
            idx = (b - 1) * blockSize + (1:blockSize); % block indices
            memIdx = idx(memory_mask(idx)); % masked samples
            logical_indices_highTurn(memIdx) = true; % mark high turn
        end
    end
end
