function [fast_pupil, slow_pupil, pupil_speed_all] = find_pupil_speed_indices( ...
    xPupil_all, yPupil_all, logical_indices, window_size, blockSize)


    if nargin < 4 || isempty(window_size)
        window_size = 1; 
    end
    
    if nargin < 5 || isempty(blockSize)
        blockSize   = 250; 
    end
    
    xPupil_all = xPupil_all(:);
    yPupil_all = yPupil_all(:);
    logical_indices = logical(logical_indices(:));
    N = numel(xPupil_all);
    
    % Pupil speed
    pupil_speed = nan(N, 1);

    if window_size == 1
        pupil_speed(1:N-1) = hypot(diff(xPupil_all), diff(yPupil_all));
        pupil_speed(N) = pupil_speed(N-1);
    else
        for i = 1:(N - window_size)
            dx = xPupil_all(i + window_size) - xPupil_all(i);
            dy = yPupil_all(i + window_size) - yPupil_all(i);
            pupil_speed(i) = hypot(dx, dy) / window_size;
        end
        pupil_speed(N - window_size + 1:N) = pupil_speed(N - window_size);
    end

    pupil_speed_all = pupil_speed;
    
    % Block split
    totalBlocks = floor(N / blockSize);
    block_speed = nan(1, totalBlocks); % mean speed per block
    maskCount_blk = zeros(1, totalBlocks); % masked samples per block
    isMaskBlock = false(1, totalBlocks); % does block touch mask (mem/mov) or not
    
    for b = 1:totalBlocks
        idx = (b - 1) * blockSize + (1:blockSize);
        block_speed(b) = mean(pupil_speed(idx), 'omitnan');
        maskCount_blk(b) = sum(logical_indices(idx));
        isMaskBlock(b) = maskCount_blk(b) > 0;
    end
    
    % Balance splits
    slow_pupil = false(N, 1);
    fast_pupil = false(N, 1);
    
    maskBlocks = find(isMaskBlock);

    spdVals = block_speed(maskBlocks); % speeds of mask blocks
    mCounts = maskCount_blk(maskBlocks); % sample counts

    [~, sortOrd] = sort(spdVals); % ascending speed
    sortedBlocks = maskBlocks(sortOrd);
    sortedCounts = mCounts(sortOrd);

    cumCounts = cumsum(sortedCounts); % cumulative samples
    halfSamp = cumCounts(end) / 2; % 50% data amount point
    splitIdx = find(cumCounts >= halfSamp, 1, 'first');

    slowBlocks = sortedBlocks(1:splitIdx);
    fastBlocks = sortedBlocks(splitIdx + 1:end);

    % Assign blocks back to recreate original length mask for each split
    for b = slowBlocks
        idx = (b - 1) * blockSize + (1:blockSize);
        maskIdx = idx(logical_indices(idx));
        slow_pupil(maskIdx) = true;
    end

    for b = fastBlocks
        idx = (b - 1) * blockSize + (1:blockSize);
        maskIdx = idx(logical_indices(idx));
        fast_pupil(maskIdx) = true;

    end
end
