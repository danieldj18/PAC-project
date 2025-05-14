function [fast_pupil, slow_pupil, pupil_speed_all] = find_pupil_speed_indices(xPupil_all, yPupil_all, ...
                                                     logical_indices, window_size)

% Set default value for window_size if not provided
if nargin < 4
    window_size = 1;
end

% Ensure column vectors
xPupil_all = xPupil_all(:);
yPupil_all = yPupil_all(:);
logical_indices = logical(logical_indices(:));


N = length(xPupil_all);
if length(yPupil_all) ~= N
    error('xPupil_all and yPupil_all must be the same length.');
end
if length(logical_indices) ~= N
    error('logical_indices must be the same length as xPupil_all.');
end
if window_size < 1 || window_size >= N
    error('window_size should be between 1 and (length of data - 1).');
end


% 1) Compute pupil speed (movement) depending on window_size
pupil_speed = nan(N, 1);

if window_size == 1
    % *Fast path* using diff for point-to-point
    dx = diff(xPupil_all);
    dy = diff(yPupil_all);
    speed_diff = sqrt(dx.^2 + dy.^2);  % size = N-1
    pupil_speed(1 : N-1) = speed_diff;
    pupil_speed(N) = pupil_speed(N-1);  % replicate last speed
else
    % General window_size approach
    for i = 1 : (N - window_size)
        dx = xPupil_all(i + window_size) - xPupil_all(i);
        dy = yPupil_all(i + window_size) - yPupil_all(i);
        pupil_speed(i) = sqrt(dx^2 + dy^2) / window_size;
    end
    pupil_speed((N - window_size + 1) : N) = pupil_speed(N - window_size);
end

% Output the full pupil speed array (no median split or logical mask)
pupil_speed_all = pupil_speed;

% 2) Median split using only data points where logical_indices==1
speed_in_mask = pupil_speed(logical_indices);
median_speed  = median(speed_in_mask, 'omitnan');

% 3) Define slow vs fast
slow_pupil = false(N, 1);
fast_pupil = false(N, 1);

slow_condition = (pupil_speed <= median_speed) & logical_indices;
fast_condition = (pupil_speed >  median_speed) & logical_indices;

slow_pupil(slow_condition) = true;
fast_pupil(fast_condition) = true;

end
