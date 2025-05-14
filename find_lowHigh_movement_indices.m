function [logical_indices_lowMovement, logical_indices_highMovement, ...
          logical_indices_lowTurn, logical_indices_highTurn, ...
          logical_indices_veryLowMovement, logical_indices_restMovement, ...
          turning_speed_all, angular_acceleration_all] = ...
          find_lowHigh_movement_indices(memory_mask, speed_all, movementDir_all, turningWindowSize, smoothingWindowSize)

% INPUTS:
%   memory_mask          - [1 x N] logical/numeric mask, 1 => in memory period.
%   speed_all            - [1 x N] array, movement speed per sample in m/s (250 Hz).
%   movementDir_all      - [1 x N] array, movement direction in degrees (0–360) per sample.
%   turningWindowSize    - (optional) #samples for turning angle calculation (default=25).
%   smoothingWindowSize  - (optional) #samples for smoothing direction (default=0, no smoothing).
%
% OUTPUTS:
%   logical_indices_lowMovement       - [1 x N] logical, speed < median (within memory_mask).
%   logical_indices_highMovement      - [1 x N] logical, speed >= median (within memory_mask).
%   logical_indices_veryLowMovement   - [1 x N] logical, speed < 0.2 m/s.
%   logical_indices_restMovement      - [1 x N] logical, speed >= 0.2 m/s (rest = rest of the, not "resting").
%   logical_indices_lowTurn           - [1 x N] logical, turn angle <= median (within memory_mask).
%   logical_indices_highTurn          - [1 x N] logical, turn angle > median (within memory_mask).
%   turning_speed_all                 - [1 x N] continuous, turning speed in degrees/sample.
%   angular_acceleration_all          - [1 x N] continuous, angular acceleration in degrees/sample².
%
% DESCRIPTION OF KEY CALCULATIONS:
%   - Movement splits: Median-based on speed within memory_mask.
%   - Very low/rest movement: < 0.2 m/s.
%   - Turn angle: Backward-looking window difference (i - turningWindowSize) to i (degrees).
%   - Turning speed: abs(diff(direction)) per sample, length-consistent.
%   - Angular acceleration: diff(turning_speed_all) with first valid acceleration replication.
%   - Balance san checks: ~5% imbalance warnings for movement and turning splits.

    if nargin < 4
        turningWindowSize = 25; 
    end
    if nargin < 5
        smoothingWindowSize = 0;
    end

    memory_mask = logical(memory_mask);
    N = length(speed_all);

    if length(movementDir_all) ~= N
        error('speed_all and movementDir_all must match in length!');
    end
    if length(memory_mask) ~= N
        error('memory_mask must match speed_all in length!');
    end

    % 1) Median split on movement speed
    memory_speed = speed_all(memory_mask);
    medianSpeed = median(memory_speed, 'omitnan');

    logical_indices_lowMovement = false(1, N);
    logical_indices_highMovement = false(1, N);
    logical_indices_lowMovement(memory_mask & speed_all < medianSpeed) = true;
    logical_indices_highMovement(memory_mask & speed_all >= medianSpeed) = true;

    % 2) Additional Output - Speed < 0.2 and Rest
    logical_indices_veryLowMovement = false(1, N);
    logical_indices_restMovement = false(1, N);
    logical_indices_veryLowMovement(memory_mask & speed_all < 0.2) = true;
    logical_indices_restMovement(memory_mask & speed_all >= 0.2) = true;

    % 3) Smooth heading data (optional)
    if smoothingWindowSize > 1
        movementDir_smooth = movmean(movementDir_all, smoothingWindowSize);
    else
        movementDir_smooth = movementDir_all;
    end

    % 4) Compute Turning Angle
    turn_angle = nan(1, N);
    for i = (turningWindowSize + 1):N
        raw_diff = movementDir_smooth(i) - movementDir_smooth(i - turningWindowSize);
        turn_angle(i) = abs(raw_diff);
    end
    
    % Handle the first 'turningWindowSize' points by replicating the first valid angle
    if turningWindowSize >= 1
        turn_angle(1:turningWindowSize) = turn_angle(turningWindowSize + 1);
    else
        turn_angle(:) = 0;
    end

    % 5) Compute Turning Speed (per sample)
    turning_speed_all = [0, abs(diff(movementDir_smooth))];
    if length(turning_speed_all) < N
        turning_speed_all(end+1:N) = turning_speed_all(end);
    end

    % 6) Compute Angular Acceleration
    angular_acceleration_all = [diff(turning_speed_all(1:2)), diff(turning_speed_all)];
    if length(angular_acceleration_all) < N
        angular_acceleration_all(end+1:N) = angular_acceleration_all(end);
    end

    % 7) Median split for turning angle within memory_mask
    memory_turn = turn_angle(memory_mask);
    medianTurn = median(memory_turn, 'omitnan');

    logical_indices_lowTurn = false(1, N);
    logical_indices_highTurn = false(1, N);
    logical_indices_lowTurn(memory_mask & turn_angle <= medianTurn) = true;
    logical_indices_highTurn(memory_mask & turn_angle > medianTurn) = true;

    % 8) Check group size imbalance (~5% rule)
    movement_counts = [sum(logical_indices_lowMovement), sum(logical_indices_highMovement)];
    turn_counts = [sum(logical_indices_lowTurn), sum(logical_indices_highTurn)];
    veryLowMovement_count = sum(logical_indices_veryLowMovement);
    restMovement_count = sum(logical_indices_restMovement);

    % Movement balance check
    if abs(movement_counts(1) - movement_counts(2)) / mean(movement_counts) > 0.05
        warning('Low vs. High Movement classification is imbalanced.');
        disp(table(["Low Movement"; "High Movement"], movement_counts'));
    end

    % Turning balance check
    if abs(turn_counts(1) - turn_counts(2)) / mean(turn_counts) > 0.05
        warning('Low vs. High Turn classification is imbalanced.');
        disp(table(["Low Turn"; "High Turn"], turn_counts'));
    end

    % Very Low vs. Rest Movement check
    if abs(veryLowMovement_count - restMovement_count) / mean([veryLowMovement_count, restMovement_count]) > 0.05
        warning('Very Low (<0.2) vs. Rest Movement classification is imbalanced.');
        disp(table(["Very Low Movement (<0.2)"; "Rest Movement"], [veryLowMovement_count; restMovement_count]));
    end
end
