function [IntoBoundary_logical, OutofBoundary_logical, StandingBoundary_logical] = ...
    find_InOut_boundary_indices(logical_indices_memory, logical_indices_navigation, ...
    logical_indices_waiting, boundary_logical)

    % Identify segments in each condition and within boundary threshold
    IntoBoundary_logical = logical_indices_navigation & boundary_logical;
    OutofBoundary_logical = logical_indices_memory & boundary_logical;
    StandingBoundary_logical = logical_indices_waiting & boundary_logical;

%     % Compute data counts
%     boundaryData = sum(boundary_logical);
%     StandingBoundary = sum(StandingBoundary_logical);
%     MemoryBoundary = sum(OutofBoundary_logical);
%     MovementBoundary = sum(IntoBoundary_logical);
%     
%     unclassifiedBoundary = boundaryData - (StandingBoundary + MemoryBoundary + MovementBoundary);
% 
%     % Ensure balance within a tolerance threshold
%     boundary_counts = [StandingBoundary, MemoryBoundary, MovementBoundary];
%     median_count = median(boundary_counts);
%     tolerance = 0.05 * median_count;  % Allow 5% deviation
% 
%     % Check if any group is significantly different
%     imbalanced_groups = boundary_counts < (median_count - tolerance) | boundary_counts > (median_count + tolerance);
% 
%     % Display warning if there are unclassified or imbalanced boundary segments
%     if unclassifiedBoundary ~= 0
%         warning('Unclassified boundary segments detected: %d', unclassifiedBoundary);
%     end
% 
%     if any(imbalanced_groups)
%         warning('Boundary condition sample sizes are significantly different.');
%         disp(table(["StandingBoundary"; "MemoryBoundary"; "MovementBoundary"], boundary_counts'));
%     else
%         disp('Boundary condition sample sizes are within a reasonable tolerance.');
%     end

end
