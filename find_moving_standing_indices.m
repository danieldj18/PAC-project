function [moving_logical, standing_logical] = find_moving_standing_indices(...
    movement_logical, waiting_logical, movement_speed, moving_threshold)

    % Initialize logical arrays for moving and standing phases
    moving_logical = false(size(movement_logical));  
    standing_logical = false(size(waiting_logical));

    % Classify moving indices: movement phase and speed above threshold
    moving_logical = movement_logical & (movement_speed > moving_threshold);

    % Classify standing indices: waiting phase and speed below or equal to threshold
    standing_logical = waiting_logical & (movement_speed <= moving_threshold);
   
    
    checkData = false;
    if checkData
        % **Check for Balanced Data Distribution**
        moving_count = sum(moving_logical);
        standing_count = sum(standing_logical);
        
        % Define a reasonable tolerance (5% deviation from median count)
        category_counts = [moving_count, standing_count];
        median_count = median(category_counts);
        tolerance = 0.05 * median_count;
    
        % Identify imbalanced groups
        imbalanced_groups = category_counts < (median_count - tolerance) | category_counts > (median_count + tolerance);
    
        % Display warnings if sample sizes are too different
        if any(imbalanced_groups)
            warning('Movement vs. Standing classification has imbalanced sample sizes.');
            disp(table(["Moving"; "Standing"], category_counts'));
        else
            disp('Movement and Standing classifications are within a reasonable balance.');
        end
    end
end
