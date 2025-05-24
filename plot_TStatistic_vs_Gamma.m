function plot_TStatistic_vs_Gamma(t_map, LF_frequencyCenters, HF_frequencyCenters, theta_range)

    % Find the indices corresponding to the theta range
    theta_idx = (LF_frequencyCenters >= theta_range(1)) & (LF_frequencyCenters <= theta_range(2));

    % Average t-values over the selected theta range for each gamma frequency
    mean_t_values = mean(t_map(theta_idx, :), 1);  % Mean across theta indices

    % Apply Gaussian smoothing
    smooth_t_values = smoothdata(mean_t_values, 'gaussian', 3); % Window size = 3

    figure;
    hold on;
    
    % Plot smoothed gamma frequencies first
    plot(HF_frequencyCenters, smooth_t_values, 'k-', 'LineWidth', 2);

    xlabel('Gamma Frequency (Hz)');
    ylabel('Mean T-Statistic');
    %title(sprintf('Visually-Guided Navigation: MI T-Stat vs. Gamma Frequency (Avg %.1f-%.1f Hz Theta)', theta_range(1), theta_range(2)));
    set(gca, 'FontSize', 12);
    grid on;
    
    % Add horizontal reference line
    yline(0, '-k', 'LineWidth', 0.1);

    % Add vertical dashed pink lines at 35 and 45 Hz
    xline(35, '--', 'Color', [1 0.4 0.6], 'LineWidth', 3.5); 
    xline(45, '--', 'Color', [1 0.4 0.6], 'LineWidth', 3.5);

    ax = gca;
    ax.FontSize   = 14; % bigger tick‑label numbers
    ax.TickLength = [0.02 0.02]; % longer ticks

    legend({'Gaussian-smoothened T-Stat'}, 'Location', 'SouthEast');
    
    hold off;
end
