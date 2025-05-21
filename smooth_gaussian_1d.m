function smoothed_signal = smooth_gaussian_1d(signal, Fs, smoothing_win)
    % smooth_gaussian_1d: Apply Gaussian smoothing to a 1D signal using sampling frequency (Fs).
    %
    % Inputs:
    % - signal: The 1D signal (e.g., IED_samples) to be smoothed.
    % - Fs: The sampling frequency (in Hz).
    % - smoothing_win: The standard deviation (width) of the Gaussian window (in seconds).
    %
    % Output:
    % - smoothed_signal: The smoothed signal after applying Gaussian filtering.

    % Calculate the number of samples that correspond to the window size in seconds
    window_samples = round(smoothing_win * Fs); 

    % Create a Gaussian kernel (window)
    kernel_half_size = window_samples;  % The size of the kernel (half)
    gauss_time = linspace(-kernel_half_size, kernel_half_size, 2 * kernel_half_size + 1);  % Time axis for the kernel
    gauss_kernel = exp(-gauss_time.^2 / (2 * (window_samples / 2)^2));  % Gaussian formula

    % Normalize the Gaussian kernel to ensure that it sums to 1
    gauss_kernel = gauss_kernel / sum(gauss_kernel);

    % Apply the Gaussian smoothing using convolution
    smoothed_signal = conv(signal, gauss_kernel, 'same');  % 'same' = output the same size as the input signal

end
