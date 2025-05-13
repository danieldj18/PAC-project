function [badData, badUnfiltered, badFiltered] = detect_IED(Fs, eegsignal, filterRange, std_unfilt, std_filt)
    %{
    detect_IED: Identify IED samples in EEG data using thresholds on unfiltered and filtered signals.

    Inputs:
    - ts: Time vector (if Fs not provided)
    - eegsignal: The EEG signal (vector of time-series data).
    - Fs: Sampling frequency (Hz).
    - filterRange: Frequency range for filtering (default [25 80] Hz).
    - std_unfilt: Standard deviation multiplier for unfiltered signal (default 3).
    - std_filt: Standard deviation multiplier for filtered signal (default 3).

    Outputs:
    - badData: Logical array marking IED samples (1) and non-IED samples (0).
    - badUnfiltered: Logical array marking IED samples based on unfiltered signal.
    - badFiltered: Logical array marking IED samples based on filtered signal.
    %}

    % Set default parameters if not provided
    if nargin < 6; std_filt = 6; end % Default standard deviation multiplier for filtered signal
    if nargin < 5; std_unfilt = 6; end % Default standard deviation multiplier for unfiltered signal
    if nargin < 4; filterRange = [15 80]; end % Default filter range

    % Step 1: Rectify the original signal and compute the envelope
    unfiltered_baseline = nanmedian(eegsignal);
    [unfiltered_envelope, ~] = envelope(eegsignal);
    unfiltered_threshold = unfiltered_baseline + std_unfilt * nanstd(unfiltered_envelope);
    badUnfiltered = unfiltered_envelope > unfiltered_threshold;

    % Step 2: Filter the signal and compute the envelope
    [b, a] = butter(4, filterRange * 2 / Fs, 'bandpass');
    filtered_signal = filtfilt(b, a, eegsignal);
    filtered_baseline = nanmedian(filtered_signal);
    [filtered_envelope, ~] = envelope(filtered_signal);
    filtered_threshold = filtered_baseline + std_filt * nanstd(filtered_envelope);
    badFiltered = filtered_envelope > filtered_threshold;

    % Step 3: Combine the bad samples from both criteria
    badData = badUnfiltered | badFiltered;
end
