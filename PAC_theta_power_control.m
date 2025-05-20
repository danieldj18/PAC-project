clear all
close all

patient_number = 5;

filename = 'PAC_thetaControl.mat';

% Excluded channels
excluded_channels = struct();
excluded_channels.subj1 = [4];
excluded_channels.subj2 = [];
excluded_channels.subj3 = [2,4]; 
excluded_channels.subj4 = [];
excluded_channels.subj5 = [4];


all_results = struct();


for subj = 1:patient_number

    data_name = sprintf('data_preprocessed_subj%d.mat', subj);
    
    % Run the trial-based PAC analysis
    all_results.(['subj', num2str(subj)]) = ...
        analyze_trialBased_MI(data_name, excluded_channels, subj, true, 250);
    
    disp(['Processed subject ', num2str(subj)]);
end


save(filename, 'all_results');


%% Regression (subj level)

load('PAC_thetaControl.mat', 'all_results');

% Get subject names
subjectNames = fieldnames(all_results);
nSubjects = length(subjectNames);

subject_beta = nan(nSubjects,1);

% Loop over subjects to compute subject-level beta coefficients
for s = 1:nSubjects
    subjName = subjectNames{s};
    subjectResults = all_results.(subjName).navigation;
    channelNames = fieldnames(subjectResults);
    
    % For each subject, collect the beta slopes computed per channel
    channel_betas = [];
    
    for c = 1:length(channelNames)
        chanKey = channelNames{c};
        
        % Retrieve trial-level MI and theta power for memory trials
        mem_MI_trials = subjectResults.(chanKey).memory_MI_trials;            
        mem_theta_trials = subjectResults.(chanKey).memory_thetaPower_trials;  
        
        if isempty(mem_MI_trials) || isempty(mem_theta_trials)
            continue;
        end
        
        % Average each trial's MI
        avg_chan_MI = mem_MI_trials(:);
        
        % Convert to column vectors: each trial is one data point
        Y = avg_chan_MI(:); % Dependent variable: MI
        X = mem_theta_trials(:); % Predictor variable: theta power
        
        % Add an intercept term for the regression model
        X_reg = [ones(length(X),1) X];
        
        % Only if at least 3 trials
        if length(X) < 3
            continue;
        end
        
        % Run linear regression to obtain beta coefficients
        beta = regress(Y, X_reg);  % beta(2) is the slope for theta power
        channel_betas(end+1) = beta(2); 
    end
    
    % Average the beta coefficients across channels for the subject
    subject_beta(s) = mean(channel_betas, 'omitnan');
    fprintf('Subject %s: Average beta = %.4f\n', subjName, subject_beta(s));
end

% Perform permutation test on subject-level beta values
nPerm = 10000;
observedMeanBeta = mean(subject_beta, 'omitnan');
permMeanBeta = nan(nPerm,1);

% Permutation test via sign-flipping subject betas
for i = 1:nPerm
    signFlip = randi([0,1], nSubjects, 1)*2 - 1;  % Randomly +1 or -1
    permutedBeta = subject_beta .* signFlip;
    permMeanBeta(i) = mean(permutedBeta, 'omitnan');
end


p_perm = sum(abs(permMeanBeta) >= abs(observedMeanBeta)) / nPerm;
fprintf('Permutation test: p = %.3f\n', p_perm);

%% Regression (Channel Level)

% Load the trial-based PAC results
load('PAC_thetaControl.mat', 'all_results');

% Collect all channel-level beta coefficients
all_channel_betas = [];

subjectNames = fieldnames(all_results);
for s = 1:numel(subjectNames)
    subjName = subjectNames{s};
    navResults = all_results.(subjName).navigation;
    channelNames = fieldnames(navResults);
    
    for c = 1:numel(channelNames)
        chanKey = channelNames{c};
        mem_MI_trials = navResults.(chanKey).memory_MI_trials;
        mem_theta_trials = navResults.(chanKey).memory_thetaPower_trials;
        
        if isempty(mem_MI_trials) || isempty(mem_theta_trials)
            continue;
        end
        
        % average across rows
        avg_MI = mem_MI_trials;

        
        Y = avg_MI(:);
        X = mem_theta_trials(:);
        if numel(X) < 3
            continue;
        end
        
        X_reg = [ones(numel(X),1), X];
        beta = regress(Y, X_reg); % beta(2) = slope
        all_channel_betas(end+1,1) = beta(2);
    end
end

% Compute observed mean across all channels
observedMean = mean(all_channel_betas, 'omitnan');

% Permutation via sign‐flips at the channel level
nPerm = 10000;
permMeans = nan(nPerm,1);
nChans = numel(all_channel_betas);

for i = 1:nPerm
    flips = (randi([0,1], nChans,1)*2 - 1);
    permutedBetas = all_channel_betas .* flips;
    permMeans(i) = mean(permutedBetas, 'omitnan');
end

% p-val
p_chan = sum(abs(permMeans) >= abs(observedMean)) / nPerm;
fprintf('Channel-level permutation test: p = %.4f\n', p_chan);


%% Example figure: Trial-Level Regression Analysis for Subject 1, Channel 1

subjectResults = all_results.subj1.navigation.channel_1;

mem_MI_trials = subjectResults.memory_MI_trials;            
mem_theta_trials = subjectResults.memory_thetaPower_trials;  

% The MI per trial is stored as a matrix (multiple values per trial),
% average across the MI vector for each trial to yield one MI val per trial.
if size(mem_MI_trials, 1) > 1
    trial_avg_MI = mean(mem_MI_trials, 1, 'omitnan'); % 1 x nTrials
else
    trial_avg_MI = mem_MI_trials; % one val per trial
end

trial_theta_power = mem_theta_trials;  

Y = trial_avg_MI(:); % Dependent variable: MI (z-scored)
X = trial_theta_power(:); % Predictor variable: Theta Power (1-3 Hz)

% Scatter plot for Subject 1, Channel 1
figure;
scatter(X, Y, 100, 'filled');
xlabel('Theta Power (1-3 Hz; µV²)');
ylabel('Modulation Index (z-scored)');
title('Trial-by-Trial Scatter: MI vs. Theta Power (Subj 1, Chan 1)');
lsline;

% Linear regression with an intercept
X_reg = [ones(length(X),1) X];
beta = regress(Y, X_reg);
fprintf('Subject 1, Channel 1: Regression slope (beta) = %.4f\n', beta(2));

% Pearson correlation coef
[r, p] = corr(X, Y, 'Type', 'Pearson');
fprintf('Subject 1, Channel 1: Pearson r = %.3f, p = %.3f\n', r, p);


str = sprintf('$$\\beta = %.4f, \\; r = %.3f, \\; p = %.3f$$', beta(2), r, p);

xPad = 0.05 * (max(X) - min(X));
yPad = 0.05 * (max(Y) - min(Y));
xPos = min(X) + xPad;
yPos = max(Y) - yPad;

h = lsline;  % least-squares regres
set(h, 'LineWidth', 3, 'Color', 'k');

text(xPos, yPos, str, ...
    'Interpreter', 'latex', ...
    'FontSize', 12, ...
    'BackgroundColor', 'w', ...
    'EdgeColor', [0.8 0.8 0.8], ...
    'Margin', 6);
