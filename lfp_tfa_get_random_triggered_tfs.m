function ecg_triggered_tfs = lfp_tfa_get_random_triggered_tfs( site_lfp, cond_trials, state, lfp_tfa_cfg )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

state_name = state{2};
width = state{4} - state{3};

ecg_triggered_tfs.powspctrm = {}; % power spectrogram
ecg_triggered_tfs.time = {}; % timebins fo spectrogram
ecg_triggered_tfs.freq = {}; % freq bins
ecg_triggered_tfs.state = [];
ecg_triggered_tfs.state_name = state_name;

for t = find(cond_trials)

    trialperiod           = site_lfp.trials(t).trialperiod;
    
    % get the LFP samples and timestamps between start and end states
    lfp_tfs_pow = site_lfp.trials(t).tfs.powspctrm(:, :, ...
        (site_lfp.trials(t).tfs.time >= trialperiod(1) & ...
        site_lfp.trials(t).tfs.time <= trialperiod(2)));
    lfp_tfs_time = site_lfp.trials(t).tfs.time(...
        (site_lfp.trials(t).tfs.time >= trialperiod(1) & ...
        site_lfp.trials(t).tfs.time <= trialperiod(2)));
    lfp_tfs_freq = site_lfp.trials(t).tfs.freq;
    
    % ecg peak times
    lfp_time = site_lfp.trials(t).time(...
        (site_lfp.trials(t).time >= trialperiod(1) & ...
        site_lfp.trials(t).time <= trialperiod(2))); 
    ecg_peaks = site_lfp.trials(t).ECG_spikes(...
        (site_lfp.trials(t).time >= trialperiod(1) & ...
        site_lfp.trials(t).time <= trialperiod(2))); 
    ecg_peak_times = lfp_time(ecg_peaks);
    
    mean_ecg_p2pt = mean(diff(ecg_peak_times));
    
    % lfp sample indices
    lfp_timebin_idx = 1:length(lfp_tfs_time);
    % width of each time bin
    lfp_ts = 1/site_lfp.trials(t).fsample;
    tbin_width = lfp_tfa_cfg.tfr.timestep * lfp_ts;
    
    % number fo time bins in the window
    ntbins_window = round(width/tbin_width);
    
    ecg_triggered_tfs.cfg = site_lfp.trials(t).tfs.cfg;
    
    % now get the midpoint of windows around ECG peak
    window_mid_idx = [];
    for p = ecg_peak_times
        if min(abs(lfp_tfs_time - p)) > tbin_width
            continue;
        end
        window_mid_idx = [window_mid_idx, find(abs(lfp_tfs_time - p) == ...
            min(abs(lfp_tfs_time - p)))];
    end
    
    % randomly shift peaks
    window_mid_idx = round(window_mid_idx + (mean_ecg_p2pt/tbin_width)*...
        (rand(size(window_mid_idx))*2 - 1));

    % loop through each window
    for w = 1:length(window_mid_idx)
        if window_mid_idx(w) - round(ntbins_window/2) < 1 || ...
                window_mid_idx(w) + round(ntbins_window/2) > length(lfp_tfs_time)
            continue;
        end
        % power spectrum for this window
        ecg_triggered_tfs.powspctrm = [ecg_triggered_tfs.powspctrm, ...
            lfp_tfs_pow(:, :, window_mid_idx(w) - round(ntbins_window/2):...
             window_mid_idx(w) + round(ntbins_window/2))];
        % timestamps
        ecg_triggered_tfs.time = ...
            lfp_tfs_time(window_mid_idx(w) - round(ntbins_window/2):...
             window_mid_idx(w) + round(ntbins_window/2));
        % set mid-timestamp to zero
        ecg_triggered_tfs.time = ecg_triggered_tfs.time - ...
            ecg_triggered_tfs.time(round(length(ecg_triggered_tfs.time)/2));

        ecg_triggered_tfs.freq = lfp_tfs_freq;
    end
       
end


% find number of time bins in power
% spectrogram
ntimebins = min(cellfun('size', ecg_triggered_tfs.powspctrm, 3));
nfreqbins = numel(ecg_triggered_tfs.freq);
% crop each tfs to the ntimebins
for k = 1:length(ecg_triggered_tfs.powspctrm)
    ecg_triggered_tfs.powspctrm{k} = ecg_triggered_tfs.powspctrm{k}(1,:,1:ntimebins);
    
end
ecg_triggered_tfs.time = ecg_triggered_tfs.time(1:ntimebins);

% average power spectrum for each state
arr_state_pow = zeros(1, nfreqbins, ntimebins);

if ~isempty(ecg_triggered_tfs.powspctrm)

    % find the average TFS for each state
    arr_state_pow = cat(1, ecg_triggered_tfs.powspctrm{:});
    ecg_triggered_tfs.powspctrm_rawmean = nanmean(arr_state_pow, 1);

    % baseline normalization
    cfg_baseline.method = lfp_tfa_cfg.baseline_method;
    if ~strcmp(cfg_baseline.method, 'none')
        if isnan(lfp_tfa_cfg.baseline_perturbation)
            baseline_cnd_idx = isnan([site_lfp.baseline.perturbation]);
        end
        if isnan(lfp_tfa_cfg.baseline_use_choice_trial)
            baseline_cnd_idx = baseline_cnd_idx & isnan([site_lfp.baseline.choice]);
        end
        if isnan(lfp_tfa_cfg.baseline_use_type)
            baseline_cnd_idx = baseline_cnd_idx & isnan([site_lfp.baseline.type]);
        end
        if isnan(lfp_tfa_cfg.baseline_use_effector)
            baseline_cnd_idx = baseline_cnd_idx & isnan([site_lfp.baseline.effector]);
        end
        if ~isnan(lfp_tfa_cfg.baseline_perturbation) && ~isnan(lfp_tfa_cfg.baseline_use_choice_trial)
            baseline_cnd_idx = ...
                [site_lfp.baseline.perturbation] == lfp_tfa_cfg.baseline_perturbation & ...
                [site_lfp.baseline.choice] == lfp_tfa_cfg.baseline_use_choice_trial & ...
                [site_lfp.baseline.type] == lfp_tfa_cfg.baseline_use_type & ...
                [site_lfp.baseline.effector] == lfp_tfa_cfg.baseline_use_effector;
        end
        cfg_baseline.mean = site_lfp.baseline(baseline_cnd_idx).pow_mean;
        cfg_baseline.std = site_lfp.baseline(baseline_cnd_idx).pow_std;
        ecg_triggered_tfs.powspctrm_normmean = lfp_tfa_baseline_normalization(...
            ecg_triggered_tfs.powspctrm_rawmean, cfg_baseline); 
    else
        ecg_triggered_tfs.powspctrm_normmean = ecg_triggered_tfs.powspctrm_rawmean;
    end
    
    
end

end

