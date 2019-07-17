function combined_tfs = lfp_tfa_get_combined_tfs( site_lfp, cond_trials, state, lfp_tfa_cfg )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

start_state = state{2}(1);
end_state = state{2}(2);
width = state{3};
nwindows = state{4};
spacing = state{5};

combined_tfs.powspctrm = {}; % power spectrogram
combined_tfs.time = {}; % timebins fo spectrogram
combined_tfs.freq = {}; % freq bins

for t = 1:find(cond_trials)

    states                = site_lfp.trials(t).states;
    state_start_onset_t   = states([states(:).id] == ...
        start_state).onset_t;
    end_state_onset_t     = states([states(:).id] == ...
        end_state).onset_t;
    
    % get the LFP samples and timestamps between start and end states
    lfp_tfs_pow = site_lfp.trials(t).tfs.powspctrm(:, :, ...
        (site_lfp.trials(t).tfs.time >= state_start_onset_t & ...
        site_lfp.trials(t).tfs.time <= end_state_onset_t));
    lfp_tfs_time = site_lfp.trials(t).tfs.time(...
        (site_lfp.trials(t).tfs.time >= state_start_onset_t & ...
        site_lfp.trials(t).tfs.time <= end_state_onset_t));
    lfp_tfs_freq = site_lfp.trials(t).tfs.freq;
    % lfp sample indices
    lfp_timebin_idx = 1:length(lfp_tfs_time);
    % width of each time bin
    lfp_ts = 1/site_lfp.trials(t).fsample;
    tbin_width = lfp_tfa_cfg.tfr.timestep * lfp_ts;
    
    % number fo time bins in each window
    ntbins_window = round(width/tbin_width);
    
    % check if LFP data is sufficiently long to extract the required number
    % of windows of specified length
    if nwindows*(ntbins_window) > length(lfp_tfs_time)
        fprintf('LFP is too short to extract the required windows. \n');
        return
    end
    
    combined_tfs.cfg = site_lfp.trials(t).tfs.cfg;
    
    % now get the windows to combine
    if strcmp(spacing, 'uniform') % uniformly spaced windows
        window_spacing = round(length(lfp_tfs_time)/nwindows);
        window_mid_idx = round(window_spacing/2):window_spacing:...
            length(lfp_tfs_time) - round(window_spacing/2);
        
        % loop through each window
        for w = 1:length(window_mid_idx)
            % evoked LFP for this state
            combined_tfs.powspctrm = [combined_tfs.powspctrm, ...
                lfp_tfs_pow(:, :, max(1, window_mid_idx(w) - round(ntbins_window/2)):...
                min(length(lfp_tfs_time), window_mid_idx(w) + round(ntbins_window/2)))];
            % timestamps
            combined_tfs.time = ...
                lfp_tfs_time(max(1, window_mid_idx(w) - round(ntbins_window/2)):...
                min(length(lfp_tfs_time), window_mid_idx(w) + round(ntbins_window/2)));
            % set mid-timestamp to zero
            combined_tfs.time = combined_tfs.time - ...
                combined_tfs.time(round(length(combined_tfs.time)/2));
            
            combined_tfs.freq = lfp_tfs_freq;
        end
        
    end
    
    if strcmp(spacing, 'random') % randomly spaced windows
        window_spacing = round(length(lfp_tfs_time)/nwindows);
        window_start_idx = 1:window_spacing:...
            length(lfp_tfs_time) - ntbins_window;
        window_start_idx = window_start_idx(1:nwindows);
                
        % loop through each window
        for w = 1:length(window_start_idx)
            window_sample_idx = lfp_timebin_idx(window_start_idx(w) + round(ntbins_window/2):min(window_start_idx(w) + ...
                window_spacing - round(ntbins_window/2), length(lfp_tfs_time) - round(ntbins_window/2)));
            % pick a random sample as window middle
            window_mid_idx = round((window_sample_idx(end) - window_sample_idx(1)) ...
                * rand) + window_sample_idx(1);%randsample(window_sample_idx, 1, true);
            %if p_window > rand % a window occurs
            % evoked LFP for this window
            combined_tfs.powspctrm = [combined_tfs.powspctrm, ...
                lfp_tfs_pow(:, :, max(1, window_mid_idx - round(ntbins_window/2)):...
                min(length(lfp_tfs_time), window_mid_idx + floor(ntbins_window/2)))];
            % timestamps
            combined_tfs.time = lfp_tfs_time(max(1, window_mid_idx - round(ntbins_window/2)):...
                min(length(lfp_tfs_time), window_mid_idx + round(ntbins_window/2)));
            % set mid-timestamp to zero
            combined_tfs.time = combined_tfs.time - ...
                combined_tfs.time(round(length(combined_tfs.time)/2));
                       
            combined_tfs.freq = lfp_tfs_freq;
            %end
            %if wnd, continue, end;
        end
    end
    combined_tfs.cfg = site_lfp.trials(t).tfs.cfg;            
end


% find number of time bins in power
% spectrogram
ntimebins = min(cellfun('size', combined_tfs.powspctrm, 3));
nfreqbins = numel(combined_tfs.freq);
% crop each tfs to the ntimebins
for k = 1:length(combined_tfs.powspctrm)
    combined_tfs.powspctrm{k} = combined_tfs.powspctrm{k}(1,:,1:ntimebins);
    
end
combined_tfs.time = combined_tfs.time(1:ntimebins);

% average power spectrum for each state
arr_state_pow = zeros(1, nfreqbins, ntimebins);

if ~isempty(combined_tfs.powspctrm)

    % find the average TFS for each state
    arr_state_pow = cat(1, combined_tfs.powspctrm{:});
    combined_tfs.powspctrm_rawmean = nanmean(arr_state_pow, 1);

    % baseline normalization
    cfg_baseline.method = lfp_tfa_cfg.baseline_method;
    baseline_cnd_idx = [site_lfp.baseline.perturbation] == ...
        lfp_tfa_cfg.baseline_perturbation & [site_lfp.baseline.choice] == ...
        lfp_tfa_cfg.baseline_use_choice_trial;
    cfg_baseline.mean = site_lfp.baseline(baseline_cnd_idx).pow_mean;
    cfg_baseline.std = site_lfp.baseline(baseline_cnd_idx).pow_std;
    combined_tfs.powspctrm_normmean = lfp_tfa_baseline_normalization(...
        combined_tfs.powspctrm_rawmean, cfg_baseline); 
    
end

end

