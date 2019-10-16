function state_tfs = lfp_tfa_get_state_tfs(site_lfp, cond_trials, state, lfp_tfa_cfg)

state_id = state{2};
state_name = state{3};
state_ref_tstart = state{4};
state_ref_tend = state{5};

state_tfs.powspctrm = {}; % power spectrogram
state_tfs.time = {}; % timebins fo spectrogram
state_tfs.freq = {}; % freq bins
state_tfs.state_id = state_id;
state_tfs.state_name = state_name;

% loop through trials
for t = find(cond_trials)
    % get the state information for this trial
    states          = site_lfp.trials(t).states;
    state_onset_t   = states([states(:).id] == ...
        state_id).onset_t;
    state_start_t   = states([states(:).id] == ...
        state_id).onset_t + state_ref_tstart;
    state_end_t     = states([states(:).id] == ...
        state_id).onset_t + state_ref_tend;
    % sampling frequency
    fs = site_lfp.trials(t).fsample;

    % crop the tfs for this state
    state_tfs.powspctrm = [state_tfs.powspctrm, ...
        site_lfp.trials(t).tfs.powspctrm(1, :, ...
        (site_lfp.trials(t).tfs.time >= state_start_t & ...
        site_lfp.trials(t).tfs.time <= state_end_t))];
    % time bins
    state_tfs.time = site_lfp.trials(t).tfs.time(1, ...
        (site_lfp.trials(t).tfs.time >= state_start_t & ...
        site_lfp.trials(t).tfs.time <= state_end_t)) - state_onset_t;
    % put onset timestamp to zero
    onset_timestamp = state_tfs.time(...
        abs(state_tfs.time) == min(abs(state_tfs.time)));
    state_tfs.time = state_tfs.time - onset_timestamp;
    % freq bins
    state_tfs.freq = site_lfp.trials(t).tfs.freq; 
    state_tfs.cfg = site_lfp.trials(t).tfs.cfg;

end

% find number of time bins in power
% spectrogram
ntimebins = min(cellfun('size', state_tfs.powspctrm, 3));
nfreqbins = numel(state_tfs.freq);
% crop each tfs to the ntimebins
for k = 1:length(state_tfs.powspctrm)
    state_tfs.powspctrm{k} = state_tfs.powspctrm{k}(1,:,1:ntimebins);
    %
end
state_tfs.time = state_tfs.time(1:ntimebins);

% average power spectrum for each state
arr_state_pow = zeros(1, nfreqbins, ntimebins);

if ~isempty(state_tfs.powspctrm)

    % find the average TFS for each state
    state_tfs.powspctrm = cat(1, state_tfs.powspctrm{:});
    state_tfs.powspctrm_rawmean = nanmean(state_tfs.powspctrm, 1);

    % baseline normalization
    cfg_baseline.method = lfp_tfa_cfg.baseline_method;
    baseline_cnd_idx = [site_lfp.baseline.perturbation] == ...
        lfp_tfa_cfg.baseline_perturbation & [site_lfp.baseline.choice] == ...
        lfp_tfa_cfg.baseline_use_choice_trial;
    cfg_baseline.mean = site_lfp.baseline(baseline_cnd_idx).pow_mean;
    cfg_baseline.std = site_lfp.baseline(baseline_cnd_idx).pow_std;
    state_tfs.powspctrm = lfp_tfa_baseline_normalization(...
        state_tfs.powspctrm, cfg_baseline); 
    state_tfs.baseline = cfg_baseline;
    
end