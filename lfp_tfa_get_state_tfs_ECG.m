function state_tfs_ecg = lfp_tfa_get_state_tfs_ECG(session_ecg, cond_trials, state, lfp_tfa_cfg)

state_id = state{1};
state_name = state{2};
state_ref_tstart = state{3};
state_ref_tend = state{4};

state_tfs_ecg.powspctrm = {}; % power spectrogram
state_tfs_ecg.time = {}; % timebins fo spectrogram
state_tfs_ecg.freq = {}; % freq bins
state_tfs_ecg.state_id = state_id;
state_tfs_ecg.state_name = state_name;

% loop through trials
for t = find(cond_trials)
    % get the state information for this trial
    states          = session_ecg.trials(t).states;
    state_onset_t   = states([states(:).id] == ...
        state_id).onset_t;
    % rearrange according to tfs timestamps
    state_onset_t   = round(session_ecg.trials(t).tfs.time(...
        abs(session_ecg.trials(t).tfs.time - state_onset_t) == ...
        min(abs(session_ecg.trials(t).tfs.time - state_onset_t))), 4);
    state_start_t   = state_onset_t + state_ref_tstart;
    state_end_t     = state_onset_t + state_ref_tend;
    % sampling frequency
    fs = session_ecg.trials(t).fsample;
    ts = session_ecg.trials(t).tfs.time(2) - session_ecg.trials(t).tfs.time(1);
    
    state_ecg_time = state_ref_tstart-ts:ts:state_ref_tend+ts;
    state_ecg_time = round(state_ecg_time - ...
        state_ecg_time(abs(state_ecg_time) == ...
        min(abs(state_ecg_time))), 4);
    state_ecg_time = state_ecg_time(...
        state_ecg_time >= state_ref_tstart & ...
        state_ecg_time <= state_ref_tend);    
    % crop the tfs for this state
    state_ecg_powspctrm = nan(...
        size(session_ecg.trials(t).tfs.powspctrm, 1), ...
        size(session_ecg.trials(t).tfs.powspctrm, 2), ...
        length(state_ecg_time));
    trial_time_idx = session_ecg.trials(t).tfs.time >= state_start_t & ...
        session_ecg.trials(t).tfs.time <= state_end_t;
    trial_ecg_time = session_ecg.trials(t).tfs.time(trial_time_idx) - state_onset_t;
    trial_ecg_time = round(trial_ecg_time - ...
        trial_ecg_time(abs(trial_ecg_time) == ...
        min(abs(trial_ecg_time))), 4);
    state_ecg_powspctrm(:, :, state_ecg_time >= trial_ecg_time(1) & ...
        state_ecg_time <= trial_ecg_time(end)) = ...
        session_ecg.trials(t).tfs.powspctrm(:, :, trial_time_idx);
%     state_tfs_ecg.powspctrm = [state_tfs_ecg.powspctrm, ...
%         session_ecg.trials(t).tfs.powspctrm(1, :, ...
%         (session_ecg.trials(t).tfs.time >= state_start_t & ...
%         session_ecg.trials(t).tfs.time <= state_end_t))];
    state_tfs_ecg.powspctrm = [state_tfs_ecg.powspctrm, ...
        state_ecg_powspctrm];
    % time bins
%     state_tfs_ecg_time = session_ecg.trials(t).tfs.time(1, ...
%         (session_ecg.trials(t).tfs.time >= state_start_t & ...
%         session_ecg.trials(t).tfs.time <= state_end_t)) - state_onset_t;
%     % put onset timestamp to zero
%     onset_timestamp = state_tfs_ecg_time(...
%         abs(state_tfs_ecg_time) == min(abs(state_tfs_ecg_time)));
%     state_tfs_ecg_time = state_tfs_ecg_time - onset_timestamp;
%     state_tfs_ecg.time = [state_tfs_ecg.time state_tfs_ecg_time];
    state_tfs_ecg.time = [state_tfs_ecg.time state_ecg_time];
    % freq bins
    state_tfs_ecg.freq = session_ecg.trials(t).tfs.freq; 
    state_tfs_ecg.cfg = session_ecg.trials(t).tfs.cfg;

end

% average power spectrum for each state
%arr_state_pow = zeros(1, nfreqbins, ntimebins);

% crop each lfp to same number of samples
% min_ecg_time = state_tfs_ecg.time{cellfun('length', state_tfs_ecg.time) == ...
%     min(cellfun('length', state_tfs_ecg.time))};
% onset_sample = find(min_ecg_time == 0);
% rel_start_sample = 1 - onset_sample;
% rel_end_sample = length(min_ecg_time) - onset_sample;
% for k = 1:length(state_tfs_ecg.time)
%     onset_sample = find(state_tfs_ecg.time{k} == 0);
%     state_tfs_ecg.powspctrm{k} = state_tfs_ecg.powspctrm{k}(:,:,...
%         onset_sample + rel_start_sample:onset_sample + rel_end_sample);
%     state_tfs_ecg.time{k} = state_tfs_ecg.time{k}(...
%         onset_sample + rel_start_sample:onset_sample + rel_end_sample);
% end    
% state_tfs_ecg.time= min_ecg_time; %state_evoked_ecg.ecg_time(1:nsamples);

% % crop each lfp to same number of samples
% nsamples = min(cellfun('size', state_tfs_ecg.powspctrm, 3));
% for k = 1:length(state_tfs_ecg.powspctrm)
%     state_tfs_ecg.powspctrm{k} = state_tfs_ecg.powspctrm{k}(:,:,1:nsamples);
% end
% state_tfs_ecg.time = state_tfs_ecg.time(1:nsamples);

if ~isempty(state_tfs_ecg.powspctrm)

    % find the average TFS for each state
    arr_state_pow = cat(1, state_tfs_ecg.powspctrm{:});
    % remove nans
    state_tfs_ecg.time = state_ecg_time;%(~any(isinf(nansum(arr_state_pow, 1)), 2));
    %arr_state_pow(:, :, any(isinf(nansum(arr_state_pow, 1)), 2)) = [];
    state_tfs_ecg.powspctrm = arr_state_pow;
    state_tfs_ecg.dimord = 'nbeats_freq_time';
    
    state_tfs_ecg.powspctrm_rawmean = nanmean(arr_state_pow, 1);
    
    % baseline normalization
    cfg_baseline.method = lfp_tfa_cfg.baseline_method;
    if ~strcmp(cfg_baseline.method, 'none')
        baseline_cnd_idx = ...
            [session_ecg.baseline.perturbation] == lfp_tfa_cfg.baseline_perturbation & ...
            [session_ecg.baseline.choice] == lfp_tfa_cfg.baseline_use_choice_trial & ...
            [session_ecg.baseline.type] == lfp_tfa_cfg.baseline_use_type & ...
            [session_ecg.baseline.effector] == lfp_tfa_cfg.baseline_use_effector;
        cfg_baseline.mean = session_ecg.baseline(baseline_cnd_idx).pow_mean;
        cfg_baseline.std = session_ecg.baseline(baseline_cnd_idx).pow_std;
        state_tfs_ecg.powspctrm_normmean = lfp_tfa_baseline_normalization(...
            state_tfs_ecg.powspctrm_rawmean, cfg_baseline); 
    else
        state_tfs_ecg.powspctrm_normmean = nanmean(arr_state_pow, 1);
    end
    
end