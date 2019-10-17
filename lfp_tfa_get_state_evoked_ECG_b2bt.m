function state_evoked = lfp_tfa_get_state_evoked_ECG_b2bt(trials_lfp, state, lfp_tfa_cfg)

state_id = state{1};
state_name = state{2};
state_reftstart = state{3};
state_reftend = state{4};

state_evoked.ecg_time       = {}; % timebins fo spectrogram
state_evoked.ecg_b2bt       = {}; % evoked LFP response
state_evoked.mean_ecg_b2bt  = [];
state_evoked.state_id       = state_id;
state_evoked.state_name     = state_name;
state_evoked.valid_trials   = [];

ft_data_R2Rt = struct();
ft_data_R2Rt.trial = {};
ft_data_R2Rt.time = {};

% find mean Rpeak interval for all the trials for normalization
R2Rt = [trials_lfp.ECG_b2btime];
Rpeaks = [trials_lfp.ECG_spikes];
mean_ECG_b2btime = nanmean(R2Rt(Rpeaks));    

for t = 1:length(trials_lfp)
        
    if isempty(trials_lfp(t).ECG_b2btime) || any(isnan(trials_lfp(t).ECG_b2btime))
        continue;
    end

    states          = trials_lfp(t).states;
    if ismember(state_id, [states(:).id])
        state_onset_t   = states([states(:).id] == ...
            state_id).onset_t;
        state_onset_sample = states([states(:).id] == ...
            state_id).onset_s;
    else
        continue;
    end
    if isnan(state_onset_t)
        continue;
    end
    
    state_evoked.valid_trials = [state_evoked.valid_trials, t];
    state_start_t   = states([states(:).id] == ...
        state_id).onset_t + state_reftstart;
    state_end_t     = states([states(:).id] == ...
        state_id).onset_t + state_reftend;

    ts = trials_lfp(t).tsample;
    %state_onset_sample = round(state_onset_t/ts) + 1;
    state_start_sample = state_onset_sample + round(state_reftstart/ts);
    state_end_sample = state_onset_sample + round(state_reftend/ts);
    nsamples = state_end_sample - state_start_sample + 1;

    trial_ecg_time = linspace(state_reftstart, state_reftend, nsamples);
    trial_ecg_time = trial_ecg_time - ...
        trial_ecg_time(abs(trial_ecg_time) == ...
        min(abs(trial_ecg_time)));
    
    trial_start_time = states([states(:).id] == ...
        lfp_tfa_cfg.trialinfo.start_state).onset_t + lfp_tfa_cfg.trialinfo.ref_tstart;
    trial_end_time = states([states(:).id] == ...
        lfp_tfa_cfg.trialinfo.end_state).onset_t + lfp_tfa_cfg.trialinfo.ref_tend;
    
    % raw ECG
    trial_evoked_ecg_b2bt = nan(1, nsamples);

    state_start_sample = max(1, state_start_sample);
    state_end_sample = min(length(trials_lfp(t).ECG_b2btime), ...
        state_end_sample);
    ref_onset_sample = find(trial_ecg_time == 0);
    trial_evoked_ecg_b2bt(ref_onset_sample - (state_onset_sample - state_start_sample):...
        ref_onset_sample + (state_end_sample - state_onset_sample)) = ...
        trials_lfp(t).ECG_b2btime(state_start_sample:state_end_sample);
    trial_mean_ecg_b2bt = mean(trials_lfp(t).ECG_b2btime(...
        trials_lfp(t).time >= trial_start_time & trials_lfp(t).time <= trial_end_time));
    
    state_evoked.ecg_time = [state_evoked.ecg_time, trial_ecg_time];
    state_evoked.ecg_b2bt = [state_evoked.ecg_b2bt, ...
        trial_evoked_ecg_b2bt];
    state_evoked.mean_ecg_b2bt = [state_evoked.mean_ecg_b2bt, ...
        trial_mean_ecg_b2bt];
    
    state_onset = false(1, size(trials_lfp(t).ECG_b2btime, 2));
    state_onset(state_onset_sample) = true;
       
    % evoked LFP for this state
    ft_data_R2Rt.trial = [ft_data_R2Rt.trial, ...
        [trials_lfp(t).ECG_b2btime; state_onset]];
    ft_data_R2Rt.label = {'ECG_R2Rt', 'state_onset'};
    %trial_evoked_ecg_b2bt];
    ft_data_R2Rt.time = [ft_data_R2Rt.time, trials_lfp(t).time];
    

end

if ~isempty(ft_data_R2Rt.trial)
    
    % evoked LFP average
%     cfg                 = [];
%     cfg.keeptrials      = 'yes';
%     cfg.timwin          = [state_reftstart state_reftend];
%     cfg.channel         = ft_data_R2Rt.label(1); % ECG R2Rt
%     cfg.spikechannel    = ft_data_R2Rt.label(2); % state onset
%     
%     state_trig_R2Rt     = ft_spiketriggeredaverage(cfg, ft_data_R2Rt);
%     
%     state_evoked.ecg_b2bt       = state_trig_R2Rt.trial;
%     state_evoked.ecg_time       = state_trig_R2Rt.time;
%     state_evoked.mean           = state_trig_R2Rt.avg;
%     state_evoked.std            = permute(nanstd(state_trig_R2Rt.trial, 0, 1), [2 3 1]);
%     
    %state_evoked.mean_ecg_b2bt = mean_ECG_b2btime;
    state_evoked.dimord = 'ntrials_time';
    
    % remove nans
    arr_state_ecg_b2bt = cat(1, state_evoked.ecg_b2bt{:});
    arr_state_ecg_b2bt(isnan(sum(arr_state_ecg_b2bt, 2)), :) = nan;
    state_evoked.ecg_time = trial_ecg_time;%(~any(isinf(arr_state_ecg_b2bt), 1));
    state_evoked.ecg_b2bt = arr_state_ecg_b2bt;
    

    
    % evoked LFP average
    if lfp_tfa_cfg.normalize_R2Rt
        norm_arr_state_ecg_b2bt = arr_state_ecg_b2bt ./ ...
            repmat(state_evoked.mean_ecg_b2bt', ...
            [1 size(arr_state_ecg_b2bt, 2)]);
        state_evoked.mean = nanmean(norm_arr_state_ecg_b2bt, 1);
        state_evoked.std = nanstd(norm_arr_state_ecg_b2bt, 0, 1);
    else
        state_evoked.mean = nanmean(arr_state_ecg_b2bt, 1);
        state_evoked.std = nanstd(arr_state_ecg_b2bt, 0, 1);
    end
    state_evoked.ntrials = sum(~isnan(sum(arr_state_ecg_b2bt, 2)));
    
end