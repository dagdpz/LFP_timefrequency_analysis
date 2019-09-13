function state_evoked = lfp_tfa_get_state_evoked_ECG_b2bt(trials_lfp, state)

state_id = state{1};
state_name = state{2};
state_reftstart = state{3};
state_reftend = state{4};

state_evoked.ecg_time       = {}; % timebins fo spectrogram
state_evoked.ecg_b2bt   = {}; % evoked LFP response
state_evoked.state_id   = state_id;
state_evoked.state_name = state_name;

% find mean Rpeak interval for all the trials for normalization
R2Rt = [trials_lfp.ECG_b2btime];
Rpeaks = [trials_lfp.ECG_spikes];
mean_ECG_b2btime = nanmean(R2Rt(Rpeaks));    

for t = 1:length(trials_lfp)
        
    if isempty(trials_lfp(t).ECG_b2btime) || any(isnan(trials_lfp(t).ECG_b2btime))
        continue;
    end

    states          = trials_lfp(t).states;
    state_onset_t   = states([states(:).id] == ...
        state_id).onset_t;
    state_start_t   = states([states(:).id] == ...
        state_id).onset_t + state_reftstart;
    state_end_t     = states([states(:).id] == ...
        state_id).onset_t + state_reftend;

    ts = trials_lfp(t).tsample;
    state_onset_sample = round(state_onset_t/ts) + 1;
    state_start_sample = state_onset_sample + round(state_reftstart/ts);
    state_end_sample = state_onset_sample + round(state_reftend/ts);
    % timestamps
    %nsamples = round((state_reftend - state_reftstart) / ts);
    nsamples = state_end_sample - state_start_sample + 1;
%     if ~mod(nsamples, 2)
%         nsamples = nsamples + 1;
%     end
    %state_end_sample = state_start_sample + nsamples - 1;    
    
    trial_ecg_time = linspace(state_reftstart, state_reftend, nsamples);
    trial_ecg_time = trial_ecg_time - ...
        trial_ecg_time(abs(trial_ecg_time) == ...
        min(abs(trial_ecg_time)));
    % raw ECG
    trial_evoked_ecg_b2bt = nan(1, nsamples);
    state_start_sample = max(1, state_start_sample);
    state_end_sample = min(length(trials_lfp(t).ECG_b2btime), ...
        state_end_sample);
    ref_onset_sample = find(trial_ecg_time == 0);
    trial_evoked_ecg_b2bt(ref_onset_sample - (state_onset_sample - state_start_sample):...
        ref_onset_sample + (state_end_sample - state_onset_sample)) = ...
        trials_lfp(t).ECG_b2btime(state_start_sample:state_end_sample);
    
%     trial_ecg_timeright = 0:ts:ts*floor(nsamples/2);
%     trial_ecg_timeleft = -ts:-ts:-ts*floor(nsamples/2);
%     trial_ecg_time = [flip(trial_ecg_timeleft) trial_ecg_timeright];
%     trial_ecg_time = trial_ecg_time - ...
%         trial_ecg_time(abs(trial_ecg_time) == ...
%         min(abs(trial_ecg_time)));
%     trial_ecg_b2bt = trials_lfp(t).ECG_b2btime(trial_samples);
%     
%     trial_evoked_ecg_b2bt(1:length(trial_samples)) = ...
%         trials_lfp(t).ECG_b2btime(trial_samples);
    
%     state_evoked_ecg_time = trials_lfp(t).time(...
%         (trials_lfp(t).time >= state_start_t & ...
%         trials_lfp(t).time <= state_end_t)) - state_onset_t;
%     % evoked LFP for this state
%     trial_evoked_ecg = trials_lfp(t).ecg_data(...
%         (trials_lfp(t).time >= state_start_t & ...
%         trials_lfp(t).time <= state_end_t));    
    
    % put onset timestamp to zero
%     onset_timestamp = state_evoked_ecg_time(...
%         abs(state_evoked_ecg_time) == min(abs(state_evoked_ecg_time)));
%     state_evoked_ecg_time = state_evoked_ecg_time - ...
%         onset_timestamp;
   
    % evoked LFP for this state
    state_evoked.ecg_b2bt = [state_evoked.ecg_b2bt, ...
        trial_evoked_ecg_b2bt];
    state_evoked.ecg_time = [state_evoked.ecg_time, trial_ecg_time];
    
    

end

if ~isempty(state_evoked.ecg_b2bt)
    
    % remove nans
    arr_state_ecg_b2bt = cat(1, state_evoked.ecg_b2bt{:});
    arr_state_ecg_b2bt(:, isnan(sum(arr_state_ecg_b2bt, 1))) = nan;
%     state_evoked.ecg_time = ...
%         state_evoked_ecg_time;
    state_evoked.ecg_time = ...
        trial_ecg_time;%(~any(isinf(arr_state_ecg_b2bt), 1));
    %arr_state_ecg_b2bt(:, any(isinf(arr_state_ecg_b2bt), 1)) = [];
    state_evoked.ecg_b2bt = arr_state_ecg_b2bt;
    state_evoked.mean_ecg_b2bt = mean_ECG_b2btime;
    state_evoked.dimord = 'ntrials_time';

    % crop each lfp to same number of samples
%     min_ecg_time = state_evoked_ecg.ecg_time{cellfun('length', state_evoked_ecg.ecg_time) == ...
%         min(cellfun('length', state_evoked_ecg.ecg_time))};
%     onset_sample = find(min_ecg_time == 0);
%     ref_start_sample = 1 - onset_sample;
%     ref_end_sample = length(min_ecg_time) - onset_sample;
% 
%     for k = 1:length(state_evoked_ecg.ecg)
%         onset_sample = find(state_evoked_ecg.ecg_time{k} == 0);
%         state_evoked_ecg.ecg{k} = state_evoked_ecg.ecg{k}(onset_sample + ref_start_sample: ...
%             onset_sample  + ref_end_sample);        
%     end    
%     state_evoked_ecg.ecg_time= min_ecg_time; %state_evoked_ecg.ecg_time(1:nsamples);
    
    % crop each lfp to same number of samples
%     nsamples = min(cellfun('length', state_evoked_ecg.ecg));
%     for k = 1:length(state_evoked_ecg.ecg)
%         state_evoked_ecg.ecg{k} = state_evoked_ecg.ecg{k}(1:nsamples);
%     end
%     state_evoked_ecg.ecg_time = state_evoked_ecg.ecg_time(1:nsamples);
    
    % evoked LFP average
    state_evoked.mean = nanmean(arr_state_ecg_b2bt / mean_ECG_b2btime, 1);
    state_evoked.std = nanstd(arr_state_ecg_b2bt / mean_ECG_b2btime, 0, 1);
    
end