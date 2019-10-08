function state_evoked_ecg = lfp_tfa_get_state_evoked_ECG(trials_lfp, state)

state_id = state{1};
state_name = state{2};
state_reftstart = state{3};
state_reftend = state{4};

state_evoked_ecg.ecg_time = {}; % timebins fo spectrogram
state_evoked_ecg.ecg = {}; % evoked LFP response
state_evoked_ecg.state_id = state_id;
state_evoked_ecg.state_name = state_name;

for t = 1:length(trials_lfp)
    
    if isempty(trials_lfp(t).ecg_data)
        continue;
    end

    states          = trials_lfp(t).states;
    if ismember(state_id, [states(:).id])
        state_onset_t   = states([states(:).id] == ...
            state_id).onset_t;
    else
        continue;
    end
    state_start_t   = states([states(:).id] == ...
        state_id).onset_t + state_reftstart;
    state_end_t     = states([states(:).id] == ...
        state_id).onset_t + state_reftend;

    ts = trials_lfp(t).tsample;
    % timestamps
    state_evoked_ecg_time = state_reftstart-ts:ts:state_reftend+ts;
    state_evoked_ecg_time = round(state_evoked_ecg_time - ...
        state_evoked_ecg_time(abs(state_evoked_ecg_time) == ...
        min(abs(state_evoked_ecg_time))), 4);
    state_evoked_ecg_time = state_evoked_ecg_time(...
        state_evoked_ecg_time >= state_reftstart & ...
        state_evoked_ecg_time <= state_reftend);
    % raw ECG
    trial_evoked_ecg = nan(size(state_evoked_ecg_time));
    trial_time_idx = trials_lfp(t).time >= state_start_t & ...
        trials_lfp(t).time <= state_end_t;
    trial_ecg_time = trials_lfp(t).time(trial_time_idx) - state_onset_t;
    trial_ecg_time = round(trial_ecg_time - ...
        trial_ecg_time(abs(trial_ecg_time) == ...
        min(abs(trial_ecg_time))), 4);
    trial_evoked_ecg(state_evoked_ecg_time >= trial_ecg_time(1) & ...
        state_evoked_ecg_time <= trial_ecg_time(end)) = ...
        trials_lfp(t).ecg_data(trial_time_idx);
    
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
    state_evoked_ecg.ecg = [state_evoked_ecg.ecg, ...
        trial_evoked_ecg];
    state_evoked_ecg.ecg_time = [state_evoked_ecg.ecg_time, state_evoked_ecg_time];
    
    

end

if ~isempty(state_evoked_ecg.ecg)
    
    % remove nans
    arr_state_ecg = cat(1, state_evoked_ecg.ecg{:});
    state_evoked_ecg.ecg_time = ...
        state_evoked_ecg_time;%(~any(isnan(arr_state_ecg), 1));
    %arr_state_ecg(:, any(isnan(arr_state_ecg), 1)) = [];
    state_evoked_ecg.ecg = arr_state_ecg;
    state_evoked_ecg.dimord = 'nbeats_time';

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
    state_evoked_ecg.mean = mean(arr_state_ecg, 1);
    state_evoked_ecg.std = std(arr_state_ecg, 0, 1);
    
end