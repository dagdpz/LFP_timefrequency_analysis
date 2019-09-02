function Rpeak_evoked_state = lfp_tfa_get_Rpeak_evoked_state_onsets(trials_ecg, state)

state_id = state{1};
state_name = state{2};
Rpeak_reftime = state{3};
%state_reftend = state{4};

%state_evoked_ecg.ecg_time = {}; % timebins fo spectrogram
Rpeak_evoked_state.ref_onset_times = ...
    nan(1, length(trials_ecg)); % evoked LFP response
Rpeak_evoked_state.state_id = state_id;
Rpeak_evoked_state.state_name = state_name;

for t = 1:length(trials_ecg)
    
    if isempty(trials_ecg(t).ECG_spikes) || ...
            isempty(trials_ecg(t).states)
        continue;
    end

    states          = trials_ecg(t).states;
    state_onset_t   = states([states(:).id] == ...
        state_id).onset_t;
    
    % Rpeaks
    trial_Rpeaks = trials_ecg(t).ECG_spikes;
    trial_timestamps = trials_ecg(t).time;
    Rpeak_times = trial_timestamps(trial_Rpeaks);
    if strcmp(Rpeak_reftime, 'afterRpeak')
        Rpeak_ref_onset_time = Rpeak_times - state_onset_t;
        Rpeak_ref_onset_time = min(Rpeak_ref_onset_time(Rpeak_ref_onset_time > 0));
    elseif strcmp(Rpeak_reftime, 'beforeRpeak')
        Rpeak_ref_onset_time = Rpeak_times - state_onset_t;
        Rpeak_ref_onset_time = max(Rpeak_ref_onset_time(Rpeak_ref_onset_time < 0)); 
    elseif strcmp(Rpeak_reftime, 'mindist')
        Rpeak_ref_onset_time = Rpeak_times - state_onset_t;
        Rpeak_ref_onset_time = Rpeak_ref_onset_time(abs(Rpeak_ref_onset_time) == ...
            min(abs(Rpeak_ref_onset_time))); 
    end
    
    if ~isempty(Rpeak_ref_onset_time)
        Rpeak_evoked_state.ref_onset_times(t) = Rpeak_ref_onset_time;
    end
    
end