function Rpeak_evoked_state = lfp_tfa_get_Rpeak_evoked_state_onsets(trials_ecg, state)

state_id = state{1};
state_name = state{2};
Rpeak_reftime = state{3};
%state_reftend = state{4};

%state_evoked_ecg.ecg_time = {}; % timebins fo spectrogram
Rpeak_evoked_state.abs_onset_times = ...
    nan(1, length(trials_ecg)); % evoked LFP response
Rpeak_evoked_state.rel_onset_times = ...
    nan(1, length(trials_ecg)); 
Rpeak_evoked_state.state_id = state_id;
Rpeak_evoked_state.state_name = state_name;

for t = 1:length(trials_ecg)
    
    if isempty(trials_ecg(t).ECG_spikes) || ...
            isempty(trials_ecg(t).states) 
        continue;
    end

    states          = trials_ecg(t).states;
    if ismember(state_id, [states(:).id])
        state_onset_t   = states([states(:).id] == ...
            state_id).onset_t;
        state_onset_sample   = states([states(:).id] == ...
            state_id).onset_s;
    else
        continue;
    end
    
    if ~trials_ecg(t).ECG_valid(state_onset_sample)
        continue;
    end
        
    % Rpeaks
    trial_Rpeaks = trials_ecg(t).ECG_spikes;
    trial_timestamps = trials_ecg(t).time;
    Rpeak_times = trial_timestamps(trial_Rpeaks);
    Rpeak_abs_onset_time = nan;
    Rpeak_rel_onset_time = nan;
    if strcmp(Rpeak_reftime, 'afterRpeak')
        Rpeak_ref_onset_times = state_onset_t - Rpeak_times;
        Rpeak_abs_onset_time = ...
            min(Rpeak_ref_onset_times(Rpeak_ref_onset_times > 0));
        % check if state onset occured after first Rpeak of this trial
        if ~isempty(Rpeak_abs_onset_time)
            % find the id of Rpeak after which state onset occured
            Rpeak_idx = find(Rpeak_ref_onset_times == ...
                Rpeak_abs_onset_time);
            if Rpeak_idx < length(Rpeak_ref_onset_times)
                Rpeak_rel_onset_time = Rpeak_abs_onset_time / ...
                    (Rpeak_times(Rpeak_idx + 1) - Rpeak_times(Rpeak_idx));
            end
        end
    elseif strcmp(Rpeak_reftime, 'beforeRpeak')
        Rpeak_ref_onset_times = Rpeak_times - state_onset_t;
        Rpeak_abs_onset_time = ...
            max(Rpeak_abs_onset_time(Rpeak_abs_onset_time < 0));
        % check if state onset occured before last Rpeak of this trial
        if ~isempty(Rpeak_abs_onset_time)
            % find the id of Rpeak before which state onset occured
            Rpeak_idx = find(Rpeak_ref_onset_times == ...
                Rpeak_abs_onset_time);
            if Rpeak_idx > 1
                Rpeak_rel_onset_time = Rpeak_abs_onset_time / ...
                    (Rpeak_times(Rpeak_idx) - Rpeak_times(Rpeak_idx - 1));
            end
        end
    elseif strcmp(Rpeak_reftime, 'mindist')
        Rpeak_abs_onset_time = Rpeak_times - state_onset_t;
        Rpeak_abs_onset_time = Rpeak_abs_onset_time(abs(Rpeak_abs_onset_time) == ...
            min(abs(Rpeak_abs_onset_time))); 
    end
    
    if ~isempty(Rpeak_abs_onset_time)
        Rpeak_evoked_state.abs_onset_times(t) = Rpeak_abs_onset_time;
        Rpeak_evoked_state.rel_onset_times(t) = Rpeak_rel_onset_time;
    end
    
end

% get histogram counts
% absolute time from Rpeak
[Rpeak_evoked_state.abs_histcounts.prob, Rpeak_evoked_state.abs_histcounts.timebins] = ...
    histcounts(Rpeak_evoked_state.abs_onset_times, 'BinWidth', 0.05, 'Normalization', 'probability');
% relative time from Rpeak
[Rpeak_evoked_state.rel_histcounts.prob, Rpeak_evoked_state.rel_histcounts.timebins] = ...
    histcounts(Rpeak_evoked_state.rel_onset_times, 'BinWidth', 0.1, 'Normalization', 'probability');

end