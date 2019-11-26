function Rpeak_evoked_state = lfp_tfa_get_Rpeak_evoked_state_onsets(trials_ecg, state)

state_id = state{1};
state_name = state{2};
Rpeak_ref_abs_time = state{3};
Rpeak_ref_rel_time = state{4};%state_reftend = state{4};

%state_evoked_ecg.ecg_time = {}; % timebins fo spectrogram
Rpeak_evoked_state.abs_onset_times = [];
Rpeak_evoked_state.rel_onset_times = [];
Rpeak_evoked_state.valid_trials = [];
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
    
    if ~all(trials_ecg(t).ECG_valid)%(state_onset_sample)
        continue;
    end
    
    % Rpeaks
    trial_Rpeaks = trials_ecg(t).ECG_spikes;
    trial_timestamps = trials_ecg(t).time;
    Rpeak_times = trial_timestamps(trial_Rpeaks);
    Rpeak_ref_onset_times = state_onset_t - Rpeak_times;
        
    Rpeak_rel_onset_time = nan;
    % check if state onset occured after first Rpeak of this trial
    abs_onset_time_after_Rpeak = ...
            min(Rpeak_ref_onset_times(Rpeak_ref_onset_times > 0));
    if ~isempty(abs_onset_time_after_Rpeak)
        % find the id of Rpeak after which state onset occured
        Rpeak_idx = find(Rpeak_ref_onset_times == ...
            abs_onset_time_after_Rpeak);
        if Rpeak_idx < length(Rpeak_ref_onset_times)
            Rpeak_rel_onset_time = abs_onset_time_after_Rpeak / ...
                (Rpeak_times(Rpeak_idx + 1) - Rpeak_times(Rpeak_idx));
        end
    end
    
    %abs_onset_time_after_Rpeak = nan;
    %if strcmp(Rpeak_ref_rel_time, 'aroundRpeak')
    abs_onset_time_around_Rpeak = ...
        Rpeak_ref_onset_times(...
        Rpeak_ref_onset_times < Rpeak_ref_abs_time(2) & ...
        Rpeak_ref_onset_times > Rpeak_ref_abs_time(1));
%     abs_onset_time_before_Rpeak = ...
%         max(Rpeak_ref_onset_times(Rpeak_ref_onset_times < 0 & ...
%         Rpeak_ref_onset_times > Rpeak_ref_abs_time(1)));
%     abs_onset_time_around_Rpeak = [abs_onset_time_before_Rpeak, ...
%         abs_onset_time_after_Rpeak];   
        
%     elseif strcmp(Rpeak_ref_rel_time, 'beforeRpeak')
%         Rpeak_ref_onset_times = Rpeak_times - state_onset_t;
%         abs_onset_time_after_Rpeak = ...
%             max(abs_onset_time_after_Rpeak(abs_onset_time_after_Rpeak < 0));
%         % check if state onset occured before last Rpeak of this trial
%         if ~isempty(abs_onset_time_after_Rpeak)
%             % find the id of Rpeak before which state onset occured
%             Rpeak_idx = find(Rpeak_ref_onset_times == ...
%                 abs_onset_time_after_Rpeak);
%             if Rpeak_idx > 1
%                 Rpeak_rel_onset_time = abs_onset_time_after_Rpeak / ...
%                     (Rpeak_times(Rpeak_idx) - Rpeak_times(Rpeak_idx - 1));
%             end
%         end
%     elseif strcmp(Rpeak_ref_rel_time, 'mindist')
%         abs_onset_time_after_Rpeak = Rpeak_times - state_onset_t;
%         abs_onset_time_after_Rpeak = abs_onset_time_after_Rpeak(abs(abs_onset_time_after_Rpeak) == ...
%             min(abs(abs_onset_time_after_Rpeak))); 
%     end
    
    if ~isempty(abs_onset_time_after_Rpeak)
        Rpeak_evoked_state.rel_onset_times = ...
            [Rpeak_evoked_state.rel_onset_times Rpeak_rel_onset_time];
        Rpeak_evoked_state.valid_trials = [Rpeak_evoked_state.valid_trials, t];
    end
    
    if ~isempty(abs_onset_time_around_Rpeak)
        Rpeak_evoked_state.abs_onset_times = ...
            [Rpeak_evoked_state.abs_onset_times abs_onset_time_around_Rpeak];
    end
end

% get histogram counts
% absolute time from Rpeak
abs_tbinwidth = 0.025;
abs_tbinedges = Rpeak_ref_abs_time(1) : abs_tbinwidth : Rpeak_ref_abs_time(2);
[Rpeak_evoked_state.abs_histcounts.prob, ...
    Rpeak_evoked_state.abs_histcounts.timebins] = ...
    histcounts(Rpeak_evoked_state.abs_onset_times, abs_tbinedges, ...
    'Normalization', 'probability');
% relative time from Rpeak
rel_tbinwidth = 1/36;
rel_tbinedges = 0:rel_tbinwidth:1;
[Rpeak_evoked_state.rel_histcounts.prob, ...
    Rpeak_evoked_state.rel_histcounts.timebins] = ...
    histcounts(Rpeak_evoked_state.rel_onset_times, rel_tbinedges, ...
    'Normalization', 'probability');

Rpeak_evoked_state.ntrials = length(Rpeak_evoked_state.rel_onset_times);

end