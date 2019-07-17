function state_evoked_lfp = lfp_tfa_get_state_evoked_lfp(trials_lfp, state)

state_id = state{2};
state_name = state{3};
state_reftstart = state{4};
state_reftend = state{5};

state_evoked_lfp.time = {}; % timebins fo spectrogram
state_evoked_lfp.freq = {}; % freq bins
state_evoked_lfp.lfp = {}; % evoked LFP response
state_evoked_lfp.state_id = state_id;
state_evoked_lfp.state_name = state_name;

for t = 1:length(trials_lfp)

    states          = trials_lfp(t).states;
    state_onset_t   = states([states(:).id] == ...
        state_id).onset_t;
    state_start_t   = states([states(:).id] == ...
        state_id).onset_t + state_reftstart;
    state_end_t     = states([states(:).id] == ...
        state_id).onset_t + state_reftend;

    % evoked LFP for this state
    state_evoked_lfp.lfp = [state_evoked_lfp.lfp, ...
        trials_lfp(t).lfp_data(...
        (trials_lfp(t).time >= state_start_t & ...
        trials_lfp(t).time <= state_end_t))];
    % timestamps
    state_evoked_lfp.lfp_time = trials_lfp(t).time(...
        (trials_lfp(t).time >= state_start_t & ...
        trials_lfp(t).time <= state_end_t)) - state_onset_t;
    % put onset timestamp to zero
    onset_timestamp = state_evoked_lfp.lfp_time(...
        abs(state_evoked_lfp.lfp_time) == min(abs(state_evoked_lfp.lfp_time)));
    state_evoked_lfp.lfp_time = state_evoked_lfp.lfp_time - ...
        onset_timestamp;

end


if ~isempty(state_evoked_lfp.lfp)

    % crop each lfp to same number of samples
    nsamples = min(cellfun('length', state_evoked_lfp.lfp));
    for k = 1:length(state_evoked_lfp.lfp)
        state_evoked_lfp.lfp{k} = state_evoked_lfp.lfp{k}(1:nsamples);
    end
    state_evoked_lfp.lfp_time = state_evoked_lfp.lfp_time(1:nsamples);
    
    % evoked LFP average
    arr_state_lfp = vertcat(state_evoked_lfp.lfp{:});
    state_evoked_lfp.mean = nanmean(arr_state_lfp, 1);
    state_evoked_lfp.std = nanstd(arr_state_lfp, 0, 1);
    
end