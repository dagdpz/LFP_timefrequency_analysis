function state_evoked_lfp = lfp_tfa_get_state_evoked_lfp(trials_lfp, state)
%lfp_tfa_get_state_evoked_lfp - Function to get the evoked LFP average 
%across trials in a time window of interest around a givem state (event)
%
% USAGE:
%	state_evoked_lfp = lfp_tfa_get_state_evoked_lfp(trials_lfp, state)
%
% INPUTS:
%       trials_lfp    - struct containing LFP signal for different trials
%       (which usually satisfy a trial condition, see
%       lfp_tfa_plot_site_evoked_LFP and lfp_tfa_get_condition_trials)
%       state         - 1x5 cell array containing the information about the 
%       state (event) which should be analysed, see
%       settings/lfp_tfa_settings_example (same as one row of
%       lfp_tfa_cfg.analyse_states)
%
%
% OUTPUTS:
%		state_evoked_lfp     - structure containing evoked LFP average 
%       across trials in a time window of interest around a givem state
%
% REQUIRES:	
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_get_condition_trials,
% lfp_tfa_site_evoked_LFP
%
% Author(s):	S.Nair, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2019-02-15:	Created function (Sarath Nair)
% 2019-03-05:	First Revision
% ...
% $Revision: 1.0 $  $Date: 2019-03-05 17:18:00 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%

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
    if ~ismember(state_id,[states(:).id])
        %% LS 2021: replace missing state alignment with NaNs
        
        state_evoked_lfp.lfp = [state_evoked_lfp.lfp, ...
            NaN(1,length(trials_lfp(t).tfs.freq),floor((state_reftend-state_reftstart)/(trials_lfp(t).time(2)-trials_lfp(t).time(1))))];
        state_evoked_lfp.lfp_time=NaN(1,floor((state_reftend-state_reftstart)/(trials_lfp(t).time(2)-trials_lfp(t).time(1))));       
        continue;
    end
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
        onset_timestamp(1); %% LS2021 funny error here: +/- the exact same number when running pulv_oculomotor

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
    state_evoked_lfp.lfp = arr_state_lfp;
    %state_evoked_lfp.std = nanstd(arr_state_lfp, 0, 1);
    
end