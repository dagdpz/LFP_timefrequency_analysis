function state_tfs = lfp_tfa_get_state_tfs(site_lfp, cond_trials, state, lfp_tfa_cfg, perturbation)
%lfp_tfa_get_state_tfs - Function to get the LFP power spectrogram average 
%across trials in a time window of interest around a givem state (event)
%
% USAGE:
%	state_tfs = lfp_tfa_get_state_tfs(site_lfp, cond_trials, state, lfp_tfa_cfg)
%
% INPUTS:
%       site_lfp      - struct containing LFP signal for all trials
%       recorded from a single site during one session
%       cond_trials   - 1xN logical array where ones represent the indices
%       of the trials which satisfy a given condition, see lfp_tfa_get_condition_trials
%       state         - 1x5 cell array containing the information about the 
%       state (event) which should be analysed, see
%       settings/lfp_tfa_settings_example (same as one row of
%       lfp_tfa_cfg.analyse_states)
%       lfp_tfa_cfg   - settings for baseline normalization of LFP power,
%       see settings/lfp_tfa_settings_example
%           Required fields:
%               baseline_method     - method used for baseline
%               normalization ('subtraction', 'division', 'zscore' or
%               'relchange'), see lfp_tfa_baseline_normalization
%               baseline_perturbation - which perturbation condition should
%               be used for baseline power calculation (0 for pre and 1 for
%               post injection)
%               baseline_use_choice_trial - whether to use choice or
%               instructed trials for baseline power computation (0 for
%               instructed, 1 for choice)
% OUTPUTS:
%		state_tfs     - structure containing LFP power spectrogram average 
%       across given trials in a time window of interest around a givem state
%
% REQUIRES:	lfp_tfa_baseline_normalization
%
% See also settings/lfp_tfa_settings_example,
% lfp_tfa_baseline_normalization, lfp_tfa_get_condition_trials,
% lfp_tfa_plot_site_average_tfr
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
    if ~ismember(state_id,[states(:).id])
        %% LS 2021: replace missing state alignment with NaNs
        
        state_tfs.powspctrm = [state_tfs.powspctrm, ...
            NaN(1,length(lfp_tfa_cfg.tfr.foi),floor((state_ref_tend-state_ref_tstart)/lfp_tfa_cfg.tfr.timestep*site_lfp.trials(t).fsample))];
        state_tfs.time=NaN(1,floor((state_ref_tend-state_ref_tstart)/lfp_tfa_cfg.tfr.timestep*site_lfp.trials(t).fsample));
        state_tfs.freq = site_lfp.trials(t).tfs.freq;
        state_tfs.cfg = site_lfp.trials(t).tfs.cfg;
        continue;
    else
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
    state_tfs.time = state_tfs.time - onset_timestamp(1); % added (1) IK 20210819
    % freq bins
    state_tfs.freq = site_lfp.trials(t).tfs.freq; 
    state_tfs.cfg = site_lfp.trials(t).tfs.cfg;
    end

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
    if  lfp_tfa_cfg.baseline_perturbation == -1
        baseline_cnd_idx = [site_lfp.baseline.perturbation] == ...
        perturbation & [site_lfp.baseline.choice] == ...
        lfp_tfa_cfg.baseline_use_choice_trial;
    else
    baseline_cnd_idx = [site_lfp.baseline.perturbation] == ...
        lfp_tfa_cfg.baseline_perturbation & [site_lfp.baseline.choice] == ...
        lfp_tfa_cfg.baseline_use_choice_trial;
    end
    cfg_baseline.mean = site_lfp.baseline(baseline_cnd_idx).pow_mean;
    cfg_baseline.std = site_lfp.baseline(baseline_cnd_idx).pow_std;
    state_tfs.powspctrm = lfp_tfa_baseline_normalization(...
        state_tfs.powspctrm, cfg_baseline); 
    state_tfs.baseline = cfg_baseline;
    
end