function cond_trials = lfp_tfa_get_condition_trials(site_lfp, condition)
%lfp_tfa_get_condition_trials - Function to get the indices of trials which
%satisfy a given trial condition
%
% USAGE:
%	cond_trials = lfp_tfa_get_condition_trials(site_lfp, condition)
%
% INPUTS:
%       site_lfp    - struct containing the condition
%       information about all trials for one site
%       condition   - struct containing the condition to be
%       analysed
%           Required fields:
%               type                - (integer) trial type
%               effector            - (integer) trial effector
%               choice              - choice (1) or instructed (1) trial
%               perturbation        - pre-injection (0) or post-injection (1)
%               perturbation_groups - blocks to be analysed for the given
%               perturbation, can be integer, integer array, 'all',
%               'allbutone'
%
% OUTPUTS:
%		cond_trials     - array of length (1 x N, N = number of
%		trials), with ones at indices of trials belonging to the given
%		condition
%
% REQUIRES:
%
% See also settings/lfp_tfa_settings, lfp_tfa_compare_conditions,
% lfp_tfa_site_average_tfr, lfp_tfa_site_evoked_LFP,
% lfp_tfa_site_powspctrum
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

cond_trials = ones(1, length(site_lfp.trials));
% filter by type
if isfield(condition, 'type') && ~isnan(condition.type) && ~isinf(condition.type) && isfield(site_lfp.trials, 'type')
    cond_trials = cond_trials & ([site_lfp.trials.type] == condition.type);
end
% filter by effector
if isfield(condition, 'effector') && ~isnan(condition.effector) && ~isinf(condition.effector) && isfield(site_lfp.trials, 'effector')
    cond_trials = cond_trials & ([site_lfp.trials.effector] == condition.effector);
end
% filter by choice
if isfield(condition, 'choice') && ~isnan(condition.choice) && ~isinf(condition.choice) && isfield(site_lfp.trials, 'choice_trial')
    cond_trials = cond_trials & ([site_lfp.trials.choice_trial] == condition.choice);
end
% filter by success
if isfield(condition, 'success') && ~isnan(condition.success) && ~isinf(condition.success) && isfield(site_lfp.trials, 'success')
    cond_trials = cond_trials & ([site_lfp.trials.success] == condition.success);
end

% filter by perturbation
cond_trials_perturb = zeros(1, length(site_lfp.trials));
if isfield(condition, 'perturbation') && ~isnan(condition.perturbation) && ~isinf(condition.perturbation) && isfield(site_lfp.trials, 'perturbation')
    perturbation_values = unique([site_lfp.trials.perturbation]);
    if ~isnan(perturbation_values)
        if condition.perturbation == 1%post-injection
            if strcmp(condition.perturbation_group, 'all')
                perturbation_values = perturbation_values(perturbation_values ~= 0);
            elseif strcmp(condition.perturbation_group, 'allbutfirst')
                perturbation_values = perturbation_values(perturbation_values ~= 0);
                perturbation_values = perturbation_values(2:end);
            else
                perturbation_values = condition.perturbation_group{1};
            end
        elseif condition.perturbation == 0 % pre-injection
            perturbation_values = condition.perturbation_group{1};
        end
        for b = perturbation_values
            cond_trials_perturb = cond_trials_perturb | ([site_lfp.trials.perturbation] == b);
        end
        cond_trials = cond_trials & cond_trials_perturb;
    end
end

%remove trials for which one state timing is not defined (e.g saccade
%is not detected so onset is NA)
for tr = 1:length(site_lfp.trials)
    state_onset_values = [site_lfp.trials(tr).states.onset_t];
    cond_trials_missing_timing(tr) = ~any(isnan(state_onset_values));
end
cond_trials = cond_trials & cond_trials_missing_timing;

end