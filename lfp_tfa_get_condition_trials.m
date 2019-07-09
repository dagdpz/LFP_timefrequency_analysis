function cond_trials = lfp_tfa_get_condition_trials(site_lfp, condition)
%lfp_tfa_get_condition_trials - Function to get the indices fo trials which
%satisfy a condition
%
% USAGE:
%	cond_trials = lfp_tfa_get_condition_trials(site_lfp, condition)
%
% INPUTS:
%       site_lfp                  - struct containing the LFP data for all
%       trials for one site
%       condition                 - struct containing the condition to be
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
%		cond_trials               - array of length (1 x N, N = number of
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
    if ~isnan(condition.type)
        cond_trials = cond_trials & ...
            ([site_lfp.trials.type] == condition.type);
    end
    % filter by effector
    if ~isnan(condition.effector)
        cond_trials = cond_trials & ...
            ([site_lfp.trials.effector] == condition.effector);
    end
    % filter by choice
    if ~isnan(condition.choice)
        cond_trials = cond_trials & ...
            ([site_lfp.trials.choice_trial] == condition.choice);
    end
    
    % filter by perturbation
    
    % commented on 08.05.2019, to be tested
%     preinj_perturb = 0;
%     postinj_perturb = unique(perturbations(perturbations ~= 0));
    cond_trials_perturb = zeros(1, length(site_lfp.trials));
    if ~isnan(condition.perturbation)
%         if condition.perturbation == 0 % pre-injection
%             perturbation_values = 0;
%         else % post-injection
%             perturbation_values = unique(perturbations(perturbations ~= 0));
%         end
        perturbation_values = unique([site_lfp.trials.perturbation]);
        
        if strcmp(condition.perturbation_group, 'all')
            for b = perturbation_values
                cond_trials_perturb = cond_trials_perturb | ...
                ([site_lfp.trials.perturbation] == b);
            end
        elseif strcmp(condition.perturbation_group, 'allbutfirst')
            for b = perturbation_values(2:end)
                cond_trials_perturb = cond_trials_perturb | ...
                ([site_lfp.trials.perturbation] == b);
            end
        else 
            for b = condition.perturbation_group{1}
                cond_trials_perturb = cond_trials_perturb | ...
                ([site_lfp.trials.perturbation] == b);
            end
        end
    end
    if isnan(condition.perturbation) % ignore perturbation
        cond_trials_perturb = ones(1, length(site_lfp.trials));
    end
    cond_trials = cond_trials & cond_trials_perturb;
