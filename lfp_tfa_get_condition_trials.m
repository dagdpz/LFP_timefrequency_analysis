function cond_trials = lfp_tfa_get_condition_trials(site_lfp, condition)

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
    % filter by hand
%                     if ~isnan(cfg_conditions(cn).reach_hand)
%                         cond_trials = cond_trials & ...
%                             ([states_lfp(i).trials.reach_hand] == ...
%                             cfg_conditions(cn).reach_hand);
%                     end
%                     % filter by space
%                     if ~isnan(cfg_conditions(cn).reach_space)
%                         cond_trials = cond_trials & ...
%                             ([states_lfp(i).trials.reach_space] == ...
%                             cfg_conditions(cn).reach_space);
%                     end
    % filter by perturbation
    perturbations = [site_lfp.trials.perturbation];
    postinj_perturb = unique(perturbations(perturbations ~= 0));
    cond_trials_perturb = zeros(1, length(site_lfp.trials));
    if ~isnan(condition.perturbation)
        
        if strcmp(condition.perturbation_group, 'all')
            for b = postinj_perturb
                cond_trials_perturb = cond_trials_perturb | ...
                ([site_lfp.trials.perturbation] == b);
            end
        elseif strcmp(condition.perturbation_group, 'allbutone')
            for b = postinj_perturb(2:end)
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
    
    cond_trials = cond_trials & cond_trials_perturb;
%                     if ~isnan(cfg_conditions(cn).block)
%                         for b = cfg_conditions(cn).block
%                             cond_trials = cond_trials | ...
%                                 ([states_lfp(i).trials.block] == b);
%                         end
%                     end