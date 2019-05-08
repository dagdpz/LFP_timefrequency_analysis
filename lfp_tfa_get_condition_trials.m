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
    
    % filter by perturbation
    perturbations = [site_lfp.trials.perturbation];
    % commented on 08.05.2019, to be tested
%     preinj_perturb = 0;
%     postinj_perturb = unique(perturbations(perturbations ~= 0));
    cond_trials_perturb = zeros(1, length(site_lfp.trials));
    if ~isnan(condition.perturbation)
        if condition.perturbation == 0 % pre-injection
            perturbation_values = 0;
        else % post-injection
            perturbation_values = unique(perturbations(perturbations ~= 0));
        end
        
        if strcmp(condition.perturbation_group, 'all')
            for b = perturbation_values
                cond_trials_perturb = cond_trials_perturb | ...
                ([site_lfp.trials.perturbation] == b);
            end
        elseif strcmp(condition.perturbation_group, 'allbutone')
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
