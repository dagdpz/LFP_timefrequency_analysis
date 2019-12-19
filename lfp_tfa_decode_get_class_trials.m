function class_trials_idx = lfp_tfa_decode_get_class_trials(site_lfp, class_condition)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
class_trials_idx = true(1, length(site_lfp.trials));
% first remove noisy trials
class_trials_idx = class_trials_idx & ~[site_lfp.trials.noisy];
% get conditions to be checked
fields = fieldnames(class_condition);
for f = 1:length(fields)
    field = fields{f};
    if ~strcmp(field, 'perturbation') && ~strcmp(field, 'hndspc_lbl') ...
            && ~strcmp(field, 'reach_hand') && ~strcmp(field, 'reach_space')
        % get trials based on each condition
        if isfield(site_lfp.trials, field) && any(~isinf(class_condition.(field)))
            class_trials_idx = class_trials_idx & ...
                ismember([site_lfp.trials.(field)], class_condition.(field));
        end
    elseif strcmp(field, 'perturbation')
        % get trials based on perturbation
        if isfield(site_lfp.trials, field) && any(~isinf(class_condition.(field)))
            if class_condition.(field) == 0 
                class_trials_idx = class_trials_idx & [site_lfp.trials.(field)] == 0;
            else
                class_trials_idx = class_trials_idx & [site_lfp.trials.(field)] ~= 0;
            end
        end
    elseif strcmp(field, 'hndspc_lbl') || strcmp(field, 'reach_hand') || strcmp(field, 'reach_space')
        % get trials based on each condition
        if isfield(site_lfp.trials, field) && any(~isinf(class_condition.(field)))
            class_trials_idx = class_trials_idx & ...
                strcmp({site_lfp.trials.(field)}, class_condition.(field));
        end
    end
end
% get trial indices
class_trials_idx = find(class_trials_idx);

% % get trials based on effector
% if isfield(class_condition, 'effector') && isfield(site_lfp.trial, 'effector')
%     class_trials_idx = class_trials_idx & ...
%         ismember([site_lfp.trial.effector], class_condition.effector);
% end
% % get trials based on perturbation
% if isfield(class_condition, 'effector') && isfield(site_lfp.trial, 'effector')
%     class_trials_idx = class_trials_idx & ...
%         ismember([site_lfp.trial.effector], class_condition.effector);
% end
end

