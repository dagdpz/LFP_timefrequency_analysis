function class_condition_idx = lfp_tfa_decode_get_class_condition(conditions,class)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

class_condition_idx = true(1, length(conditions));
% get conditions to be checked
fields = fieldnames(class);
for f = 1:length(fields)
    field = fields{f};
    if isempty(class.(field))
        continue;
    end
    if ~strcmp(field, 'label') && ~strcmp(field, 'perturbation') && ~strcmp(field, 'hs_label') ...
            && ~strcmp(field, 'reach_hand') && ~strcmp(field, 'reach_space') && ~strcmp(field, 'choice_trial')
        % get trials based on each condition
        if isfield(conditions, field) && any(~isinf(class.(field)))
            class_condition_idx = class_condition_idx & ...
                ismember([conditions.(field)], class.(field));
        end
    elseif strcmp(field, 'perturbation')
        % get trials based on perturbation
        if isfield(conditions, field) && any(~isinf(class.(field)))
            if class.(field) == 0 
                class_condition_idx = class_condition_idx & [conditions.(field)] == 0;
            else
                class_condition_idx = class_condition_idx & [conditions.(field)] ~= 0;
            end
        end
    elseif strcmp(field, 'hs_label') || strcmp(field, 'reach_hand') || strcmp(field, 'reach_space')
        % get trials based on each condition
        if isfield(conditions, field) && any(~isinf(class.(field)))
            class_condition_idx = class_condition_idx & ...
                strcmp({conditions.(field)}, class.(field));
        end
    elseif strcmp(field, 'choice_trial')
        % get trials based on each condition
        if isfield(conditions, 'choice') && any(~isinf(class.(field)))
            class_condition_idx = class_condition_idx & ...
                ismember([conditions.choice], class.(field));
        end
    end
end

class_condition_idx = find(class_condition_idx);


end

