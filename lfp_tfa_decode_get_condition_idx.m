function [cnd_idx, st_idx, hs_idx] = lfp_tfa_decode_get_condition_idx(site_evoked, class_condition)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here    
    trial_condition = [site_evoked.condition.cfg_condition];
    class_condition_fieldnames = fieldnames(class_condition)';
    cnd_idx = true(1, length(trial_condition)); 
    % find condition index 
    for f = 1:length(class_condition_fieldnames)
        field = class_condition_fieldnames{f};
        if ~strcmp(field, 'label')
            if isfield(class_condition, field) && ...
                    isfield(trial_condition, field)
                cnd_idx = cnd_idx & ...
                    ismember([trial_condition.(field)], class_condition.(field));
            end
        end
    end
    cnd_idx = find(cnd_idx);
    % find state and handspace index
    hs_tuned_evoked = site_evoked.condition(cnd_idx).hs_tuned_evoked;
    states = [hs_tuned_evoked(:,1).state];
    st_idx = true(1, length(states));
    if isfield(class_condition, 'state')
        st_idx = find(states == class_condition.state);
    end
    st_idx = find(st_idx);
    hs_labels = [hs_tuned_evoked(1, :).hs_label];
    hs_idx = true(1, length(hs_labels));
    if isfield(class_condition, 'hs_label')
        hs_idx = (ismember(hs_labels, class_condition.hs_label));
    end
    hs_idx = find(hs_idx);
    
end

