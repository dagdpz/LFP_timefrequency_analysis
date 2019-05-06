function [ cmp_conditions ] = lfp_tfa_compare_conditions( lfp_tfa_cfg, varargin )
%lfp_tfa_compare_conditions  - Create conditions to compare by permuting
%the compare conditions in the settings
%
% USAGE:
%	cmp_conditions = lfp_tfa_compare_conditions( lfp_tfa_cfg )
%
% INPUTS:
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields: see lfp_tfa_settings
%               1. compare.type                 - trial types to be compared
%               2. compare.effector             - trial effectors to be compared
%               3. compare.choice               - trial choices to be compared
%               (0 = instructed, 1 = choice trial)
%               4. compare.perturbation         - perturbations to be compared
%               (0 = preinjection, 1 = postinjection)
%               5. compare.perturbation_groups  - perturbation_groups for
%               pre and post injection (typically 0 for pre and same as
%               block number for post)
%               6. compare.reach_hands          - hand labels to compare
%               ('R' = right, 'L' = left)
%               7. compare.reach_spaces         - space labels to compare
%               ('R' = right, 'L' = left)
%               8. ref_hemisphere               - reference hemisphere for
%               hand-space labelling ('R' or 'L', typically, the inactivated
%               hemisphere)
% OUTPUTS:
%		cmp_conditions      - structure containing conditions to compare, a
%		permutation of all given comparison conditions
%
% REQUIRES:	
%
% See also lfp_tfa_settings, lfp_tfa_define_settings
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

    task_types = lfp_tfa_cfg.compare.types;
    effectors = lfp_tfa_cfg.compare.effectors;
    targets = lfp_tfa_cfg.compare.targets;
    if nargin > 1
        targets = varargin{1};
    end
    %unique([states_lfp.recorded_hemispace]);
    if lfp_tfa_cfg.compare.choice_trials
        choices = unique([states_lfp(1).trials.choice_trial]);
    else
        choices = lfp_tfa_cfg.compare.choice_trials;
    end
    perturbations = lfp_tfa_cfg.compare.perturbations;
    perturbation_groups = lfp_tfa_cfg.compare.perturbation_groups;
    
    hands = lfp_tfa_cfg.compare.reach_hands;
    spaces = lfp_tfa_cfg.compare.reach_spaces;  
    % assign hand space labels
    hs_labels = cell(1, length(hands)*length(spaces));
    reach_hands = cell(1, length(hands)*length(spaces));
    reach_spaces = cell(1, length(hands)*length(spaces));
    for h = 1:length(hands)
        if strcmp(hands{h},'R') || strcmp(hands{h},'L')
            if strcmp(hands{h},lfp_tfa_cfg.ref_hemisphere)
                hand_label = 'IH';
            else
                hand_label = 'CH';
            end
        else
            hand_label = [hands{h}, 'H'];
        end
        for s = 1:length(spaces)
            if strcmp(spaces{s},'R') || strcmp(spaces{s},'L')
                if strcmp(spaces{s},lfp_tfa_cfg.ref_hemisphere)
                    space_label = 'IS';
                else
                    space_label = 'CS';
                end
            else
                space_label = [spaces{s}, 'S'];
            end
            hs_idx = (h-1)*length(spaces) + s;
            reach_hands{hs_idx} = hands{h};
            reach_spaces{hs_idx} = spaces{s};
            hs_labels{hs_idx} = [hand_label ' ' space_label];
        end
    end
    
    % create conditions
    cmp_conditions = struct();
    
    i = 0;
    %for target = targets
        %target_label = target{1};
    for type = task_types
        type_label = ['Type_' num2str(type)];
        for eff = effectors
            eff_label = ['Eff_' num2str(eff)];            
            for ch = choices
                if ch == 0
                    ch_label = 'Instr';
                elseif choice == 1
                    ch_label = 'Choice';
                else
                    ch_label = [];
                end
                for p = 1:length(perturbations)
                    if perturbations(p) == 0
                        p_label = 'Pre';
                    elseif perturbations(p) == 1
                        p_label = 'Post';
                    else
                        p_label = [];
                    end 

                    i = i + 1;
                    condition_label = [type_label, '_', eff_label, '_', ...
                        ch_label, '_', p_label];
                    cmp_conditions(i).type = type;
                    cmp_conditions(i).effector = eff;
                    %cmp_conditions(i).target = target{1};
                    cmp_conditions(i).choice = ch;
                    cmp_conditions(i).perturbation = perturbations(p);
                    cmp_conditions(i).perturbation_group = perturbation_groups(p);
                    cmp_conditions(i).hs_labels = hs_labels;
                    cmp_conditions(i).reach_hands = reach_hands;
                    cmp_conditions(i).reach_spaces = reach_spaces;
                    cmp_conditions(i).label = condition_label;
%                             end
%                         end
                end
            end
        end
    end
    % end
    

end

