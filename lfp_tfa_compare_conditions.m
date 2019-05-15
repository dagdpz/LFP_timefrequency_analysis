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
    
    choices = lfp_tfa_cfg.compare.choice_trials;
    perturbations = lfp_tfa_cfg.compare.perturbations;
    % if different sessions have different perturbation groups
    if nargin > 1
        perturbation_groups = varargin{1};

        % default perturbation groups based on perturbation input
    elseif sum (perturbations == [0, 1])
        perturbation_groups = {0, 'all'};
    elseif perturbations == 0
        perturbation_groups = {0};
    end
    % commented on 08052019, each session has its own perturbation group
    %perturbation_groups = lfp_tfa_cfg.compare.perturbation_groups;
    
    hands = lfp_tfa_cfg.compare.reach_hands;
    spaces = lfp_tfa_cfg.compare.reach_spaces; 
    
    % hand-space labels to be excluded
    exclude_handpsace = {};
    if isfield(lfp_tfa_cfg.compare, 'exclude_handspace')
        exclude_handpsace = lfp_tfa_cfg.compare.exclude_handspace;
    end
    % assign hand space labels
    hs_labels = cell(1, length(hands)*length(spaces) - length(exclude_handpsace));
    reach_hands = cell(1, length(hands)*length(spaces) - length(exclude_handpsace));
    reach_spaces = cell(1, length(hands)*length(spaces) - length(exclude_handpsace));
    hs_idx = 0;
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
            % check if this hand space label should be excluded
            hs_label = [hands{h}, spaces{s}];
            if ~isempty(exclude_handpsace) && ...
                    any(strcmp(lfp_tfa_cfg.compare.exclude_handspace, hs_label))
                continue;
            end
            if strcmp(spaces{s},'R') || strcmp(spaces{s},'L')
                if strcmp(spaces{s},lfp_tfa_cfg.ref_hemisphere)
                    space_label = 'IS';
                else
                    space_label = 'CS';
                end
            else
                space_label = [spaces{s}, 'S'];
            end
            hs_idx = hs_idx + 1;
            reach_hands{hs_idx} = hands{h};
            reach_spaces{hs_idx} = spaces{s};
            hs_labels{hs_idx} = [hand_label ' ' space_label];
        end
    end
    
    % create conditions
    cmp_conditions = struct();
    
    i = 0;
    % should clarify if target belongs to condition
    %for target = targets
        %target_label = target{1};
    for type = task_types
        type_label = ['Type_' num2str(type)];
        for eff = effectors
            eff_label = ['Eff_' num2str(eff)];            
            for ch = choices
                if ch == 0
                    ch_label = 'Instr';
                elseif ch == 1
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
                    cmp_conditions(i).choice = ch;
                    cmp_conditions(i).perturbation = perturbations(p);
                    cmp_conditions(i).perturbation_group = perturbation_groups(p);
                    cmp_conditions(i).hs_labels = hs_labels;
                    cmp_conditions(i).reach_hands = reach_hands;
                    cmp_conditions(i).reach_spaces = reach_spaces;
                    cmp_conditions(i).label = condition_label;
                end
            end
        end
    end
    % end
    

end

