function [ cmp_conditions ] = lfp_tfa_compare_conditions( lfp_tfa_cfg, varargin )
%lfp_tfa_compare_conditions  - Create conditions to compare by permuting
%the compare conditions in the settings
%
% USAGE:
%	cmp_conditions = lfp_tfa_compare_conditions( lfp_tfa_cfg )
%   cmp_conditions = lfp_tfa_compare_conditions( lfp_tfa_cfg, {0, [2, 3]})
%   cmp_conditions = lfp_tfa_compare_conditions( lfp_tfa_cfg, {0, 'all'})
%   cmp_conditions = lfp_tfa_compare_conditions( lfp_tfa_cfg, {0, 'allbutfirst'})
%
% INPUTS:
%		lfp_tfa_cfg     - struct containing the required settings
%           Allowed Fields: see lfp_tfa_settings
%               1. compare.type                 - trial types to be compared
%               2. compare.effector             - trial effectors to be compared
%               3. compare.choice               - trial choices to be compared
%               (0 = instructed, 1 = choice trial)
%               4. compare.success              - trial successes to be compared
%               (0 = unsuccessful trial, 1 = successful trial), this is
%               used only by the ECG analysis pipeline currently
%               4. compare.perturbation         - perturbations to be compared
%               (0 = preinjection, 1 = postinjection)
%               5. compare.reach_hands          - hand labels to compare
%               ('R' = right, 'L' = left)
%               6. compare.reach_spaces         - space labels to compare
%               ('R' = right, 'L' = left)
%               7. ref_hemisphere               - reference hemisphere for
%               hand-space labelling ('R' or 'L', typically, the inactivated
%               hemisphere)
%               8. compare.exclude_handspace    - hand-space labels to be
%               excluded from analysis (this is usually done when there
%               exists not enough trials for some hand-space conditions,
%               for example IHCS and CHIS during choice trials)
%
%       varargin - should be a 1x2 cell array containing the blocks to be
%       considered as pre- and post- injection respectively.
% OUTPUTS:
%		cmp_conditions      - structure containing conditions to compare, a
%		permutation of all given comparison conditions
%
% REQUIRES:
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_define_settings
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

task_types = inf;
if isfield(lfp_tfa_cfg.compare, 'types')
    task_types = lfp_tfa_cfg.compare.types;
end
effectors = inf;
if isfield(lfp_tfa_cfg.compare, 'effectors')
    effectors = lfp_tfa_cfg.compare.effectors;
end

choices = inf;
if isfield(lfp_tfa_cfg.compare, 'choice_trials')
    choices = lfp_tfa_cfg.compare.choice_trials;
end
perturbations = inf;
if isfield(lfp_tfa_cfg.compare, 'perturbations')
    perturbations = lfp_tfa_cfg.compare.perturbations;
end
trial_success = inf;
if isfield(lfp_tfa_cfg.compare, 'success')
    trial_success = lfp_tfa_cfg.compare.success;
end

% if different sessions have different perturbation groups
if nargin > 1
    perturbation_groups = varargin{1};
    % default perturbation groups based on perturbation input
elseif perturbations == 0
    perturbation_groups = {0};
else%if sum (perturbations == [0, 1])
    perturbation_groups = {};%{0, 'all'};
end

hands = lfp_tfa_cfg.compare.reach_hands;
spaces = lfp_tfa_cfg.compare.reach_spaces;

% hand-space labels to be excluded
exclude_handspace = {};
if isfield(lfp_tfa_cfg.compare, 'exclude_handspace')
    exclude_handspace = lfp_tfa_cfg.compare.exclude_handspace;
end
% assign hand space labels
hs_labels = cell(1, length(hands)*length(spaces) - length(exclude_handspace));
reach_hands = cell(1, length(hands)*length(spaces) - length(exclude_handspace));
reach_spaces = cell(1, length(hands)*length(spaces) - length(exclude_handspace));
hs_idx = 0;
for h = 1:length(hands)
    %         if strcmp(hands{h},'R') || strcmp(hands{h},'L')
    %             if strcmp(hands{h},lfp_tfa_cfg.ref_hemisphere)
    %                 hand_label = 'IH';
    %             else
    %                 hand_label = 'CH';
    %             end
    %         else
    %             hand_label = [hands{h} 'H'];
    %         end
    for s = 1:length(spaces)
        % check if this hand space label should be excluded
        hs_label = [hands{h}, spaces{s}];
        if ~isempty(exclude_handspace) && ...
                any(strcmp(lfp_tfa_cfg.compare.exclude_handspace, hs_label))
            continue;
        end
        %             if (strcmp(spaces{s},'R') || strcmp(spaces{s},'L'))
        %                 if strcmp(spaces{s},lfp_tfa_cfg.ref_hemisphere)
        %                     space_label = 'IS';
        %                 else
        %                     space_label = 'CS';
        %                 end
        %             else
        %                 space_label = [spaces{s} 'S'];
        %             end
        hs_idx = hs_idx + 1;
        reach_hands{hs_idx} = hands{h};
        reach_spaces{hs_idx} = spaces{s};
        %hs_labels{hs_idx} = [hand_label ' ' space_label];
        hs_labels{hs_idx} = [hands{h} 'H ' spaces{s} 'S'];
    end
end

% create conditions
cmp_conditions = struct();

% get all combinations of given conditions
conditions = {task_types,effectors,choices,trial_success,perturbations};
tmp = conditions;
[tmp{:}] = ndgrid(conditions{:});
combinations = cell2mat(cellfun(@(m)m(:),tmp,'uni',0));

for i = 1:size(combinations, 1)
    cmp_conditions(i).type = combinations(i, 1);
    cmp_conditions(i).effector = combinations(i, 2);
    cmp_conditions(i).choice = combinations(i, 3);
    cmp_conditions(i).success = combinations(i, 4);
    cmp_conditions(i).perturbation = combinations(i, 5);
    condition_label = lfp_tfa_get_condition_label(cmp_conditions(i), 'long');
    % pre-injection
    if ~isempty(perturbation_groups)
        if cmp_conditions(i).perturbation == 0
            if ~isempty(perturbation_groups(1))
                cmp_conditions(i).perturbation_group = ...
                    perturbation_groups(1);
            end
        end
        if cmp_conditions(i).perturbation == 1
            if ~isempty(perturbation_groups(2))
                cmp_conditions(i).perturbation_group = ...
                    perturbation_groups(2);
            end
        end
    end
    cmp_conditions(i).hs_labels = hs_labels;
    cmp_conditions(i).reach_hands = reach_hands;
    cmp_conditions(i).reach_spaces = reach_spaces;
    cmp_conditions(i).label = condition_label;
end

end

