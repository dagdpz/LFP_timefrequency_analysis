function [ diff_evoked ] = lfp_tfa_compute_difference_condition_ev( lfp_evoked, diff_condition, stat_test, lfp_tfa_cfg )
%lfp_tfa_compute_difference_condition_evoked - function to compute the difference
%of LFP evoked amplitude averages between different conditions
%
% USAGE:
%   diff_evoked = lfp_tfa_compute_difference_condition_ev( lfp_evoked,
%	diff_condition)
%	diff_evoked = lfp_tfa_compute_difference_condition_ev( lfp_evoked,
%	diff_condition, stat_test, lfp_tfa_cfg )
%
% INPUTS:
%       lfp_evoked      - struct containing the LFP time lfp amplitude
%       averages for different conditions as average across
%       multiple trials within a site, or average across sites within
%       single or multiple sessions, see lfp_tfa_plot_site_average_evoked,
%       lfp_tfa_avg_evoked_across_sites, lfp_tfa_avg_evoked_across_sessions
%       diff_condition  - the conditions between which difference has to be
%       calculated, see settings/lfp_tfa_settings_example
%       stat_test       - flag which indicate whether to perform a
%       statistical significance test of differences, set to true while
%       computing difference of averages across site averages of multiple
%       sessions
%		lfp_tfa_cfg     - struct containing the required settings (required
%		only if stat_test = true), see settings/lfp_tfa_settings_example
%           Required fields:
%           fd_rate: Desired false discovery rate
%           fdr_method: method to be used for statistical significance test
%
% OUTPUTS:
%		diff_evoked        - struct containing the LFP time lfp amplitude
%       difference average between different conditions
%
% REQUIRES:
%
% See also lfp_tfa_plot_site_average_evoked,
% lfp_tfa_avg_evoked_across_sites, lfp_tfa_avg_evoked_across_sessions
%
% Author(s):	S.Nair, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2019-02-15:	Created function (Sarath Nair)
% 2019-03-05:	First Revision
% ...
% $Revision: 1.0 $  $Date: 2019-03-05 17:18:00 $
% 2022-08-29:   changed from tfr function to evoked LFP

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%

diff_evoked = [];

% defaults
fd_rate = 0.05;
fdr_method = 'pdep';
if nargin < 3
    stat_test = false;
elseif nargin > 3
    stat_test = lfp_tfa_cfg.plot_significant;
    if isfield(lfp_tfa_cfg, 'fd_rate')
        fd_rate = lfp_tfa_cfg.fd_rate;
    end
    if isfield(lfp_tfa_cfg, 'fd_rate')
        fdr_method = lfp_tfa_cfg.fdr_method;
    end
end

if ~isempty([lfp_evoked.cfg_condition])
    conditions = [lfp_evoked.cfg_condition];
    % elseif ~isempty(
else
    return;
end

for i = 1:length(diff_condition)/2
    diff_evoked = struct();
    compare.field = diff_condition{1, 2*(i-1) + 1};
    compare.values = diff_condition{1, 2*(i-1) + 2};
    % check if conditions to compare exist
    if strcmp(compare.field, 'perturbation')
        if sum([compare.values{:}] == unique([conditions.perturbation])) < 2
            continue;
        end
    elseif strcmp(compare.field, 'choice')
        if sum([compare.values{:}] == unique([conditions.choice])) < 2
            continue;
        end
    elseif strcmp(compare.field, 'type_eff')
        if sum(ismember(vertcat(compare.values{:}), ...
                unique([conditions.type; conditions.effector]', 'rows'), 'rows')) < 2
            continue;
        end
    elseif strcmp(compare.field, 'reach_hands')
        if length(unique([conditions.reach_hands])) < 2
            continue
        end
    elseif strcmp(compare.field, 'reach_spaces')
        if length(unique([conditions.reach_spaces])) < 2
            continue
        end
    end
    
    dcn = 0;
    traversed_idx = [];
    for cn = 1:length(lfp_evoked)
        condition_found = false;
        if strcmp(compare.field, 'choice')
            condition_found = lfp_evoked(cn).cfg_condition.choice == compare.values{1};
            
        elseif strcmp(compare.field, 'perturbation')
            condition_found = lfp_evoked(cn).cfg_condition.perturbation == compare.values{1};
            
        elseif strcmp(compare.field, 'type_eff')
            condition_found = lfp_evoked(cn).cfg_condition.type == compare.values{1}(1) ...
                & lfp_evoked(cn).cfg_condition.effector == compare.values{1}(2);
            
        elseif strcmp(compare.field, 'reach_hands')
            condition_found = strcmp(unique(lfp_evoked(cn).cfg_condition.reach_hands), compare.values);
            
        elseif strcmp(compare.field, 'reach_spaces')
            condition_found = strcmp(unique(lfp_evoked(cn).cfg_condition.reach_spaces), compare.values);
        end
        % initially load the pre-injection data structure
        if condition_found
            traversed_idx = [traversed_idx cn];
        else
            continue;
        end
        
        if strcmp(compare.field, 'choice') || strcmp(compare.field, 'perturbation') || strcmp(compare.field, 'type_eff') % the logic is different for hand or space comparisons
            for d = 1:length(lfp_evoked)
                if any(traversed_idx == d), continue; end
                comparison_pair_found = false;
                
                if strcmp(compare.field, 'choice')
                    comparison_pair_found = lfp_evoked(d).cfg_condition.type == lfp_evoked(cn).cfg_condition.type ...
                        & lfp_evoked(d).cfg_condition.effector == lfp_evoked(cn).cfg_condition.effector ...
                        & lfp_evoked(d).cfg_condition.choice == compare.values{2} ...
                        & lfp_evoked(d).cfg_condition.perturbation == lfp_evoked(cn).cfg_condition.perturbation;
                    
                elseif strcmp(compare.field, 'perturbation')
                    comparison_pair_found = lfp_evoked(d).cfg_condition.type == lfp_evoked(cn).cfg_condition.type ...
                        & lfp_evoked(d).cfg_condition.effector == lfp_evoked(cn).cfg_condition.effector ...
                        & lfp_evoked(d).cfg_condition.choice == lfp_evoked(cn).cfg_condition.choice ...
                        & lfp_evoked(d).cfg_condition.perturbation == compare.values{2};
                    
                elseif strcmp(compare.field, 'type_eff')
                    comparison_pair_found = lfp_evoked(d).cfg_condition.type == compare.values{2}(1) ...
                        & lfp_evoked(d).cfg_condition.effector == compare.values{2}(2) ...
                        & lfp_evoked(d).cfg_condition.choice == lfp_evoked(cn).cfg_condition.choice ...
                        & lfp_evoked(d).cfg_condition.perturbation == lfp_evoked(cn).cfg_condition.perturbation;
                end
                if comparison_pair_found
                    
                    dcn = dcn + 1;
                    %diff_tfr.difference(dcn) = struct();
                    % pre injection
                    preinj_sync = lfp_evoked(cn);
                    if isempty(fieldnames(preinj_sync.avg_across_sessions) )
                        continue;
                    end
                    % post injection
                    postinj_evoked = lfp_evoked(d);
                    if isempty(fieldnames(postinj_evoked.avg_across_sessions))
                        continue;
                    end
                    
                    if ~isfield(postinj_evoked.avg_across_sessions, 'lfp') || ...
                            ~isfield(preinj_sync.avg_across_sessions, 'lfp')
                        continue;
                    end
                    
                    diff_evoked.difference(dcn) = postinj_evoked; %HERE IT IS ASSIGNED AVG ACROSS SESS
                    
                    % change the condition label
                    diff_evoked.difference(dcn).label = ['( ' postinj_evoked.label, ' - ', ...
                        preinj_sync.label ' )'];
                    
                    diff_evoked.difference(dcn).cfg_condition = postinj_evoked.cfg_condition;
                    if strcmp(compare.field, 'choice')
                        %                         diff_tfr.difference(dcn).cfg_condition.choice = ['diff' num2str(i)];
                        diff_evoked.difference(dcn).cfg_condition.diff = 'choice';
                    elseif strcmp(compare.field, 'perturbation')
                        %                         diff_tfr.difference(dcn).cfg_condition.perturbation = ['diff' num2str(i)];
                        diff_evoked.difference(dcn).cfg_condition.diff = 'perturbation';
                    elseif strcmp(compare.field, 'type_eff')
                        %                         diff_tfr.difference(dcn).cfg_condition.type_eff = ['diff' num2str(i)];
                        diff_evoked.difference(dcn).cfg_condition.diff = 'type_eff';
                    end
                    
                    % loop through handspace tunings
                    diff_evoked.difference(dcn).hs_tuned_evoked = postinj_evoked.avg_across_sessions;
                    for hs = 1:size(postinj_evoked.avg_across_sessions, 2)
                        for st = 1:size(postinj_evoked.avg_across_sessions, 1)
                            
                            if isfield(preinj_sync.avg_across_sessions(st, hs).lfp, 'powspctrm') && ...
                                    isfield(postinj_evoked.avg_across_sessions(st, hs).lfp, 'powspctrm') && ...
                                    ~isempty(preinj_sync.avg_across_sessions(st, hs).lfp) && ...
                                    ~isempty(postinj_evoked.avg_across_sessions(st, hs).lfp)
                                ntimebins = min([size(postinj_evoked.avg_across_sessions(st, hs).lfp, 3), ...
                                    size(preinj_sync.avg_across_sessions(st, hs).lfp, 3)]);
                                diff_evoked.difference(dcn).hs_tuned_evoked(st, hs).lfp.time = ...
                                    postinj_evoked.avg_across_sessions(st, hs).lfp.time(1:ntimebins);
                            else
                                ntimebins = length(postinj_evoked.avg_across_sessions(st, hs).time);
                                diff_evoked.difference(d).hs_tuned_evoked = lfp_evoked(d).avg_across_sessions;
                                diff_evoked.difference(d).label = [lfp_evoked(d).label];
                                %                     diff_tfr.difference(d).cfg_condition.reach_spaces = ['diff' num2str(i)];
                                diff_evoked.difference(d).cfg_condition = lfp_evoked(d).cfg_condition;
                                diff_evoked.difference(d).cfg_condition.diff = 'choice';
                                % calculate the difference
                                if isfield(postinj_evoked.avg_across_sessions(st, hs), 'ntrials')
                                    diff_evoked.difference(dcn).hs_tuned_evoked(st, hs).lfp = ...
                                        nanmean(postinj_evoked.avg_across_sessions(st, hs).lfp(:,1:ntimebins), 1) - ...
                                        nanmean(preinj_sync.avg_across_sessions(st, hs).lfp(:,1:ntimebins), 1);
                                    diff_evoked.difference(dcn).hs_tuned_evoked(st, hs).ntrials = ...
                                        [];
                                else
                                    diff_evoked.difference(dcn).hs_tuned_evoked(st, hs).lfp = ...
                                        postinj_evoked.avg_across_sessions(st, hs).lfp(:,1:ntimebins) - ...
                                        preinj_sync.avg_across_sessions(st, hs).lfp(:,1:ntimebins);
                                    % statistical significance test
                                    if stat_test == true
                                        % paired ttest
                                        [~, p] = ttest(...
                                            diff_evoked.difference(dcn).hs_tuned_evoked(st, hs).lfp);
                                        % multiple comparison correction
                                        if strcmp(lfp_tfa_cfg.correction_method,'FDR');
                                            [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p, fd_rate,...
                                                fdr_method, 'yes');
                                            diff_evoked.difference(dcn).hs_tuned_evoked(st, hs).lfp.stat_test.h = h;
                                            diff_evoked.difference(dcn).hs_tuned_evoked(st, hs).lfp.stat_test.crit_p = crit_p;
                                            diff_evoked.difference(dcn).hs_tuned_evoked(st, hs).lfp.stat_test.adj_ci_cvrg = adj_ci_cvrg;
                                            diff_evoked.difference(dcn).hs_tuned_evoked(st, hs).lfp.stat_test.adj_p = adj_p;
                                        elseif strcmp(lfp_tfa_cfg.correction_method,'Bonferroni');
                                            corrected_p_val_sig = 0.05/(size(preinj_sync.avg_across_sessions(st, hs).lfp.time,2)*size(preinj_sync.avg_across_sessions(st, hs).lfp.lfp,2));
                                            diff_evoked.difference(dcn).hs_tuned_evoked(st, hs).lfp.stat_test.h = p < corrected_p_val_sig;
                                        end
                                    end
                                end
                                
                                %                             else
                                %                                 diff_evoked.difference(dcn).hs_tuned_evoked(st, hs).lfp = [];
                                %                                 diff_evoked.difference(dcn).hs_tuned_evoked(st, hs).lfp.time = [];
                                %                                 diff_evoked.difference(dcn).hs_tuned_evoked(st, hs).lfp.lfp = [];
                            end
                            
                        end
                    end
                else
                    continue;
                end
            end
            
            
        elseif strcmp(compare.field, 'reach_hands') || strcmp(compare.field, 'reach_spaces') % need to average same hand or same space for each lfp_tfr(d), then do the difference between the 2 averages
            
            
            for d = 1:length(lfp_evoked)
                %                 diff_evoked.difference(d) = lfp_evoked(d);
                if strcmp(compare.field, 'reach_hands')
                    diff_evoked.difference(d).hs_tuned_evoked = lfp_evoked(d).avg_across_sessions;
                    diff_evoked.difference(d).label = [lfp_evoked(d).label, 'CH - IH'];
                    %                     diff_tfr.difference(d).cfg_condition.reach_hands = ['diff' num2str(i)];
                    diff_evoked.difference(d).cfg_condition = lfp_evoked(d).cfg_condition;
                    diff_evoked.difference(d).cfg_condition.diff = 'hands';
                    
                    for st = 1:size(lfp_evoked(d).avg_across_sessions, 1) %st is windows here
                        %calculate difference between same space,
                        %opposite hands condition
                        ntimebins_CS = min([size(lfp_evoked(d).avg_across_sessions(st,1).lfp, 2)...
                            size(lfp_evoked(d).avg_across_sessions(st,3).lfp, 2)]);
                        ntimebins_IS = min([size(lfp_evoked(d).avg_across_sessions(st,2).lfp, 2)...
                            size(lfp_evoked(d).avg_across_sessions(st,4).lfp, 2)]);
                        
                        if isfield(lfp_evoked(d).avg_across_sessions, 'ntrials')
                            diff_evoked.difference(d).hs_tuned_evoked(st, 1).lfp = nanmean(lfp_evoked(d).avg_across_sessions(st,1).lfp(:,1:ntimebins_CS))...
                                - nanmean(lfp_evoked(d).avg_across_sessions(st,3).lfp(:,1:ntimebins_CS));
                            diff_evoked.difference(d).hs_tuned_evoked(st, 2).lfp = nanmean(lfp_evoked(d).avg_across_sessions(st,2).lfp(:,1:ntimebins_IS))...
                                - nanmean(lfp_evoked(d).avg_across_sessions(st,4).lfp(:,1:ntimebins_IS));
                        else
                            diff_evoked.difference(d).hs_tuned_evoked(st, 1).lfp = lfp_evoked(d).avg_across_sessions(st,1).lfp(:,1:ntimebins_CS)...
                                - lfp_evoked(d).avg_across_sessions(st,3).lfp(:,1:ntimebins_CS);
                            diff_evoked.difference(d).hs_tuned_evoked(st, 2).lfp = lfp_evoked(d).avg_across_sessions(st,2).lfp(:,1:ntimebins_IS)...
                                - lfp_evoked(d).avg_across_sessions(st,4).lfp(:,1:ntimebins_IS);
                            
                        end
                        %change hand space label
                        diff_evoked.difference(d).hs_tuned_evoked(st, 1).hs_label = {'CHCS - IHCS'};
                        diff_evoked.difference(d).hs_tuned_evoked(st, 2).hs_label = {'CHIS - IHIS'};
                        %empty other hand space condition since not needed
                        %here
                        diff_evoked.difference(d).hs_tuned_evoked(st,3).lfp = [];
                        diff_evoked.difference(d).hs_tuned_evoked(st,3).lfp.time = [];
                        diff_evoked.difference(d).hs_tuned_evoked(st,3).lfp.lfp = [];
                        diff_evoked.difference(d).hs_tuned_evoked(st,4).lfp = [];
                        diff_evoked.difference(d).hs_tuned_evoked(st,4).lfp.time = [];
                        diff_evoked.difference(d).hs_tuned_evoked(st,4).lfp.lfp = [];
                        
                        
                        if stat_test == true
                            for hs = 1:2
                                % paired ttest
                                [~, p] = ttest(...
                                    diff_evoked.difference(d).hs_tuned_evoked(st, hs).lfp);
                                % multiple comparison correction
                                if strcmp(lfp_tfa_cfg.correction_method,'FDR');
                                    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p, fd_rate,...
                                        fdr_method, 'yes');
                                    diff_evoked.difference(d).hs_tuned_evoked(st, hs).lfp.stat_test.h = h;
                                    diff_evoked.difference(d).hs_tuned_evoked(st, hs).lfp.stat_test.crit_p = crit_p;
                                    diff_evoked.difference(d).hs_tuned_evoked(st, hs).lfp.stat_test.adj_ci_cvrg = adj_ci_cvrg;
                                    diff_evoked.difference(d).hs_tuned_evoked(st, hs).lfp.stat_test.adj_p = adj_p;
                                elseif strcmp(lfp_tfa_cfg.correction_method,'Bonferroni');
                                    corrected_p_val_sig = 0.05/size(lfp_evoked(d).avg_across_sessions(st,hs).lfp,1);
                                    diff_evoked.difference(d).hs_tuned_evoked(st, hs).lfp.stat_test.h = p < corrected_p_val_sig;
                                end
                            end
                        end
                    end
                    
                elseif strcmp(compare.field, 'reach_spaces')
                    
                    diff_evoked.difference(d).hs_tuned_evoked = lfp_evoked(d).avg_across_sessions;
                    diff_evoked.difference(d).label = [lfp_evoked(d).label, 'CS - IS'];
                    %                     diff_tfr.difference(d).cfg_condition.reach_spaces = ['diff' num2str(i)];
                    diff_evoked.difference(d).cfg_condition = lfp_evoked(d).cfg_condition;
                    diff_evoked.difference(d).cfg_condition.diff = 'spaces';
                    
                    for st = 1:size(lfp_evoked(d).avg_across_sessions, 1) %st is windows here
                        %calculate difference between same hands,
                        %opposite space condition
                        if ismember(lfp_evoked(d).cfg_condition.effector,[1 2 3 4 6])
                            
                            ntimebins_CH = min([size(lfp_evoked(d).avg_across_sessions(st,1).lfp, 3)...
                                size(lfp_evoked(d).avg_across_sessions(st,2).lfp, 3)]);
                            ntimebins_IH = min([size(lfp_evoked(d).avg_across_sessions(st,3).lfp, 3)...
                                size(lfp_evoked(d).avg_across_sessions(st,4).lfp, 3)]);
                            
                            
                            if isfield(lfp_evoked(d).avg_across_sessions, 'ntrials')
                                diff_evoked.difference(d).hs_tuned_evoked(st, 1).lfp = nanmean(lfp_evoked(d).avg_across_sessions(st,1).lfp(:,:,1:ntimebins_CH))...
                                    - nanmean(lfp_evoked(d).avg_across_sessions(st,2).lfp(:,:,1:ntimebins_CH));
                                diff_evoked.difference(d).hs_tuned_evoked(st, 2).lfp = nanmean(lfp_evoked(d).avg_across_sessions(st,3).lfp(:,:,1: ntimebins_IH))...
                                    - nanmean(lfp_evoked(d).avg_across_sessions(st,4).lfp(:,:,1: ntimebins_IH));
                                
                            else
                                diff_evoked.difference(d).hs_tuned_evoked(st, 1).lfp = lfp_evoked(d).avg_across_sessions(st,1).lfp(:,:,1:ntimebins_CH)...
                                    - lfp_evoked(d).avg_across_sessions(st,2).lfp(:,:,1:ntimebins_CH);
                                diff_evoked.difference(d).hs_tuned_evoked(st, 2).lfp = lfp_evoked(d).avg_across_sessions(st,3).lfp(:,:,1: ntimebins_IH)...
                                    - lfp_evoked(d).avg_across_sessions(st,4).lfp(:,:,1: ntimebins_IH);
                            end
                            %change hand space label
                            diff_evoked.difference(d).hs_tuned_evoked(st, 1).hs_label = {'CHCS - CHIS'};
                            diff_evoked.difference(d).hs_tuned_evoked(st, 2).hs_label = {'IHCS - IHIS'};
                            %empty other hand space condition since not needed
                            %here
                            diff_evoked.difference(d).hs_tuned_evoked(st,3).lfp = [];
                            diff_evoked.difference(d).hs_tuned_evoked(st,3).lfp.time = [];
                            diff_evoked.difference(d).hs_tuned_evoked(st,3).lfp.lfp = [];
                            diff_evoked.difference(d).hs_tuned_evoked(st,4).lfp = [];
                            diff_evoked.difference(d).hs_tuned_evoked(st,4).lfp.time = [];
                            diff_evoked.difference(d).hs_tuned_evoked(st,4).lfp.lfp = [];
                            
                            
                            if stat_test == true
                                for hs = 1:2
                                    % paired ttest
                                    [~, p] = ttest(...
                                        diff_evoked.difference(d).hs_tuned_evoked(st, hs).lfp);
                                    % multiple comparison correction
                                    if strcmp(lfp_tfa_cfg.correction_method,'FDR');
                                        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p, fd_rate,...
                                            fdr_method, 'yes');
                                        diff_evoked.difference(d).hs_tuned_evoked(st, hs).lfp.stat_test.h = h;
                                        diff_evoked.difference(d).hs_tuned_evoked(st, hs).lfp.stat_test.crit_p = crit_p;
                                        diff_evoked.difference(d).hs_tuned_evoked(st, hs).lfp.stat_test.adj_ci_cvrg = adj_ci_cvrg;
                                        diff_evoked.difference(d).hs_tuned_evoked(st, hs).lfp.stat_test.adj_p = adj_p;
                                    elseif strcmp(lfp_tfa_cfg.correction_method,'Bonferroni');
                                        corrected_p_val_sig = 0.05/size(lfp_evoked(d).avg_across_sessions(st,hs).lfp,1);
                                        diff_evoked.difference(d).hs_tuned_evoked(st, hs).lfp.stat_test.h = p < corrected_p_val_sig;
                                    end
                                    
                                end
                            end
                        elseif lfp_evoked(d).cfg_condition.effector == 0
                            ntimebins_space = min([size(lfp_evoked(d).avg_across_sessions(st,1).lfp, 2)...
                                size(lfp_evoked(d).avg_across_sessions(st,2).lfp, 2)]);
                            
                            if isfield(lfp_evoked(d).avg_across_sessions, 'ntrials')
                                diff_evoked.difference(d).hs_tuned_evoked(st, 1).lfp = nanmean(lfp_evoked(d).avg_across_sessions(st,1).lfp(:,1:ntimebins_space),1)...
                                    - nanmean(lfp_evoked(d).avg_across_sessions(st,2).lfp(:,1:ntimebins_space),1);
                                
                                
                            else
                                diff_evoked.difference(d).hs_tuned_evoked(st, 1).lfp = lfp_evoked(d).avg_across_sessions(st,1).lfp(:,1:ntimebins_space)...
                                    - lfp_evoked(d).avg_across_sessions(st,2).lfp(:,1:ntimebins_space);
                                
                            end
                            %change hand space label
                            diff_evoked.difference(d).hs_tuned_evoked(st, 1).hs_label = {'anyH_CS - anyH_IS'};
                            %empty other hand space condition since not needed
                            %here
                            diff_evoked.difference(d).hs_tuned_evoked(st,2).lfp = [];
                            %                             diff_evoked.difference(d).hs_tuned_evoked(st,2).lfp.time = [];
                            %                             diff_evoked.difference(d).hs_tuned_evoked(st,2).lfp.lfp = [];
                            
                            
                            if stat_test == true
                                
                                % paired ttest
                                [~, p] = ttest(...
                                    diff_evoked.difference(d).hs_tuned_evoked(st, 1).lfp);
                                % multiple comparison correction
                                if strcmp(lfp_tfa_cfg.correction_method,'FDR');
                                    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p, fd_rate,...
                                        fdr_method, 'yes');
                                    diff_evoked.difference(d).hs_tuned_evoked(st, 1).lfp.stat_test.h = h;
                                    diff_evoked.difference(d).hs_tuned_evoked(st, 1).lfp.stat_test.crit_p = crit_p;
                                    diff_evoked.difference(d).hs_tuned_evoked(st, 1).lfp.stat_test.adj_ci_cvrg = adj_ci_cvrg;
                                    diff_evoked.difference(d).hs_tuned_evoked(st, 1).lfp.stat_test.adj_p = adj_p;
                                elseif strcmp(lfp_tfa_cfg.correction_method,'Bonferroni');
                                    corrected_p_val_sig = 0.05/size(lfp_evoked(d).avg_across_sessions(st,1).lfp,1);
                                    diff_evoked.difference(d).hs_tuned_evoked(st, 1).lfp.stat_test.h = p < corrected_p_val_sig;
                                end
                                
                                
                            end
                            
                            
                        end
                    end
                    
                end
                
                
            end
            
            
            
        end
    end
    if isfield(diff_evoked, 'difference')
        lfp_evoked = diff_evoked.difference;
    end
end
% generate return variable

if isfield(diff_evoked, 'difference')
    diff_evoked = diff_evoked.difference;
    if isfield(diff_evoked, 'avg_across_sessions')
        diff_evoked = rmfield(diff_evoked, 'avg_across_sessions');
    end
else
    diff_evoked = [];
end
end

