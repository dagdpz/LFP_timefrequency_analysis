function [ diff_sync ] = lfp_tfa_compute_diff_condition_tfsync( sitepair_sync, diff_condition, stat_test, lfp_tfa_cfg )
%lfp_tfa_compute_diff_condition_tfsync - function to compute the difference
%of LFP-LFP phase synchronization averages between different conditions
%
% USAGE:
%   diff_sync = lfp_tfa_compute_diff_condition_tfsync( sitepair_sync,
%	diff_condition)
%	diff_sync = lfp_tfa_compute_diff_condition_tfsync( sitepair_sync,
%	diff_condition, stat_test, lfp_tfa_cfg )
%
% INPUTS:
%       sitepair_sync   - struct containing the LFP-LFP phase sync
%       spectrogram averages for different conditions as average for single
%       sitepair, or average across sitepairs within single or multiple sessions
%       see lfp_tfa_sitepair_averaged_sync,
%       lfp_tfa_avg_sitepairs_sync, lfp_tfa_avg_sessions_sync
%       diff_condition  - the conditions between which difference has to be
%       calculated, see settings/lfp_tfa_settings_example
%       stat_test       - flag which indicate whether to perform a statistical significance
%       test of differences, set to true while computing difference of
%       averages across sitepair averages of multiple sessions
%		lfp_tfa_cfg     - struct containing the required settings (required
%		only if stat_test = true), see settings/lfp_tfa_settings_example
%           Required fields:
%           fd_rate: Desired false discovery rate
%           fdr_method: method to be used for statistical significance test
%
% OUTPUTS:
%		diff_sync       - struct containing the LFP-LFP sync spectrum
%       difference average between different conditions
%
% REQUIRES:
%
% See also lfp_tfa_sitepair_averaged_sync, lfp_tfa_avg_sitepairs_sync,
% lfp_tfa_avg_sessions_sync
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

diff_sync = [];

if nargin < 3
    stat_test = false;
end

if ~isempty([sitepair_sync.cfg_condition])
    conditions = [sitepair_sync.cfg_condition];
else
    return;
end

for i = 1:length(diff_condition)/2
    diff_sync = struct();
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
    for cn = 1:length(sitepair_sync)
        condition_found = false;
        if strcmp(compare.field, 'choice')
            condition_found = sitepair_sync(cn).cfg_condition.choice == compare.values{1};
            
        elseif strcmp(compare.field, 'perturbation')
            condition_found = sitepair_sync(cn).cfg_condition.perturbation == compare.values{1};
            
        elseif strcmp(compare.field, 'type_eff')
            condition_found = sitepair_sync(cn).cfg_condition.type == compare.values{1}(1) ...
                & sitepair_sync(cn).cfg_condition.effector == compare.values{1}(2);
        elseif strcmp(compare.field, 'reach_hands')
            condition_found = strcmp(unique(sitepair_sync(cn).cfg_condition.reach_hands), compare.values);
            
        elseif strcmp(compare.field, 'reach_spaces')
            condition_found = strcmp(unique(sitepair_sync(cn).cfg_condition.reach_spaces), compare.values);
        end
        % initially load the pre-injection data structure
        if condition_found
            traversed_idx = [traversed_idx cn];
        else
            continue;
        end
        if strcmp(compare.field, 'choice') || strcmp(compare.field, 'perturbation') || strcmp(compare.field, 'type_eff') % the logic is different for hand or space comparisons
            for d = 1:length(sitepair_sync)
                if any(traversed_idx == d), continue; end
                comparison_pair_found = false;
                
                if strcmp(compare.field, 'choice')
                    comparison_pair_found = sitepair_sync(d).cfg_condition.type == sitepair_sync(cn).cfg_condition.type ...
                        & sitepair_sync(d).cfg_condition.effector == sitepair_sync(cn).cfg_condition.effector ...
                        & sitepair_sync(d).cfg_condition.choice == compare.values{2} ...
                        & sitepair_sync(d).cfg_condition.perturbation == sitepair_sync(cn).cfg_condition.perturbation;
                    
                elseif strcmp(compare.field, 'perturbation')
                    comparison_pair_found = sitepair_sync(d).cfg_condition.type == sitepair_sync(cn).cfg_condition.type ...
                        & sitepair_sync(d).cfg_condition.effector == sitepair_sync(cn).cfg_condition.effector ...
                        & sitepair_sync(d).cfg_condition.choice == sitepair_sync(cn).cfg_condition.choice ...
                        & sitepair_sync(d).cfg_condition.perturbation == compare.values{2};
                    
                elseif strcmp(compare.field, 'type_eff')
                    comparison_pair_found = sitepair_sync(d).cfg_condition.type == compare.values{2}(1) ...
                        & sitepair_sync(d).cfg_condition.effector == compare.values{2}(2) ...
                        & sitepair_sync(d).cfg_condition.choice == sitepair_sync(cn).cfg_condition.choice ...
                        & sitepair_sync(d).cfg_condition.perturbation == sitepair_sync(cn).cfg_condition.perturbation;
                end
                if comparison_pair_found
                    
                    % pre injection
                    preinj_sync = sitepair_sync(cn);
                    if ~isfield(preinj_sync.hs_tuned_sync, 'ppc')
                        continue;
                    end
                    % post injection
                    postinj_tfr = sitepair_sync(d);
                    if ~isfield(postinj_tfr.hs_tuned_sync, 'ppc')
                        continue;
                    end
                    
                    %                     if ~isfield(postinj_tfr.hs_tuned_sync.ppc, 'ppcspctrm') || ...
                    %                             ~isfield(preinj_sync.hs_tuned_sync.ppc, 'ppcspctrm')
                    %                         continue;
                    %                     end
                    
                    % if all tests are passed
                    dcn = dcn + 1;
                    
                    diff_sync.difference(dcn) = postinj_tfr;
                    % calculate the number of time frequency bins in all
                    % windows for later bonferroni correction
                    bon_time_sum = 0;
                    for stb = 1:size(preinj_sync.hs_tuned_sync, 1)
                        bon_time_sum = bon_time_sum + size(preinj_sync.hs_tuned_sync(stb, 1).time,2);
                    end
                    bon_factor = bon_time_sum*size(preinj_sync.hs_tuned_sync(1, 1).freq,2);
                    
                    % change the condition label
                    diff_sync.difference(dcn).label = ['( ' postinj_tfr.label, ' - ', ...
                        preinj_sync.label ' )'];
                    
                    diff_sync.difference(dcn).cfg_condition = postinj_tfr.cfg_condition;
                    if strcmp(compare.field, 'choice')
                        diff_sync.difference(dcn).cfg_condition.choice = ['diff' num2str(i)];
                    elseif strcmp(compare.field, 'perturbation')
                        diff_sync.difference(dcn).cfg_condition.perturbation = ['diff' num2str(i)];
                    elseif strcmp(compare.field, 'type_eff')
                        diff_sync.difference(dcn).cfg_condition.type_eff = ['diff' num2str(i)];
                    end
                    
                    % loop through handspace tunings
                    diff_sync.difference(dcn).hs_tuned_sync = postinj_tfr.hs_tuned_sync;
                    for hs = 1:size(postinj_tfr.hs_tuned_sync, 2)
                        for st = 1:size(postinj_tfr.hs_tuned_sync, 1)
                            
                            if isfield(preinj_sync.hs_tuned_sync(st, hs).ppc, 'ppcspctrm') && ...
                                    isfield(postinj_tfr.hs_tuned_sync(st, hs).ppc, 'ppcspctrm') && ...
                                    ~isempty(preinj_sync.hs_tuned_sync(st, hs).ppc.ppcspctrm) && ...
                                    ~isempty(postinj_tfr.hs_tuned_sync(st, hs).ppc.ppcspctrm)
                                ntimebins = min([size(postinj_tfr.hs_tuned_sync(st, hs).ppc.ppcspctrm, 3), ...
                                    size(preinj_sync.hs_tuned_sync(st, hs).ppc.ppcspctrm, 3)]);
                                % calculate the difference
                                diff_sync.difference(dcn).hs_tuned_sync(st, hs).ppc.ppcspctrm = ...
                                    postinj_tfr.hs_tuned_sync(st, hs).ppc.ppcspctrm(:,:,1:ntimebins) - ...
                                    preinj_sync.hs_tuned_sync(st, hs).ppc.ppcspctrm(:,:,1:ntimebins);
                                % statistical test
                                if stat_test
                                    [~, p] = ttest(...
                                        diff_sync.difference(dcn).hs_tuned_sync(st, hs).ppc.ppcspctrm);
                                    if strcmp(lfp_tfa_cfg.correction_method,'FDR');
                                        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p,lfp_tfa_cfg.fd_rate,...
                                            lfp_tfa_cfg.fdr_method,'yes');
                                        diff_sync.difference(dcn).hs_tuned_sync(st, hs).ppc.stat_test.h = h;
                                    elseif strcmp(lfp_tfa_cfg.correction_method,'Bonferroni');
                                        
                                        corrected_p_val_sig = 0.05/bon_factor;
                                        diff_sync.difference(dcn).hs_tuned_sync(st, hs).ppc.stat_test.h = p < corrected_p_val_sig;
                                    end
                                end
                                diff_sync.difference(dcn).hs_tuned_sync(st, hs).ppc.time = ...
                                    postinj_tfr.hs_tuned_sync(st, hs).ppc.time(1:ntimebins);
                                if isfield(postinj_tfr.hs_tuned_sync(st, hs), 'ntrials')
                                    diff_sync.difference(dcn).hs_tuned_sync(st, hs).ntrials = ...
                                        [];
                                end
                                
                            else
                                diff_sync.difference(dcn).hs_tuned_sync(st, hs).ppc.ppcspctrm = [];
                                diff_sync.difference(dcn).hs_tuned_sync(st, hs).ppc.time = [];
                                diff_sync.difference(dcn).hs_tuned_sync(st, hs).ppc.freq = [];
                            end
                        end
                    end
                else
                    continue;
                end
            end
        elseif strcmp(compare.field, 'reach_hands') || strcmp(compare.field, 'reach_spaces') % need to average same hand or same space for each sitepair_sync(d), then do the difference between the 2 averages
            
            
            % calculate the number of time frequency bins in all
            % windows for later bonferroni correction
            bon_time_sum = 0;
            for stb = 1:size(sitepair_sync(1).hs_tuned_sync, 1)
                bon_time_sum = bon_time_sum + size(sitepair_sync(1).hs_tuned_sync(stb, 1).time,2);
            end
            bon_factor = bon_time_sum*size(sitepair_sync(1).hs_tuned_sync(1, 1).freq,2);
            
            for d = 1:length(sitepair_sync)
                diff_sync.difference(d) = sitepair_sync(d);
                if strcmp(compare.field, 'reach_hands')
                    diff_sync.difference(d).hs_tuned_sync = sitepair_sync(d).hs_tuned_sync;
                    diff_sync.difference(d).label = [sitepair_sync(d).label, 'CH - IH'];
                    %                     diff_sync.difference(d).cfg_condition.reach_hands = ['diff' num2str(i)];
                    diff_sync.difference(d).cfg_condition = sitepair_sync(d).cfg_condition;
                    diff_sync.difference(d).cfg_condition.diff = 'hands';
                    
                    
                    
                    for st = 1:size(sitepair_sync(d).hs_tuned_sync, 1) %st is windows here
                        %calculate difference between same space,
                        %opposite hands condition
                        ntimebins_CS = min([size(sitepair_sync(d).hs_tuned_sync(st,1).ppc.ppcspctrm,3)...
                            size(sitepair_sync(d).hs_tuned_sync(st,3).ppc.ppcspctrm,3)]);
                        ntimebins_IS = min([size(sitepair_sync(d).hs_tuned_sync(st,2).ppc.ppcspctrm,3)...
                            size(sitepair_sync(d).hs_tuned_sync(st,4).ppc.ppcspctrm,3)]);
                        
                        if isfield(sitepair_sync(d).hs_tuned_sync, 'ntrials')
                            diff_sync.difference(d).hs_tuned_sync(st, 1).ppc.ppcspctrm= nanmean(sitepair_sync(d).hs_tuned_sync(st,1).ppc.ppcspctrm(:,:,1:ntimebins_CS))...
                                - nanmean(sitepair_sync(d).hs_tuned_sync(st,3).ppc.ppcspctrm(:,:,1:ntimebins_CS));
                            diff_sync.difference(d).hs_tuned_sync(st, 2).ppc.ppcspctrm= nanmean(sitepair_sync(d).hs_tuned_sync(st,2).ppc.ppcspctrm(:,:,1:ntimebins_IS))...
                                - nanmean(sitepair_sync(d).hs_tuned_sync(st,4).ppc.ppcspctrm(:,:,1:ntimebins_IS));
                            
                            
                        else
                            
                            
                            
                            diff_sync.difference(d).hs_tuned_sync(st, 1).ppc.ppcspctrm= sitepair_sync(d).hs_tuned_sync(st,1).ppc.ppcspctrm(:,:,1:ntimebins_CS)...
                                - sitepair_sync(d).hs_tuned_sync(st,3).ppc.ppcspctrm(:,:,1:ntimebins_CS);
                            diff_sync.difference(d).hs_tuned_sync(st, 2).ppc.ppcspctrm= sitepair_sync(d).hs_tuned_sync(st,2).ppc.ppcspctrm(:,:,1:ntimebins_IS)...
                                - sitepair_sync(d).hs_tuned_sync(st,4).ppc.ppcspctrm(:,:,1:ntimebins_IS);
                            
                        end
                        %change hand space label
                        diff_sync.difference(d).hs_tuned_sync(st, 1).hs_label = {'CHCS - IHCS'};
                        diff_sync.difference(d).hs_tuned_sync(st, 2).hs_label = {'CHIS - IHIS'};
                        %empty other hand space condition since not needed
                        %here
                        diff_sync.difference(d).hs_tuned_sync(st,3).ppc.ppcspctrm= [];
                        diff_sync.difference(d).hs_tuned_sync(st,3).ppc.time = [];
                        diff_sync.difference(d).hs_tuned_sync(st,3).ppc.freq = [];
                        diff_sync.difference(d).hs_tuned_sync(st,4).ppc.ppcspctrm= [];
                        diff_sync.difference(d).hs_tuned_sync(st,4).ppc.time = [];
                        diff_sync.difference(d).hs_tuned_sync(st,4).ppc.freq = [];
                        
                        
                        if stat_test == true
                            for hs = 1:2
                                % paired ttest
                                [~, p] = ttest(...
                                    diff_sync.difference(d).hs_tuned_sync(st, hs).ppc.ppcspctrm);
                                % multiple comparison correction
                                if strcmp(lfp_tfa_cfg.correction_method,'FDR');
                                    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p, fd_rate,...
                                        fdr_method, 'yes');
                                    diff_sync.difference(d).hs_tuned_sync(st, hs).ppc.stat_test.h = h;
                                    diff_sync.difference(d).hs_tuned_sync(st, hs).ppc.stat_test.crit_p = crit_p;
                                    diff_sync.difference(d).hs_tuned_sync(st, hs).ppc.stat_test.adj_ci_cvrg = adj_ci_cvrg;
                                    diff_sync.difference(d).hs_tuned_sync(st, hs).ppc.stat_test.adj_p = adj_p;
                                elseif strcmp(lfp_tfa_cfg.correction_method,'Bonferroni');
                                    corrected_p_val_sig = 0.05/bon_factor;
                                    diff_sync.difference(d).hs_tuned_sync(st, hs).ppc.stat_test.h = p < corrected_p_val_sig;
                                end
                            end
                        end
                    end
                    
                elseif strcmp(compare.field, 'reach_spaces')
                    
                    diff_sync.difference(d).hs_tuned_sync = sitepair_sync(d).hs_tuned_sync;
                    diff_sync.difference(d).label = [sitepair_sync(d).label, 'CS - IS'];
                    %                     diff_sync.difference(d).cfg_condition.reach_spaces = ['diff' num2str(i)];
                    diff_sync.difference(d).cfg_condition = sitepair_sync(d).cfg_condition;
                    diff_sync.difference(d).cfg_condition.diff = 'spaces';
                    
                    for st = 1:size(sitepair_sync(d).hs_tuned_sync, 1) %st is windows here
                        %calculate difference between same hands,
                        %opposite space condition
                        if ismember(sitepair_sync(d).cfg_condition.effector,[1 2 3 4 6])
                            
                            ntimebins_CH = min([size(sitepair_sync(d).hs_tuned_sync(st,1).ppc.ppcspctrm,3)...
                                size(sitepair_sync(d).hs_tuned_sync(st,2).ppc.ppcspctrm,3)]);
                            ntimebins_IH = min([size(sitepair_sync(d).hs_tuned_sync(st,3).ppc.ppcspctrm,3)...
                                size(sitepair_sync(d).hs_tuned_sync(st,4).ppc.ppcspctrm,3)]);
                            
                            
                            if isfield(sitepair_sync(d).hs_tuned_sync, 'ntrials')
                                diff_sync.difference(d).hs_tuned_sync(st, 1).ppc.ppcspctrm= nanmean(sitepair_sync(d).hs_tuned_sync(st,1).ppc.ppcspctrm(:,:,1:ntimebins_CH))...
                                    - nanmean(sitepair_sync(d).hs_tuned_sync(st,2).ppc.ppcspctrm(:,:,1:ntimebins_CH));
                                diff_sync.difference(d).hs_tuned_sync(st, 2).ppc.ppcspctrm= nanmean(sitepair_sync(d).hs_tuned_sync(st,3).ppc.ppcspctrm(:,:,1: ntimebins_IH))...
                                    - nanmean(sitepair_sync(d).hs_tuned_sync(st,4).ppc.ppcspctrm(:,:,1: ntimebins_IH));
                                
                            else
                                diff_sync.difference(d).hs_tuned_sync(st, 1).ppc.ppcspctrm= sitepair_sync(d).hs_tuned_sync(st,1).ppc.ppcspctrm(:,:,1:ntimebins_CH)...
                                    - sitepair_sync(d).hs_tuned_sync(st,2).ppc.ppcspctrm(:,:,1:ntimebins_CH);
                                diff_sync.difference(d).hs_tuned_sync(st, 2).ppc.ppcspctrm= sitepair_sync(d).hs_tuned_sync(st,3).ppc.ppcspctrm(:,:,1: ntimebins_IH)...
                                    - sitepair_sync(d).hs_tuned_sync(st,4).ppc.ppcspctrm(:,:,1: ntimebins_IH);
                            end
                            %change hand space label
                            diff_sync.difference(d).hs_tuned_sync(st, 1).hs_label = {'CHCS - CHIS'};
                            diff_sync.difference(d).hs_tuned_sync(st, 2).hs_label = {'IHCS - IHIS'};
                            %empty other hand space condition since not needed
                            %here
                            diff_sync.difference(d).hs_tuned_sync(st,3).ppc.ppcspctrm= [];
                            diff_sync.difference(d).hs_tuned_sync(st,3).ppc.time = [];
                            diff_sync.difference(d).hs_tuned_sync(st,3).ppc.freq = [];
                            diff_sync.difference(d).hs_tuned_sync(st,4).ppc.ppcspctrm= [];
                            diff_sync.difference(d).hs_tuned_sync(st,4).ppc.time = [];
                            diff_sync.difference(d).hs_tuned_sync(st,4).ppc.freq = [];
                            
                            
                            if stat_test == true
                                for hs = 1:2
                                    % paired ttest
                                    [~, p] = ttest(...
                                        diff_sync.difference(d).hs_tuned_sync(st, hs).ppc.ppcspctrm);
                                    % multiple comparison correction
                                    if strcmp(lfp_tfa_cfg.correction_method,'FDR');
                                        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p, fd_rate,...
                                            fdr_method, 'yes');
                                        diff_sync.difference(d).hs_tuned_sync(st, hs).ppc.stat_test.h = h;
                                        diff_sync.difference(d).hs_tuned_sync(st, hs).ppc.stat_test.crit_p = crit_p;
                                        diff_sync.difference(d).hs_tuned_sync(st, hs).ppc.stat_test.adj_ci_cvrg = adj_ci_cvrg;
                                        diff_sync.difference(d).hs_tuned_sync(st, hs).ppc.stat_test.adj_p = adj_p;
                                    elseif strcmp(lfp_tfa_cfg.correction_method,'Bonferroni');
                                        corrected_p_val_sig = 0.05/bon_factor;
                                        diff_sync.difference(d).hs_tuned_sync(st, hs).ppc.stat_test.h = p < corrected_p_val_sig;
                                    end
                                    
                                end
                            end
                        elseif sitepair_sync(d).cfg_condition.effector == 0
                            ntimebins_space = min([size(sitepair_sync(d).hs_tuned_sync(st,1).ppc.ppcspctrm,3)...
                                size(sitepair_sync(d).hs_tuned_sync(st,2).ppc.ppcspctrm,3)]);
                            
                            if isfield(sitepair_sync(d).hs_tuned_sync, 'ntrials')
                                diff_sync.difference(d).hs_tuned_sync(st, 1).ppc.ppcspctrm= nanmean(sitepair_sync(d).hs_tuned_sync(st,1).ppc.ppcspctrm(:,:,1:ntimebins_space),1)...
                                    - nanmean(sitepair_sync(d).hs_tuned_sync(st,2).ppc.ppcspctrm(:,:,1:ntimebins_space),1);
                                
                                
                            else
                                diff_sync.difference(d).hs_tuned_sync(st, 1).ppc.ppcspctrm= sitepair_sync(d).hs_tuned_sync(st,1).ppc.ppcspctrm(:,:,1:ntimebins_space)...
                                    - sitepair_sync(d).hs_tuned_sync(st,2).ppc.ppcspctrm(:,:,1:ntimebins_space);
                                
                            end
                            %change hand space label
                            diff_sync.difference(d).hs_tuned_sync(st, 1).hs_label = {'anyH_CS - anyH_IS'};
                            %empty other hand space condition since not needed
                            %here
                            diff_sync.difference(d).hs_tuned_sync(st,2).ppc.ppcspctrm= [];
                            diff_sync.difference(d).hs_tuned_sync(st,2).ppc.time = [];
                            diff_sync.difference(d).hs_tuned_sync(st,2).ppc.freq = [];
                            diff_sync.difference(d).hs_tuned_sync(st,2).ppc.ppcspctrm= [];
                            diff_sync.difference(d).hs_tuned_sync(st,2).ppc.time = [];
                            diff_sync.difference(d).hs_tuned_sync(st,2).ppc.freq = [];
                            
                            if stat_test == true
                                
                                % paired ttest
                                [~, p] = ttest(...
                                    diff_sync.difference(d).hs_tuned_sync(st, 1).ppc.ppcspctrm);
                                % multiple comparison correction
                                if strcmp(lfp_tfa_cfg.correction_method,'FDR');
                                    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p, fd_rate,...
                                        fdr_method, 'yes');
                                    diff_sync.difference(d).hs_tuned_sync(st, 1).ppc.stat_test.h = h;
                                    diff_sync.difference(d).hs_tuned_sync(st, 1).ppc.stat_test.crit_p = crit_p;
                                    diff_sync.difference(d).hs_tuned_sync(st, 1).ppc.stat_test.adj_ci_cvrg = adj_ci_cvrg;
                                    diff_sync.difference(d).hs_tuned_sync(st, 1).ppc.stat_test.adj_p = adj_p;
                                elseif strcmp(lfp_tfa_cfg.correction_method,'Bonferroni');
                                    corrected_p_val_sig = 0.05/bon_factor;
                                    diff_sync.difference(d).hs_tuned_sync(st, 1).ppc.stat_test.h = p < corrected_p_val_sig;
                                end
                                
                                
                            end
                            
                            
                        end
                    end
                    
                end
                
                
            end
            
            
            
        end
    end
    if isfield(diff_sync, 'difference')
        sitepair_sync = diff_sync.difference;
    end
end
% generate return variable
if isfield(diff_sync, 'difference')
    diff_sync = diff_sync.difference;
else
    diff_sync = [];
end
end

