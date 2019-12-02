function [ diff_tfr ] = lfp_tfa_compute_difference_condition_tfr( lfp_tfr, diff_condition, stat_test, lfp_tfa_cfg )
%lfp_tfa_compute_difference_condition_tfr - function to compute the difference
%of LFP time frequency spectrogram averages between different conditions
%
% USAGE:
%   diff_tfr = lfp_tfa_compute_difference_condition_tfr( lfp_tfr,
%	diff_condition)
%	diff_tfr = lfp_tfa_compute_difference_condition_tfr( lfp_tfr,
%	diff_condition, stat_test, lfp_tfa_cfg )
%
% INPUTS:
%       lfp_tfr   - struct containing the LFP time frequency power
%       spectrogram averages for different conditions as average across 
%       multiple trials within a site, or average across sites within 
%       single or multiple sessions, see lfp_tfa_plot_site_average_tfr,
%       lfp_tfa_avg_tfr_across_sites, lfp_tfa_avg_tfr_across_sessions
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
%		diff_tfr        - struct containing the LFP time freq spectrogram
%       difference average between different conditions
%
% REQUIRES:	
%
% See also lfp_tfa_plot_site_average_tfr,
% lfp_tfa_avg_tfr_across_sites, lfp_tfa_avg_tfr_across_sessions
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

    diff_tfr = [];
    
    % defaults
    fd_rate = 0.05;
    fdr_method = 'pdep';
    if nargin < 3
        stat_test = false;        
    elseif nargin > 3
        if isfield(lfp_tfa_cfg, 'fd_rate')
            fd_rate = lfp_tfa_cfg.fd_rate;
        end
        if isfield(lfp_tfa_cfg, 'fd_rate')
            fdr_method = lfp_tfa_cfg.fdr_method;
        end        
    end

    if ~isempty([lfp_tfr.cfg_condition])
        conditions = [lfp_tfr.cfg_condition];
    else
        return;
    end
    
    for i = 1:length(diff_condition)/2
        diff_tfr = struct();
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
        end
    
        dcn = 0;
        traversed_idx = [];
        for cn = 1:length(lfp_tfr)
            condition_found = false;
            if strcmp(compare.field, 'choice')
                condition_found = lfp_tfr(cn).cfg_condition.choice == compare.values{1};

            elseif strcmp(compare.field, 'perturbation')
                condition_found = lfp_tfr(cn).cfg_condition.perturbation == compare.values{1};
                
            elseif strcmp(compare.field, 'type_eff')
                condition_found = lfp_tfr(cn).cfg_condition.type == compare.values{1}(1) ...
                    & lfp_tfr(cn).cfg_condition.effector == compare.values{1}(2);

            end
            % initially load the pre-injection data structure
            if condition_found
                traversed_idx = [traversed_idx cn];            
            else
                continue;
            end
            for d = 1:length(lfp_tfr)
                if any(traversed_idx == d), continue; end
                comparison_pair_found = false;

                if strcmp(compare.field, 'choice')
                    comparison_pair_found = lfp_tfr(d).cfg_condition.type == lfp_tfr(cn).cfg_condition.type ...
                        & lfp_tfr(d).cfg_condition.effector == lfp_tfr(cn).cfg_condition.effector ...
                        & lfp_tfr(d).cfg_condition.choice == compare.values{2} ...
                        & lfp_tfr(d).cfg_condition.perturbation == lfp_tfr(cn).cfg_condition.perturbation;

                elseif strcmp(compare.field, 'perturbation')
                    comparison_pair_found = lfp_tfr(d).cfg_condition.type == lfp_tfr(cn).cfg_condition.type ...
                        & lfp_tfr(d).cfg_condition.effector == lfp_tfr(cn).cfg_condition.effector ...
                        & lfp_tfr(d).cfg_condition.choice == lfp_tfr(cn).cfg_condition.choice ...
                        & lfp_tfr(d).cfg_condition.perturbation == compare.values{2};
                    
                elseif strcmp(compare.field, 'type_eff')
                    comparison_pair_found = lfp_tfr(d).cfg_condition.type == compare.values{2}(1) ...
                        & lfp_tfr(d).cfg_condition.effector == compare.values{2}(2) ...
                        & lfp_tfr(d).cfg_condition.choice == lfp_tfr(cn).cfg_condition.choice ...
                        & lfp_tfr(d).cfg_condition.perturbation == lfp_tfr(cn).cfg_condition.perturbation;
                end
                if comparison_pair_found

                    dcn = dcn + 1;
                    %diff_tfr.difference(dcn) = struct();
                    % pre injection
                    preinj_sync = lfp_tfr(cn);
                    if isempty(fieldnames(preinj_sync.hs_tuned_tfs) )
                        continue;
                    end
                    % post injection
                    postinj_tfr = lfp_tfr(d);
                    if isempty(fieldnames(postinj_tfr.hs_tuned_tfs))
                        continue;
                    end

                    if ~isfield(postinj_tfr.hs_tuned_tfs, 'freq') || ...
                            ~isfield(preinj_sync.hs_tuned_tfs, 'freq')
                        continue;
                    end
                    
                    diff_tfr.difference(dcn) = postinj_tfr;
                    
                    % change the condition label
                    diff_tfr.difference(dcn).label = ['( ' postinj_tfr.label, ' - ', ...
                                preinj_sync.label ' )'];  

                    diff_tfr.difference(dcn).cfg_condition = postinj_tfr.cfg_condition;
                    if strcmp(compare.field, 'choice')                        
                        diff_tfr.difference(dcn).cfg_condition.choice = ['diff' num2str(i)];
                    elseif strcmp(compare.field, 'perturbation')
                        diff_tfr.difference(dcn).cfg_condition.perturbation = ['diff' num2str(i)];
                    elseif strcmp(compare.field, 'type_eff')
                        diff_tfr.difference(dcn).cfg_condition.type_eff = ['diff' num2str(i)];
                    end

                    % loop through handspace tunings
                    diff_tfr.difference(dcn).hs_tuned_tfs = postinj_tfr.hs_tuned_tfs;
                    for hs = 1:size(postinj_tfr.hs_tuned_tfs, 2)
                        for st = 1:size(postinj_tfr.hs_tuned_tfs, 1)
                            
                            if isfield(preinj_sync.hs_tuned_tfs(st, hs).freq, 'powspctrm') && ...
                                    isfield(postinj_tfr.hs_tuned_tfs(st, hs).freq, 'powspctrm') && ...
                                    ~isempty(preinj_sync.hs_tuned_tfs(st, hs).freq.powspctrm) && ...
                                ~isempty(postinj_tfr.hs_tuned_tfs(st, hs).freq.powspctrm)
                                ntimebins = min([size(postinj_tfr.hs_tuned_tfs(st, hs).freq.powspctrm, 3), ...
                                    size(preinj_sync.hs_tuned_tfs(st, hs).freq.powspctrm, 3)]);
                                diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.time = ...
                                        postinj_tfr.hs_tuned_tfs(st, hs).freq.time(1:ntimebins);  
                                % calculate the difference
                                if isfield(postinj_tfr.hs_tuned_tfs(st, hs), 'ntrials')
                                    diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
                                        nanmean(postinj_tfr.hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins), 1) - ...
                                        nanmean(preinj_sync.hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins), 1);
                                    diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).ntrials = ...
                                        [];
                                else
                                    diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
                                        postinj_tfr.hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins) - ...
                                        preinj_sync.hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins);
                                    % statistical significance test
                                    if stat_test == true
                                        % paired ttest
                                        [~, p] = ttest(...
                                            diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.powspctrm);
                                        % multiple comparison correction
                                        [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(p, fd_rate,...
                                            fdr_method, 'yes');
                                        diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.stat_test.h = h;
                                        diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.stat_test.crit_p = crit_p;
                                        diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.stat_test.adj_ci_cvrg = adj_ci_cvrg;
                                        diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.stat_test.adj_p = adj_p;
                                    end
                                end
                                
                            else
                                diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.powspctrm = [];
                                diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.time = [];
                                diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).freq.freq = [];
                            end
                        end
                    end
                else
                    continue;
                end
            end
        end
        if isfield(diff_tfr, 'difference')
            lfp_tfr = diff_tfr.difference;
        end
    end
    % generate return variable
    if isfield(diff_tfr, 'difference')
        diff_tfr = diff_tfr.difference;
    else
        diff_tfr = [];
    end
end

