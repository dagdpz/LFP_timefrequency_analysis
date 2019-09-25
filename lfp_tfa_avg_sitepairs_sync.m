function [result_matfile, sitepairs_avg] = lfp_tfa_avg_sitepairs_sync(sessions_info, lfp_tfa_cfg, varargin)
%lfp_tfa_avg_tfr_across_sites  - Condition-based LFP time frequency
%response average across many site averages
%
% USAGE:
%	sites_avg = lfp_tfa_avg_tfr_across_sites(lfp_tfr, lfp_tfa_cfg)
%
% INPUTS:
%		lfp_tfr     	- struct containing the condition-based LFP time freq spectrogram for
%		indiviual sites, output of lfp_tfa_plot_site_average_tfr.m
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields:
%               1. conditions          - trial conditions to compare, see
%               lfp_tfa_settings.m and lfp_tfa_compare_conditions.m
%               2. root_results_fldr   - root folder where results are saved
%               3. compare.targets     - targets to compare, see lfp_tfa_settings.m
%               4. 
% OUTPUTS:
%		sites_avg      - structure containing condition-based evoked LFP
%		response averaged across multiple sites
%
% REQUIRES:	lfp_tfa_plot_hs_tuned_tfr, lfp_tfa_compute_diff_tfr
%
% See also lfp_tfa_settings, lfp_tfa_define_settings, lfp_tfa_compare_conditions, 
% lfp_tfa_plot_site_average_tfr, lfp_tfa_compute_diff_tfr
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

    close all;

    % results folder
    if nargin > 2
        results_fldr = varargin{3};
    elseif length(sessions_info) == 1
        results_fldr = fullfile(sessions_info.lfp_sync_results_fldr, 'session_avg');
    else
        results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sites', 'LFP_LFP_Synch');
    end
    
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end   
    
    
    % Average TFR across sites
    sitepairs_avg = struct();
    
    for t = 1:length(lfp_tfa_cfg.compare.target_pairs)
        fprintf('Target pair: %s - %s\n', lfp_tfa_cfg.compare.target_pairs{t}{:});
        sitepairs_avg(t).targets = lfp_tfa_cfg.compare.target_pairs{t};
        % initialise struct for all conditions
        for cn = 1:length(lfp_tfa_cfg.conditions)
            sitepairs_avg(t).condition(cn).hs_tuned_sync = struct();
            sitepairs_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
            sitepairs_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
            % initialize number of site pairs for each handspace
            % label
            for st = 1:size(lfp_tfa_cfg.analyse_states, 1)
                for hs = 1:size(lfp_tfa_cfg.conditions(1).hs_labels, 2)
                    sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).nsites = 0;
                    sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm = [];
                end
            end                        
        end
        for i = 1:length(sessions_info) 
            sitepairs_fldr = dir(sessions_info(i).lfp_sync_results_fldr);
            sitepairs_fldr(1:2) = [];
            for j = 1:length(sitepairs_fldr)
%                 if ~all(strcmp(sort(sessions_info(i).sitepair_info(j).targets), sort(lfp_tfa_cfg.compare.target_pairs{t})))
%                     continue;
%                 end
%                     if ~sessions_info(i).sitepair_info(j).use_for_avg
%                         continue;
%                     end
                % load the analysis results for this sitepair, variable
                % name is sitepair_sync
                %fprintf('Loading site pair %s - %s\n', sessions_info(i).sitepair_info(j).sites{:});
                if strcmp(sitepairs_fldr(j).name, 'session_avg')
                    continue;
                end
                sitepair_sync_file = dir(fullfile(...
                    sessions_info(i).lfp_sync_results_fldr, sitepairs_fldr(j).name, '*.mat'));
                sitepair_sync_filepath = fullfile(sessions_info(i).lfp_sync_results_fldr, ...
                    sitepairs_fldr(j).name, sitepair_sync_file(1).name);
                load(sitepair_sync_filepath, 'sitepair_sync');
                
                if ~all(strcmp(sort(sitepair_sync.targets), sort(lfp_tfa_cfg.compare.target_pairs{t})))
                    continue;
                end
                
                sitepair_sync.use_for_avg = 1;
                for cn = 1:length(lfp_tfa_cfg.conditions)
                    sitepair_sync.use_for_avg = ...
                        ~any(sitepair_sync.condition(cn).ntrials - sitepair_sync.condition(cn).noisytrials < 5);
                end
                if sitepair_sync.use_for_avg == 0
                    continue;
                end
                % compute condition-based average across sitepairs
                for cn = 1:length(lfp_tfa_cfg.conditions)
                    fprintf('Condition %s\n', lfp_tfa_cfg.conditions(cn).label);     
                                
                    if ~isempty(sitepair_sync.condition(cn).hs_tuned_sync) && ... 
                        isfield(sitepair_sync.condition(cn).hs_tuned_sync, 'ppc') 
                        
                        for st = 1:size(sitepair_sync.condition(cn).hs_tuned_sync, 1)
                            for hs = 1:size(sitepair_sync.condition(cn).hs_tuned_sync, 2)
                                if isfield(sitepair_sync.condition(cn).hs_tuned_sync(st, hs), 'ppc') ...
                                        && ~isempty(sitepair_sync.condition(cn).hs_tuned_sync(st, hs).ppc) ...
                                        && isfield(sitepair_sync.condition(cn).hs_tuned_sync(st, hs).ppc, 'ppcspctrm') ...
                                        && ~isempty(sitepair_sync.condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm)
                                    % increment number of site pairs
                                    sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).nsites = ...
                                        sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).nsites + 1;
                                    if sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).nsites == 1%~isfield(sessions_avg.cond_based_tfs(cn).tfs_across_sessions, 'powspctrm')
                                        sitepairs_avg(t).condition(cn).hs_tuned_sync(st,hs).hs_label ...
                                            = sitepair_sync.condition(cn).hs_tuned_sync(st, hs).hs_label;
                                        sitepairs_avg(t).condition(cn).hs_tuned_sync(st,hs).state ...
                                            = sitepair_sync.condition(cn).hs_tuned_sync(st, hs).state;
                                        sitepairs_avg(t).condition(cn).hs_tuned_sync(st,hs).state_name ...
                                            = sitepair_sync.condition(cn).hs_tuned_sync(st, hs).state_name;
                                        sitepairs_avg(t).condition(cn).label = ...
                                            sitepair_sync.condition(cn).label;
                                        sitepairs_avg(t).condition(cn).cfg_condition = ...
                                            sitepair_sync.condition(cn).cfg_condition;
                                        sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.time ...
                                            = sitepair_sync.condition(cn).hs_tuned_sync(st, hs).ppc.time;
                                        sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.freq ...
                                            = sitepair_sync.condition(cn).hs_tuned_sync(st, hs).ppc.freq;                                        
                                        sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.cfg ...
                                            = sitepair_sync.condition(cn).hs_tuned_sync(st, hs).ppc.cfg;
                                        sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm ...
                                            = cat(1, sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm, sitepair_sync.condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm);
                                        
%                                         if isfield(sitepair_sync.condition(cn).hs_tuned_tfs(st, hs), 'nsites')
%                                             sitepairs_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsites ...
%                                                 = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs).nsites;
%                                         end
                                    else
                                        ntimebins = size(sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm, 3);
                                        % average same number of time bins
                                        if ntimebins > length(sitepair_sync.condition(cn).hs_tuned_sync(st, hs).ppc.time)
                                            ntimebins = length(sitepair_sync.condition(cn).hs_tuned_sync(st, hs).ppc.time);
                                        end
                                        sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm ...
                                            = cat(1, sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm(:,:,1:ntimebins), ...
                                            (sitepair_sync.condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm(:,:,1:ntimebins)) ...
                                            );
%                                         if isfield(sitepairs_avg(t).condition(cn).hs_tuned_tfs(st,hs), 'nsites')
%                                             sitepairs_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsites ...
%                                                 = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs).nsites + ...
%                                                 sitepairs_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsites;
%                                         end                               
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end

        % plot
        for cn = 1:length(sitepairs_avg(t).condition)
%             if ~isempty(sitepairs_avg(t).condition(cn).hs_tuned_sync)
%                 if isfield(sitepairs_avg(t).condition(cn).hs_tuned_sync, 'ppc')
%                     for st = 1:size(sitepairs_avg(t).condition(cn).hs_tuned_sync, 1)
%                         for hs = 1:size(sitepairs_avg(t).condition(cn).hs_tuned_sync, 2)
%                             if isempty(sitepairs_avg(t).condition(cn).hs_tuned_sync(st,hs).ppc)
%                                 continue;
%                             end
%                             sitepairs_avg(t).condition(cn).hs_tuned_sync(st,hs).ppc.ppcspctrm_std = ...
%                                 nanstd( ...
%                                 sitepairs_avg(t).condition(cn).hs_tuned_sync(st,hs).ppc.ppcspctrm, 0, 1);
%                             sitepairs_avg(t).condition(cn).hs_tuned_sync(st,hs).ppc.ppcspctrm = ...
%                                 nanmean( ...
%                                 sitepairs_avg(t).condition(cn).hs_tuned_sync(st,hs).ppc.ppcspctrm, 1);
%                         end
%                     end
%                 end
%             end


            if ~isempty(sitepairs_avg(t).condition(cn).hs_tuned_sync)
                if isfield(sitepairs_avg(t).condition(cn).hs_tuned_sync,... 
                        'ppc')
                    plottitle = sprintf('%s - %s (ref: %s) %s', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                        lfp_tfa_cfg.ref_hemisphere, sitepairs_avg(t).condition(cn).label);
                    result_file = fullfile(results_fldr, ...
                        sprintf('%s-%s_%s.png', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                        sitepairs_avg(t).condition(cn).label));
                    lfp_tfa_plot_hs_tuned_sync(sitepairs_avg(t).condition(cn).hs_tuned_sync, ...
                        lfp_tfa_cfg, plottitle, result_file, 'imscale', [0, 1]);
                end
            end        
        end

        % average difference between conditions
        sitepairs_avg(t).difference = [];
        for diff = 1:length(lfp_tfa_cfg.diff_condition)
            diff_condition = lfp_tfa_cfg.diff_condition{diff};
            sitepairs_avg(t).difference = [sitepairs_avg(t).difference, ...
                lfp_tfa_compute_diff_condition_tfsync(sitepairs_avg(t).condition, diff_condition, 1)];
        end
        % plot Difference TFR
        for dcn = 1:length(sitepairs_avg(t).difference)
            if ~isempty(sitepairs_avg(t).difference(dcn).hs_tuned_sync)
                if isfield(sitepairs_avg(t).difference(dcn).hs_tuned_sync,... 
                        'ppc')
                    plottitle = sprintf('%s - %s (ref: %s) %s', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                        lfp_tfa_cfg.ref_hemisphere, sitepairs_avg(t).difference(dcn).label);
                    result_file = fullfile(results_fldr, ...
                        sprintf('%s-%s_%s.png', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                        sitepairs_avg(t).difference(dcn).label));
                    lfp_tfa_plot_hs_tuned_sync(sitepairs_avg(t).difference(dcn).hs_tuned_sync, ...
                        lfp_tfa_cfg, plottitle, result_file, 'cmap', 'bluewhitered', 'imscale', [-0.3, 0.3], ...
                        'significant', lfp_tfa_cfg.stat_test.plot_significant);
                end
            end
        end
        
    end
    
    % save session average tfs
    result_matfile = fullfile(results_fldr, ['Avg_LFP_LFP_sync.mat']);
    save(result_matfile, 'sitepairs_avg');
    
end