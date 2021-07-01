function [result_matfile, sitepairs_avg] = lfp_tfa_avg_sitepairs_syncspctrm(sessions_info, lfp_tfa_cfg, varargin)
%lfp_tfa_avg_sessions_syncspctrm  - Condition-based LFP-LFP phase synchronization
%spectral average across many sitepair averages (A sitepair average is inturn an
%average across multiple trials recorded in a pair of sites)
%
% USAGE:
%	[result_matfile, sitepairs_avg] = lfp_tfa_avg_sessions_syncspctrm
%   (sessions_info, lfp_tfa_cfg)
%   [result_matfile, sitepairs_avg] = lfp_tfa_avg_sessions_syncspctrm
%   (sessions_info, lfp_tfa_cfg, results_folder)
%
% INPUTS:
%		sessions_info   - 1xM struct containing information about the sessions
%		to be analysed, see output of lfp_tfa_define_settings
%       Required fields: 
%           lfp_syncspctrm_results_fldr - absolute path to the folder containing LFP-LFP
%           phase synchronization spectral average results for each of the sitepairs within
%           a session, see lfp_tfa_sitepair_averaged_syncspctrm 
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields:
%               1. conditions               - trial conditions to compare, 
%               see settings/lfp_tfa_settings_example.m and
%               lfp_tfa_compare_conditions.m
%               2. root_results_fldr        - root folder where results are
%               saved. Results will be stored under 
%               [lfp_tfa_cfg.root_results_fldr
%               '/Avg_across_sessions/LFP-LFP Syncspectrum']
%               3. compare.target_pairs     - targets to compare, see
%               settings/lfp_tfa_settings_example
%               4. analyse_epochs           - epochs of interest, see
%               settings/lfp_tfa_settings_example
%               5. ref_hemisphere           - reference hemisphere for 
%               contra- and ipsi- labelling, see
%               settings/lfp_tfa_settings_example
% OUTPUTS:
%       result_matfile   - absolute path to the file where the LFP-LFP sync
%       spectral average across multiple sitepairs would be stored
%		sitepairs_avg   - structure containing LFP-LFP sync spectral average across
%       multiple sitepair averages
%
% REQUIRES:	lfp_tfa_plot_hs_tuned_syncsp
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_define_settings, 
% lfp_tfa_compare_conditions, lfp_tfa_sitepair_averaged_syncspctrm
% lfp_tfa_plot_hs_tuned_syncsp, lfp_tfa_avg_sessions_syncspctrm
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
        results_fldr = fullfile(sessions_info.lfp_syncspctrm_results_fldr, 'session_avg');
    else
        results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sites', 'LFP_LFP_Syncspectrum');
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
            sitepairs_avg(t).condition(cn).hs_tuned_syncsp = struct();
            sitepairs_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
            sitepairs_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
            % initialize number of site pairs for each handspace
            % label
            for ep = 1:size(lfp_tfa_cfg.analyse_epochs, 1)
                for hs = 1:size(lfp_tfa_cfg.conditions(1).hs_labels, 2)
                    sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep, hs).nsites = 0;
                    sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep, hs).ppc.ppcspctrm = [];
                end
            end                        
        end
    end
    for i = 1:length(sessions_info) 
        sitepairs_fldr = dir(sessions_info(i).lfp_syncspctrm_results_fldr);
        sitepairs_fldr(1:2) = [];
        for j = 1:length(sitepairs_fldr)
%                     if ~sessions_info(i).sitepair_info(j).use_for_avg
%                         continue;
%                     end
            % load the analysis results for this sitepair, variable
            % name is sitepair_syncspctrm
            %fprintf('Loading site pair %s - %s\n', sessions_info(i).sitepair_info(j).sites{:});
            if strcmp(sitepairs_fldr(j).name, 'session_avg')
                continue;
            end
            sitepair_sync_file = dir(fullfile(...
                sessions_info(i).lfp_syncspctrm_results_fldr, sitepairs_fldr(j).name, '*.mat'));
            sitepair_sync_filepath = fullfile(sessions_info(i).lfp_syncspctrm_results_fldr, ...
                sitepairs_fldr(j).name, sitepair_sync_file(1).name);
            load(sitepair_sync_filepath, 'sitepair_syncspctrm');
            
            t = [];
            for t_idx = 1:length(sitepairs_avg)
                if find(all(strcmp(sort(sitepair_syncspctrm.targets), ...
                        sort(sitepairs_avg(t_idx).targets))))
                    t = t_idx;
                    fprintf('Target pair: %s - %s\n', sitepair_syncspctrm.targets{1}, sitepair_syncspctrm.targets{2});
                    break;
                end
            end

            if isempty(t)%~all(strcmp(sort(sitepair_syncspctrm.targets), sort(lfp_tfa_cfg.compare.target_pairs{t})))
                continue;
            end

            sitepair_syncspctrm.use_for_avg = 1;
            for cn = 1:length(lfp_tfa_cfg.conditions)
                sitepair_syncspctrm.use_for_avg = ...
                    ~any(sitepair_syncspctrm.condition(cn).ntrials - sitepair_syncspctrm.condition(cn).noisytrials < 5);
            end
            if sitepair_syncspctrm.use_for_avg == 0
                continue;
            end
            % compute condition-based average across sitepairs
            for cn = 1:length(lfp_tfa_cfg.conditions)
                fprintf('Condition %s\n', lfp_tfa_cfg.conditions(cn).label);     

                if ~isempty(sitepair_syncspctrm.condition(cn).hs_tuned_syncsp) && ... 
                    isfield(sitepair_syncspctrm.condition(cn).hs_tuned_syncsp, 'ppc') 

                    for ep = 1:size(sitepair_syncspctrm.condition(cn).hs_tuned_syncsp, 1)
                        for hs = 1:size(sitepair_syncspctrm.condition(cn).hs_tuned_syncsp, 2)
                            if isfield(sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs), 'ppc') ...
                                    && ~isempty(sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).ppc) ...
                                    && isfield(sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).ppc, 'ppcspctrm') ...
                                    && ~isempty(sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).ppc.ppcspctrm)
                                % increment number of site pairs
                                sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep, hs).nsites = ...
                                    sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep, hs).nsites + 1;
                                if sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep, hs).nsites == 1%~isfield(sessions_avg.cond_based_tfs(cn).tfs_across_sessions, 'powspctrm')
                                    sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep,hs).hs_label ...
                                        = sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).hs_label;
                                    sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep,hs).epoch_name ...
                                        = sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).epoch_name;
                                    sitepairs_avg(t).condition(cn).label = ...
                                        sitepair_syncspctrm.condition(cn).label;
                                    sitepairs_avg(t).condition(cn).cfg_condition = ...
                                        sitepair_syncspctrm.condition(cn).cfg_condition;
%                                         sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep, hs).ppc.time ...
%                                             = sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).ppc.time;
                                    sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep, hs).ppc.freq ...
                                        = sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).ppc.freq;                                        
                                    sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep, hs).ppc.cfg ...
                                        = sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).ppc.cfg;
                                    sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep, hs).ppc.ppcspctrm ...
                                        = [sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep, hs).ppc.ppcspctrm; sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).ppc.ppcspctrm];

                                else
                                    sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep, hs).ppc.ppcspctrm ...
                                        = [sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep, hs).ppc.ppcspctrm; ...
                                        (sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).ppc.ppcspctrm) ...
                                        ];
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

    for t = 1:length(sitepairs_avg)
        % compute average and plot
        for cn = 1:length(sitepairs_avg(t).condition)
            if ~isempty(sitepairs_avg(t).condition(cn).hs_tuned_syncsp)
                if isfield(sitepairs_avg(t).condition(cn).hs_tuned_syncsp, 'ppc')
                    for ep = 1:size(sitepairs_avg(t).condition(cn).hs_tuned_syncsp, 1)
                        for hs = 1:size(sitepairs_avg(t).condition(cn).hs_tuned_syncsp, 2)
                            if isempty(sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep,hs).ppc)
                                continue;
                            end
%                             sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep,hs).ppc.ppcspctrm_std = ...
%                                 nanstd( ...
%                                 sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep,hs).ppc.ppcspctrm, 0, 1);
                            [sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep,hs).ppc.mean, ...
                                sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep,hs).ppc.ppcspctrm_std] = ...
                                lfp_tfa_compute_statistics( ...
                                sitepairs_avg(t).condition(cn).hs_tuned_syncsp(ep,hs).ppc.ppcspctrm, lfp_tfa_cfg.error_measure);
                        end
                    end
                end
            end


            if ~isempty(sitepairs_avg(t).condition(cn).hs_tuned_syncsp)
                if isfield(sitepairs_avg(t).condition(cn).hs_tuned_syncsp,... 
                        'ppc')
                    plottitle = sprintf(lfp_tfa_cfg.monkey, 'LFP-LFP sync spectrum %s - %s (ref: %s) %s', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                        lfp_tfa_cfg.ref_hemisphere, sitepairs_avg(t).condition(cn).label);
                    result_file = fullfile(results_fldr, ...
                        sprintf(lfp_tfa_cfg.monkey, '%s-%s_%s', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                        sitepairs_avg(t).condition(cn).label));
                    lfp_tfa_plot_hs_tuned_syncsp(sitepairs_avg(t).condition(cn).hs_tuned_syncsp, ...
                        lfp_tfa_cfg, plottitle, result_file);
                end
            end        
        end
      
        
    end
    
    % save session average tfs
    result_matfile = fullfile(results_fldr, [lfp_tfa_cfg.monkey 'Avg_LFP_LFP_syncspctrm.mat']);
    save(result_matfile, sitepairs_avg');
    
end