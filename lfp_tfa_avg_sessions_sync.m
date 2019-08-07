function [results_fldr, sessions_avg] = lfp_tfa_avg_sessions_sync(sessions_info, lfp_tfa_cfg, varargin)
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
        results_fldr = varargin{1};
    else
        results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sessions', 'LFP-LFP Sync');
    end
    
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % Average TFR across sites
    sessions_avg = struct();
    % initialize struct for all target pairs and conditions
    for t = 1:length(lfp_tfa_cfg.compare.target_pairs)

        sessions_avg(t).targets = lfp_tfa_cfg.compare.target_pairs{t};

        % initialize session average struct for each condition
        for cn = 1:length(lfp_tfa_cfg.conditions)
            sessions_avg(t).condition(cn).hs_tuned_sync = struct();
            sessions_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
            sessions_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
            % initialize number of site pairs for each handspace
            % label
            for st = 1:size(lfp_tfa_cfg.analyse_states, 1)
                for hs = 1:size(lfp_tfa_cfg.conditions(1).hs_labels, 2)
                    sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).nsites = 0;
                    sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).nsessions = 0;
                    sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm = [];
                    sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.crsspctrm_abs_mean_norm = [];
                    sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.crsspctrm_ang_mean = [];
                end
            end  
        end
    end
    
    for i = 1:length(sessions_info)  
        %for j = 1:length(sessions_info(i).sitepair_info)

        % load the analysis results for this sitepair, variable
        % name is sitepairs_avg
        %fprintf('Loading average of sites across session %s\n', sessions_info(i).session);
        load(sessions_info(i).avg_sync_results, 'sitepairs_avg');  

        
%                     if ~sessions_info(i).sitepair_info(j).use_for_avg
%                         continue;
%                     end
        for t = 1:length(lfp_tfa_cfg.compare.target_pairs)
            
            if ~all(strcmp(sort(sitepairs_avg(t).targets), sort(lfp_tfa_cfg.compare.target_pairs{t})))
                continue;
            end
            
            fprintf('Target pair: %s - %s\n', lfp_tfa_cfg.compare.target_pairs{t}{:});
                 
        

            for cn = 1:length(lfp_tfa_cfg.conditions)
                fprintf('Condition %s\n', lfp_tfa_cfg.conditions(cn).label);
                
                
                    if ~isempty(sitepairs_avg(t).condition(cn).hs_tuned_sync) && ... 
                        isfield(sitepairs_avg(t).condition(cn).hs_tuned_sync, 'ppc') 
                        
                        for st = 1:size(sitepairs_avg(t).condition(cn).hs_tuned_sync, 1)
                            for hs = 1:size(sitepairs_avg(t).condition(cn).hs_tuned_sync, 2)
                                if isfield(sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs), 'ppc') ...
                                        && ~isempty(sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc) ...
                                        && isfield(sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc, 'ppcspctrm') ...
                                        && ~isempty(sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm)
                                    % increment number of session
                                    sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).nsessions = ...
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).nsessions + 1;   
                                    if sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).nsessions == 1%~isfield(sessions_avg.cond_based_tfs(cn).tfs_across_sessions, 'powspctrm')
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).hs_label ...
                                            = sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).hs_label;
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).state ...
                                            = sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).state;
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).state_name ...
                                            = sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).state_name;
                                        sessions_avg(t).condition(cn).label = ...
                                            sitepairs_avg(t).condition(cn).label;
                                        sessions_avg(t).condition(cn).cfg_condition = ...
                                            sitepairs_avg(t).condition(cn).cfg_condition;
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.time ...
                                            = sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.time;
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.freq ...
                                            = sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.freq;
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.cfg ...
                                            = sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.cfg;
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm ...
                                            = cat(1, ...
                                            sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm, ...
                                            sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm) ;
                                        % store cross spectruma and phase
                                        % difference average
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.time ...
                                            = sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.time;
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.freq ...
                                            = sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.freq;
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.cfg ...
                                            = sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.cfg;
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.crsspctrm_abs_mean_norm ...
                                            = cat(1, ...
                                            sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.crsspctrm_abs_mean_norm, ...
                                            sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.crsspctrm_abs_mean_norm) ;
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.crsspctrm_ang_mean ...
                                            = cat(1, ...
                                            sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.crsspctrm_ang_mean, ...
                                            sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.crsspctrm_ang_mean) ;                                       
                                        
                                        
                                                                                
                                        if isfield(sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs), 'nsites')
                                            sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).nsites ...
                                                = sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).nsites;
                                        end
                                    else
                                        ntimebins = size(sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm, 3);
                                        % average same number of time bins
                                        if ntimebins > length(sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.time)
                                            ntimebins = length(sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.time);
                                        end
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm ...
                                            = cat(1, sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm(:,:,1:ntimebins), (sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm(1,:,1:ntimebins))...
                                            );
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.time = sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.time(1:ntimebins);
                                        
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.crsspctrm_abs_mean_norm ...
                                            = cat(1, sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.crsspctrm_abs_mean_norm(:,:,1:ntimebins), ...
                                            (sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.crsspctrm_abs_mean_norm(1,:,1:ntimebins))...
                                            );
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.crsspctrm_ang_mean ...
                                            = cat(1, sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.crsspctrm_ang_mean(:,:,1:ntimebins), ...
                                            (sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.crsspctrm_ang_mean(1,:,1:ntimebins))...
                                            );
                                        sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).csd.time = ...
                                            sessions_avg(t).condition(cn).hs_tuned_sync(st, hs).ppc.time(1:ntimebins);
                                        if isfield(sitepairs_avg(t).condition(cn).hs_tuned_sync(st,hs), 'nsites')
                                            sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).nsites ...
                                                = sitepairs_avg(t).condition(cn).hs_tuned_sync(st, hs).nsites + ...
                                                sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).nsites;
                                        end                               
                                    end
                                end
                            end
                        end
                    end
                %end
            end
        end
    end

    for t = 1:length(lfp_tfa_cfg.compare.target_pairs)
        % compute average and plot
        for cn = 1:length(sessions_avg(t).condition)
            if ~isempty(sessions_avg(t).condition(cn).hs_tuned_sync)
                if isfield(sessions_avg(t).condition(cn).hs_tuned_sync, 'ppc')
                    for st = 1:size(sessions_avg(t).condition(cn).hs_tuned_sync, 1)
                        for hs = 1:size(sessions_avg(t).condition(cn).hs_tuned_sync, 2)
                            % average synchronization measure
                            sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).ppc.ppcspctrm_all = ...
                                sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).ppc.ppcspctrm;
                            sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).ppc.ppcspctrm_std = ...
                                nanstd( ...
                                sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).ppc.ppcspctrm_all, 0, 1);
                            sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).ppc.ppcspctrm = ...
                                nanmean( ...
                                sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).ppc.ppcspctrm_all, 1); 
                            % average cross spectrum and phase difference
                            sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).csd.crsspctrm_abs_mean_norm_all = ...
                                sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).csd.crsspctrm_abs_mean_norm;
                            sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).csd.crsspctrm_abs_mean_norm = ...
                                nanmean( ...
                                sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).csd.crsspctrm_abs_mean_norm_all, 1);
                            sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).csd.crsspctrm_ang_mean = ...
                                nanmean( ...
                                sessions_avg(t).condition(cn).hs_tuned_sync(st,hs).csd.crsspctrm_ang_mean, 1);
                        end
                    end
                end
            end


            if ~isempty(sessions_avg(t).condition(cn).hs_tuned_sync)
                if isfield(sessions_avg(t).condition(cn).hs_tuned_sync,... 
                        'ppc')
                    plottitle = sprintf('%s - %s (ref: %s) %s', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                        lfp_tfa_cfg.ref_hemisphere, sessions_avg(t).condition(cn).label);
                    result_file = fullfile(results_fldr, ...
                        sprintf('%s-%s_%s.png', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                        sessions_avg(t).condition(cn).label));
                    lfp_tfa_plot_hs_tuned_sync(sessions_avg(t).condition(cn).hs_tuned_sync, ...
                        lfp_tfa_cfg, plottitle, result_file, 'imscale', [0, 1]);
                    
                    plottitle_crs = sprintf('LFP-LFP crossspectrum %s - %s (ref: %s) %s', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                        lfp_tfa_cfg.ref_hemisphere, sessions_avg(t).condition(cn).label);
                    result_file_crs = fullfile(results_fldr, ...
                        sprintf('LFP-LFP_crsspctrm_%s-%s_%s.png', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                        sessions_avg(t).condition(cn).label));
                    lfp_tfa_plot_hs_tuned_crsspctrm(sessions_avg(t).condition(cn).hs_tuned_sync, ...
                        lfp_tfa_cfg, plottitle_crs, result_file_crs, 'imscale', [-1, 1]);
                    
                    plottitle = sprintf('LFP-LFP phase diff %s - %s (ref: %s) %s', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                        lfp_tfa_cfg.ref_hemisphere, sessions_avg(t).condition(cn).label);
                    result_file = fullfile(results_fldr, ...
                        sprintf('LFP-LFP_phase_%s-%s_%s.png', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                        sessions_avg(t).condition(cn).label));
                    lfp_tfa_plot_hs_tuned_phasediff(sessions_avg(t).condition(cn).hs_tuned_sync, ...
                        lfp_tfa_cfg, plottitle, result_file, 'imscale', [-pi/2, pi/2]);
                end
            end        
        end

        % average difference between conditions
        %if sum(lfp_tfa_cfg.compare.perturbations == [0, 1]) > 1
            %sites_avg(t).difference = lfp_tfa_compute_diff_condition_tfr(sites_avg(t), lfp_tfa_cfg.difference);
            sessions_avg(t).difference = [];
            for diff = 1:length(lfp_tfa_cfg.diff_condition)
                diff_condition = lfp_tfa_cfg.diff_condition{diff};
                sessions_avg(t).difference = [sessions_avg(t).difference, ...
                    lfp_tfa_compute_diff_condition_tfsync(sessions_avg(t).condition, diff_condition)];
            end
            
            % plot Difference TFR
            for dcn = 1:length(sessions_avg(t).difference)
                if ~isempty(sessions_avg(t).difference(dcn).hs_tuned_sync)
                    if isfield(sessions_avg(t).difference(dcn).hs_tuned_sync,... 
                            'ppc')
                        plottitle = sprintf('LFP-LFP sync %s - %s (ref: %s) %s', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                            lfp_tfa_cfg.ref_hemisphere, sessions_avg(t).difference(dcn).label);
                        result_file = fullfile(results_fldr, ...
                            sprintf('LFP-LFP_Diffsync_%s-%s_%s.png', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                            sessions_avg(t).difference(dcn).label));
                        lfp_tfa_plot_hs_tuned_sync(sessions_avg(t).difference(dcn).hs_tuned_sync, ...
                            lfp_tfa_cfg, plottitle, result_file, 'cmap', 'bluewhitered', 'imscale', [-0.3, 0.3]);
                        
                        plottitle_crs = sprintf('LFP-LFP cross spectrum %s - %s (ref: %s) %s', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                            lfp_tfa_cfg.ref_hemisphere, sessions_avg(t).difference(dcn).label);
                        result_file_crs = fullfile(results_fldr, ...
                            sprintf('LFP-LFP_Diffcrsspctrm_%s-%s_%s.png', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                            sessions_avg(t).difference(dcn).label));
                        lfp_tfa_plot_hs_tuned_crsspctrm(sessions_avg(t).difference(dcn).hs_tuned_sync, ...
                            lfp_tfa_cfg, plottitle_crs, result_file_crs, 'cmap', 'bluewhitered', 'imscale', [-1, 1]);
                        
                        plottitle_phase = sprintf('LFP-LFP phase diff %s - %s (ref: %s) %s', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                            lfp_tfa_cfg.ref_hemisphere, sessions_avg(t).difference(dcn).label);
                        result_file_phase = fullfile(results_fldr, ...
                            sprintf('LFP-LFP_Diffphase_%s-%s_%s.png', lfp_tfa_cfg.compare.target_pairs{t}{:}, ...
                            sessions_avg(t).difference(dcn).label));
                        lfp_tfa_plot_hs_tuned_phasediff(sessions_avg(t).difference(dcn).hs_tuned_sync, ...
                            lfp_tfa_cfg, plottitle_phase, result_file_phase, 'cmap', 'bluewhitered', 'imscale', [-pi/2, pi/2]);
                    end
                end
            end
        %end
    end
    
    % save session average tfs
    save(fullfile(results_fldr, ['Avg_LFP_LFP_sync.mat']), 'sessions_avg');
    
end