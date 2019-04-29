%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading processed LFP data, rejection of noise trials
% and task specific analysis using TFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 

% file containing settings for LFP analysis
settings_filepath = 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\LFP_timefrequency_analysis\settings\lfp_tfa_settings_v1.m';
% folder containing LFP data for analysis
%data_folder = 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data';
% maxsites 
maxsites = inf;

%% INITIALIZATION
close all;

lfp_tfa_cfg = lfp_tfa_define_settings(settings_filepath, maxsites);

%% Function calls for LFP TFA per site and per session

% Read the file list
file_list = lfp_tfa_cfg.file_list;
lfp_datafiles = file_list(:,3);

session_lfp = struct();
session_proc_lfp = struct();
lfp_tfr = struct();
lfp_evoked = struct();
lfp_pow = struct();
try
    % loop through each session for processing lfp
    for i = 1:length(lfp_datafiles)
        session_name = [file_list{i,1} '_' file_list{i,2}];
        fprintf('Processing LFP for session %s\n', session_name);
        lfp_tfa_cfg.session = session_name;
        lfp_tfa_cfg.session_results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, session_name);
        lfp_tfa_cfg.data_filepath = lfp_datafiles{i};
        % load the LFP for one session
        session_lfp(i).sites = load(lfp_tfa_cfg.data_filepath);
%         session_proc_lfp(i).monkey = file_list{i,1};
%         session_proc_lfp(i).date = file_list{i,2};
%         session_proc_lfp(i).session = session_name;
        session_proc_lfp(i).sites = lfp_tfa_process_LFP(session_lfp(i).sites, lfp_tfa_cfg);
        session_proc_lfp(i).sites = lfp_tfa_reject_noisy_lfp(session_proc_lfp(i).sites, lfp_tfa_cfg.noise);
        session_proc_lfp(i).sites = lfp_tfa_compute_baseline_power(session_proc_lfp(i).sites, lfp_tfa_cfg);
    end
    % loop through each processed session for analysis
    for i = 1:length(session_proc_lfp)
        session_name = [file_list{i,1} '_' file_list{i,2}];
        fprintf('Processing LFP for session %s\n', session_name);
        lfp_tfa_cfg.session = session_name;
        lfp_tfa_cfg.session_results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, session_name);
        lfp_tfa_cfg.data_filepath = lfp_datafiles{i};
        
        lfp_tfr.session(i) = ...
            lfp_tfa_plot_site_average_tfr( session_proc_lfp(i).sites, lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
        lfp_evoked.session(i) = ...
            lfp_tfa_plot_site_evoked_LFP( session_proc_lfp(i).sites, lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
        lfp_pow.session(i) = ...
            lfp_tfa_plot_site_powspctrum( session_proc_lfp(i).sites, lfp_tfa_cfg ) ;        
        
    end
catch e
    error(e.message());%save session_proc_lfp session_proc_lfp;
end

%% Average across sessions
if length(session_proc_lfp) > 1
    if any(strcmp(lfp_tfa_cfg.compute_avg_across, 'sessions'))
    lfp_tfr.sessions_avg = lfp_tfa_avg_tfr_across_sessions(lfp_tfr, lfp_tfa_cfg);%struct();
    lfp_evoked.sessions_avg = lfp_tfa_avg_evoked_LFP_across_sessions(lfp_evoked, lfp_tfa_cfg);
    lfp_pow.sessions_avg = lfp_tfa_avg_pow_across_sessions(lfp_pow, lfp_tfa_cfg);%struct();    
    end
    if any(strcmp(lfp_tfa_cfg.compute_avg_across, 'sites'))
        lfp_tfr.sites_avg = lfp_tfa_avg_tfr_across_sites(lfp_tfr, lfp_tfa_cfg);%struct();
        lfp_evoked.sites_avg = lfp_tfa_avg_evoked_LFP_across_sites(lfp_evoked, lfp_tfa_cfg);
        lfp_pow.sites_avg = lfp_tfa_avg_pow_across_sites(lfp_pow, lfp_tfa_cfg);
    end
end

% nsessions = length(session_proc_lfp);
% for cn = 1:length(session_proc_lfp(1).cond_based_tfs)
%     fprintf('Condition %s\n', session_proc_lfp(1).cond_based_tfs(cn).label);
%     %sessions_avg.cond_based_tfs(cn).tfs_across_sessions = struct();
%     for i = 1:length(session_proc_lfp)    
%         sessions_avg.cond_based_tfs(cn).tfs_across_sessions = struct();
%         for st = 1:size(session_proc_lfp(i).cond_based_tfs(cn).tfs_avg_session, 1)
%             for hs = 1:size(session_proc_lfp(i).cond_based_tfs(cn).tfs_avg_session, 2)
%                 if i == 1 
%                     if isfield(session_proc_lfp(i).cond_based_tfs(cn).tfs_avg_session(st, hs), 'powspctrm') ...
%                         && ~isempty(session_proc_lfp(i).cond_based_tfs(cn).tfs_avg_session(st, hs).powspctrm)
%                         sessions_avg.cond_based_tfs(cn).tfs_across_sessions(st,hs).powspctrm ...
%                             = (1/nsessions) * ...
%                             session_proc_lfp(i).cond_based_tfs(cn).tfs_avg_session(st, hs).powspctrm;
%                         sessions_avg.cond_based_tfs(cn).tfs_across_sessions(st,hs).nsites ...
%                             = session_proc_lfp(i).cond_based_tfs(cn).tfs_avg_session(st, hs).nsites;
%                         sessions_avg.cond_based_tfs(cn).tfs_across_sessions(st,hs).nsessions = nsessions;
%                     end
%                 else
%                     if isequaln(session_proc_lfp(i).cond_based_tfs(cn).cfg_condition, ...
%                             session_proc_lfp(i-1).cond_based_tfs(cn).cfg_condition)
%                         if isfield(session_proc_lfp(i).cond_based_tfs(cn).tfs_avg_session(st, hs), 'powspctrm') ...
%                             && ~isempty(session_proc_lfp(i).cond_based_tfs(cn).tfs_avg_session(st, hs).powspctrm)
%                             sessions_avg.cond_based_tfs(cn).tfs_across_sessions(st,hs).powspctrm ...
%                                 = ((1/nsessions) * ...
%                                 session_proc_lfp(i).cond_based_tfs(cn).tfs_avg_session(st, hs).powspctrm) + ...
%                                 sessions_avg.cond_based_tfs(cn).tfs_across_sessions(st,hs).powspctrm;
%                             if isfield(sessions_avg.cond_based_tfs(cn).tfs_across_sessions(st,hs), 'nsites')
%                                 sessions_avg.cond_based_tfs(cn).tfs_across_sessions(st,hs).nsites ...
%                                     = session_proc_lfp(i).cond_based_tfs(cn).tfs_avg_session(st, hs).nsites + ...
%                                     sessions_avg.cond_based_tfs(cn).tfs_across_sessions(st,hs).nsites;
%                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
%     % plot
% %     plottitle = ['Target = ' session_proc_lfp(i).sites(1).target ', Block ' ...
% %         num2str(cfg_conditions(cn).block) ', '];
% %     if cfg_conditions(cn).choice == 0
% %         plottitle = [plottitle 'Instructed trials'];
% %     else
% %         plottitle = [plottitle 'Choice trials'];
% %     end
%     plottitle = session_proc_lfp(1).cond_based_tfs(cn).label;
%     result_file = fullfile(lfp_tfa_cfg.root_results_fldr, ...
%                     ['Avg_Spectrogram_' session_proc_lfp(1).cond_based_tfs(cn).label '.png']);
%     lfp_tfa_plot_hs_tuned_tfr(sessions_avg.cond_based_tfs(cn).tfs_across_sessions, ...
%                 lfp_tfa_cfg, plottitle, result_file);
% 
%     % save session average tfs
%     %save(fullfile(lfp_tfa_cfg.root_results_fldr, 'session_average_tfs.mat'), 'cond_based_tfs');
% end
% 
% % Detect noisy trials
% %[noisy_trials] = lfp_tfa_reject_noisy_lfp(sites_lfp_folder,
% %lfp_tfa_cfg.noise);
% 
% % Compute baseline power
% %[sites_lfp_folder, baseline] = lfp_tfa_compute_baseline_power(sites_lfp_folder, lfp_tfa_cfg);
% 
% % Average TFR
% %[cond_based_tfs] = lfp_tfa_plot_average_tfr(sites_lfp_folder, lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
% 
% % Evoked LFP analysis
% %[ cond_based_evoked ] = lfp_tfa_plot_average_evoked_LFP( sites_lfp_folder, lfp_tfa_cfg.analyse_states, lfp_tfa_cfg ) ;
% 
% % LFP power spectrum
% %[ cond_based_psd ] = lfp_tfa_plot_average_powspctrum( lfp_tfa_cfg.session_results_folder, lfp_tfa_cfg);
% 
% 
