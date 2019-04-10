%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading processed LFP data, rejection of noise trials
% and task specific analysis using TFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION

clear;

close all; 

version = 1;

lfp_tfa_cfg = lfp_tfa_define_settings(version);

%% Function calls for LFP TFA per site and per session

% Read and process LFP data
lfp_datafiles = dir(fullfile(lfp_tfa_cfg.data_folder, 'sites_*.mat'));

session_lfp = struct();
session_proc_lfp = struct();
try
    % loop through each session
    for i = 1:length(lfp_datafiles)
        session_name = lfp_datafiles(i).name(7:end-4);
        fprintf('Processing LFP for session %s\n', session_name);
        lfp_tfa_cfg.session = session_name;
        lfp_tfa_cfg.session_results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, session_name);
        lfp_tfa_cfg.data_filepath = fullfile(lfp_tfa_cfg.data_folder, lfp_datafiles(i).name);
        % load the LFP for one session
        session_lfp(i).sites = load(lfp_tfa_cfg.data_filepath);
        session_proc_lfp(i).sites = lfp_tfa_process_LFP(session_lfp(i).sites, lfp_tfa_cfg);
        session_proc_lfp(i).sites = lfp_tfa_reject_noisy_lfp(session_proc_lfp(i).sites, lfp_tfa_cfg.noise);
        session_proc_lfp(i).sites = lfp_tfa_compute_baseline_power(session_proc_lfp(i).sites, lfp_tfa_cfg);
        [session_proc_lfp(i).sites, session_proc_lfp(i).cond_based_tfs] = ...
            lfp_tfa_plot_site_average_tfr( session_proc_lfp(i).sites, lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
        [ session_proc_lfp(i).sites, session_proc_lfp(i).cond_based_evoked ] = ...
            lfp_tfa_plot_site_evoked_LFP( session_proc_lfp(i).sites, lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
        [ session_proc_lfp(i).sites, session_proc_lfp(i).cond_based_psd ] = ...
            lfp_tfa_plot_site_powspctrum( session_proc_lfp(i).sites, lfp_tfa_cfg ) ;
    end
catch e
    error(e.msg());%save session_proc_lfp session_proc_lfp;
end

%% Average across sessions
%for i = 1:length(session_proc_lfp)
    

% Detect noisy trials
%[noisy_trials] = lfp_tfa_reject_noisy_lfp(sites_lfp_folder,
%lfp_tfa_cfg.noise);

% Compute baseline power
%[sites_lfp_folder, baseline] = lfp_tfa_compute_baseline_power(sites_lfp_folder, lfp_tfa_cfg);

% Average TFR
%[cond_based_tfs] = lfp_tfa_plot_average_tfr(sites_lfp_folder, lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );

% Evoked LFP analysis
%[ cond_based_evoked ] = lfp_tfa_plot_average_evoked_LFP( sites_lfp_folder, lfp_tfa_cfg.analyse_states, lfp_tfa_cfg ) ;

% LFP power spectrum
%[ cond_based_psd ] = lfp_tfa_plot_average_powspctrum( lfp_tfa_cfg.session_results_folder, lfp_tfa_cfg);


