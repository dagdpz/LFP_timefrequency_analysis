%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading processed LFP data, rejection of noise trials
% and task specific analysis using TFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 

% file containing settings for LFP analysis
settings_filepath = 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\LFP_timefrequency_analysis\settings\lfp_tfa_settings_v1.m';

%% INITIALIZATION
close all;

lfp_tfa_cfg = lfp_tfa_define_settings(settings_filepath);

%% Function calls for LFP TFA per site and per session

% Read the info about sessions to analyse
sessions_info = lfp_tfa_cfg.session_info;
% Get the path to the mat file containing LFP data
lfp_datafiles = {sessions_info.Input};

% initialize structs to store intermediate results
% struct to read in LFP data
session_lfp = struct();
% struct to store processed LFP data
session_proc_lfp = struct();
% struct to store average LFP TFR for different conditions
lfp_tfr = struct();
% struct to store average LFP evoked response for different conditions
lfp_evoked = struct();
% struct to store average LFP power spectrum for different conditions
lfp_pow = struct();

try
    %% loop through each session for processing lfp
    for i = 1:length(sessions_info)
        % name of session = [Monkey name '_' Recording date]
        session_name = [sessions_info(i).Monkey '_' sessions_info(i).Date];
        fprintf('Processing LFP for session %s\n', session_name);
        lfp_tfa_cfg.session = session_name;
        % folder to which results of analysis of this session should be
        % stored
        lfp_tfa_cfg.session_results_fldr = ...
            fullfile(lfp_tfa_cfg.root_results_fldr, session_name);
        % absolute path of file containing LFP data for this session
        lfp_tfa_cfg.data_filepath = lfp_datafiles{i};
        % load the processed LFP data for one session
        session_lfp(i).sites = load(lfp_tfa_cfg.data_filepath);
        % read LFP data for each site and each trial and calculate the 
        % trial-wise time frequency spectrogram
        session_proc_lfp(i).sites = ...
            lfp_tfa_process_LFP(session_lfp(i).sites, lfp_tfa_cfg);
        % Detect the noisy trials for each site of a session
        if lfp_tfa_cfg.noise.detect
            session_proc_lfp(i).sites = ...
                lfp_tfa_reject_noisy_lfp(session_proc_lfp(i).sites, ...
                lfp_tfa_cfg.noise);
        end
        % Calculate the baseline power for each site
        session_proc_lfp(i).sites = ...
            lfp_tfa_compute_baseline_power(session_proc_lfp(i).sites, ...
            lfp_tfa_cfg);
    end
    %% loop through each processed session for analysis
    for i = 1:length(session_proc_lfp)
        session_name = [sessions_info(i).Monkey '_' sessions_info(i).Date];
        fprintf('Processing LFP for session %s\n', session_name);
        lfp_tfa_cfg.session = session_name;
        % folder to which results of analysis of this session should be
        % stored
        lfp_tfa_cfg.session_results_fldr = ...
            fullfile(lfp_tfa_cfg.root_results_fldr, session_name);
        % absolute path of file containing LFP data for this session
        lfp_tfa_cfg.data_filepath = lfp_datafiles{i};
        % preinjection blocks for this session
        lfp_tfa_cfg.perturbation_groups = cell(1, 2);
        if isfield(sessions_info, 'Preinj_blocks') && ...
                ~isempty(sessions_info(i).Preinj_blocks)
            lfp_tfa_cfg.perturbation_groups{1} = sessions_info(i).Preinj_blocks;
        elseif isfield(lfp_tfa_cfg.compare, 'perturbation_groups')
            lfp_tfa_cfg.perturbation_groups{1} = lfp_tfa_cfg.compare.perturbation_groups(1);
        else
            lfp_tfa_cfg.perturbation_groups{1} = 0;            
        end
        % postinjection blocks for this session
        if isfield(sessions_info, 'Postinj_blocks') && ...
                ~isempty(sessions_info(i).Postinj_blocks)
            lfp_tfa_cfg.perturbation_groups{2} = sessions_info(i).Postinj_blocks;
        elseif isfield(lfp_tfa_cfg.compare, 'perturbation_groups')
            lfp_tfa_cfg.perturbation_groups{2} = lfp_tfa_cfg.compare.perturbation_groups(2);
        else
            lfp_tfa_cfg.perturbation_groups{2} = 'all';            
        end
        % Calculate and plot the site-wise and session average TFR, 
        % evoked response and power spectral density for 
        %different conditions and hand-space labels 
        lfp_tfr.session(i) = ...
            lfp_tfa_plot_site_average_tfr( session_proc_lfp(i).sites, ...
            lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
        lfp_evoked.session(i) = ...
            lfp_tfa_plot_site_evoked_LFP( session_proc_lfp(i).sites, ...
            lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
        lfp_pow.session(i) = ...
            lfp_tfa_plot_site_powspctrum( session_proc_lfp(i).sites, ...
            lfp_tfa_cfg ) ;        
        
    end
catch e
    error(e.message());
end

%% Average across sessions
if length(session_proc_lfp) > 1
    % average session averages
    if any(strcmp(lfp_tfa_cfg.compute_avg_across, 'sessions'))
        % Average of session averages of LFP TFR, LFP evoked responce and
        % LFP power spectrum response
    lfp_tfr.sessions_avg = ...
        lfp_tfa_avg_tfr_across_sessions(lfp_tfr, lfp_tfa_cfg);
    lfp_evoked.sessions_avg = ...
        lfp_tfa_avg_evoked_LFP_across_sessions(lfp_evoked, lfp_tfa_cfg);
    lfp_pow.sessions_avg = ...
        lfp_tfa_avg_pow_across_sessions(lfp_pow, lfp_tfa_cfg);
    end
    % average site averages
    if any(strcmp(lfp_tfa_cfg.compute_avg_across, 'sites'))
        % Average of site averages of LFP TFR, LFP evoked response and LFP
        % power spectrum
        lfp_tfr.sites_avg = ...
            lfp_tfa_avg_tfr_across_sites(lfp_tfr, lfp_tfa_cfg);
        lfp_evoked.sites_avg = ...
            lfp_tfa_avg_evoked_LFP_across_sites(lfp_evoked, lfp_tfa_cfg);
        lfp_pow.sites_avg = ...
            lfp_tfa_avg_pow_across_sites(lfp_pow, lfp_tfa_cfg);
    end
end

