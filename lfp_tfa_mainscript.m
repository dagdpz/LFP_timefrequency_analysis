%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading LFP data, rejection of noise trials
% and condition specific analysis using TFR, evoked, spectra, sync
% spectrograms and sync spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;

% file containing settings for LFP analysis
% should have the same format as settings/lfp_tfa_settings_example.m
%   settings_filepath = 'C:\Users\mpachoud\Documents\GitHub\LFP_timefrequency_analysis\settings\PPC_pulv_eye_hand\Linus\Linus_dPul_LIP_inactivation_combined.m';
   %settings_filepath = 'C:\Users\lschneider\GitHub\Settings\LFP_time_frequency_analysis\Pulv_eye_hand\Interleaved\lfp_tfa_settings.m';
   settings_filepath = 'C:\Users\lschneider\GitHub\Settings\LFP_time_frequency_analysis\Pulv_oculomotor\paper\lfp_tfa_settings.m';

%% INITIALIZATION
close all;

% read settings file
lfp_tfa_cfg = lfp_tfa_define_settings(settings_filepath);

%% Get info about sessions to be analysed

% Read the info about sessions to analyse
sessions_info = lfp_tfa_cfg.session_info;
% Get the path to the mat file containing LFP data
lfp_datafiles = {sessions_info.Input};

%% initialize structs to store intermediate results
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


%% LFP processing
% loop through each session to process
for i = 1:length(sessions_info)
    % name of session = [Monkey name '_' Recording date]
    session_name = [sessions_info(i).Monkey '_' sessions_info(i).Date];
    lfp_tfa_cfg.session = session_name;
    % folder to which results of analysis of this session should be
    % stored
    sessions_info(i).proc_results_fldr = ...
        fullfile(lfp_tfa_cfg.proc_lfp_folder, session_name);
    % absolute path of file containing LFP data for this session
    % check if folder exists and if there are processed LFP files inside
    if ~exist(sessions_info(i).proc_results_fldr, 'dir') || ...
            isempty(dir(fullfile(sessions_info(i).proc_results_fldr, '*.mat')))
        warning('No mat files with processed LFP found for session %s\n'...
            , session_name);
        fprintf('Processing LFP for session %s\n', session_name);
        sessions_info(i) = ...
            lfp_tfa_process_LFP(sessions_info(i), lfp_tfa_cfg);
    elseif lfp_tfa_cfg.process_LFP
        % read LFP data for each site and each trial and calculate the
        % trial-wise time frequency spectrogram
        fprintf('Processing LFP for session %s\n', session_name);
        sessions_info(i) = ...
            lfp_tfa_process_LFP(sessions_info(i), lfp_tfa_cfg);
    end
    
end
lfp_tfa_cfg.sessions_info = sessions_info;
%% loop through each processed session for analysis
for i = 1:length(sessions_info)
    
    % clear variables
    %clear session_proc_lfp session_lfp;
    
    session_name = [sessions_info(i).Monkey '_' sessions_info(i).Date];
    fprintf('Analysing LFP for session %s\n', session_name);
    lfp_tfa_cfg.session = session_name;
    
    % folder to which results of analysis of this session should be
    % stored
    sessions_info(i).analyse_lfp_fldr = ...
        fullfile(lfp_tfa_cfg.analyse_lfp_folder, session_name);
    lfp_tfa_cfg.session_results_fldr = ...
        fullfile(lfp_tfa_cfg.analyse_lfp_folder, session_name);
    
    
    % read the processed lfp mat files for all sites of this session
    sites_lfp_files = dir(fullfile(sessions_info(i).proc_results_fldr, '*.mat'));
    session_proc_lfp = [];
    % first check if the site averages need to be calculated
    if (any(strcmp(lfp_tfa_cfg.analyses, 'tfs')) && ...
            (~exist(sessions_info(i).lfp_tfs_results_fldr, 'dir') || ...
            isempty(dir(fullfile(sessions_info(i).lfp_tfs_results_fldr, 'LFP_TFR*.mat'))))) || ...
            (any(strcmp(lfp_tfa_cfg.analyses, 'evoked')) && ...
            (~exist(sessions_info(i).lfp_evoked_results_fldr, 'dir')|| ...
            isempty(dir(fullfile(sessions_info(i).lfp_evoked_results_fldr, 'LFP_Evoked*.mat'))))) || ...
            (any(strcmp(lfp_tfa_cfg.analyses, 'pow')) && ...
            (~exist(sessions_info(i).lfp_pow_results_fldr, 'dir') || ...
            isempty(dir(fullfile(sessions_info(i).lfp_pow_results_fldr, 'LFP_Power*.mat')))))
        lfp_tfa_cfg.compute_site_average = true;
    end
    % if site averages need to be computed, load the processed LFP
    if lfp_tfa_cfg.compute_site_average
        for file = {sites_lfp_files.name}
            if ~strcmp(file{1}, 'allsites_lfp.mat')
                fprintf('Reading LFP data from %s\n', file{1});
                load(fullfile(sessions_info(i).proc_results_fldr, file{1}))
                session_proc_lfp = [session_proc_lfp site_lfp];
            end
        end
    end
    
    lfp_tfa_cfg.data_filepath = lfp_datafiles{i};
    
    % trial conditions to analyse
    conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg, ...
        {sessions_info(i).Preinj_blocks, sessions_info(i).Postinj_blocks});
    
    % Calculate and plot the site-wise and session average TFR,
    % evoked response and power spectral density for
    % different conditions and hand-space labels
    if any(strcmp(lfp_tfa_cfg.analyses, 'tfs'))
        % check if site-wise average needs to be computed
        if lfp_tfa_cfg.compute_site_average
            lfp_tfr.session(i) = ...
                lfp_tfa_plot_site_average_tfr( session_proc_lfp, ...
                conditions, lfp_tfa_cfg );
        else
            % load pre-computed results if available
            lfp_tfr_files = dir(fullfile(...
                sessions_info(i).lfp_tfs_results_fldr, 'LFP_TFR*.mat'));
            load(fullfile(sessions_info(i).lfp_tfs_results_fldr, ...
                lfp_tfr_files(1).name), 'session_tfs');
            lfp_tfr.session(i) = session_tfs;
        end
    end
    if any(strcmp(lfp_tfa_cfg.analyses, 'evoked'))
        % check if site-wise average needs to be computed
        if lfp_tfa_cfg.compute_site_average
            lfp_evoked.session(i) = ...
                lfp_tfa_plot_site_evoked_LFP( session_proc_lfp, ...
                conditions, lfp_tfa_cfg );
        else
            % load pre-computed results if available
            lfp_evoked_files = dir(fullfile(...
                sessions_info(i).lfp_evoked_results_fldr, 'LFP_Evoked*.mat'));
            load(fullfile(sessions_info(i).lfp_evoked_results_fldr, ...
                lfp_evoked_files(1).name), 'session_evoked');
            lfp_evoked.session(i) = session_evoked;
        end
    end
    if any(strcmp(lfp_tfa_cfg.analyses, 'pow'))
        % check if site-wise average needs to be computed
        if lfp_tfa_cfg.compute_site_average
            lfp_pow.session(i) = ...
                lfp_tfa_plot_site_powspctrum( session_proc_lfp, ...
                conditions, lfp_tfa_cfg ) ;
        else
            % load pre-computed results if available
            lfp_pow_files = dir(fullfile(...
                sessions_info(i).lfp_pow_results_fldr, 'LFP_Power*.mat'));
            load(fullfile(sessions_info(i).lfp_pow_results_fldr, ...
                lfp_pow_files(1).name), 'session_pow');
            lfp_pow.session(i) = session_pow;
        end
    end
    
    % clear variables
    clear session_proc_lfp session_lfp;
    
    if any(strcmp(lfp_tfa_cfg.analyses, 'sync')) ...
            || any(strcmp(lfp_tfa_cfg.analyses, 'syncsp'))
        % session cross power spctrum folder
        session_csd_folder = ...
            fullfile(sessions_info(i).proc_results_fldr, 'crossspectrum');
        % loop through each sitepair
        for sitepair_csd_file = dir(fullfile(session_csd_folder, '*.mat'))'
            clear sitepair_crosspow;
            sitepair_sync_folder = ...
                fullfile(sessions_info(i).lfp_sync_results_fldr, ...
                sitepair_csd_file.name(16:end-4));
            sitepair_syncsp_folder = ...
                fullfile(sessions_info(i).lfp_syncspctrm_results_fldr, ...
                sitepair_csd_file.name(16:end-4));
            if lfp_tfa_cfg.compute_site_average || ...
                    ~exist(sitepair_sync_folder, 'dir') || ...
                    isempty(dir(fullfile(sitepair_sync_folder, 'sitepair_sync*.mat'))) || ...
                    ~exist(sitepair_syncsp_folder, 'dir') || ...
                    isempty(dir(fullfile(sitepair_syncsp_folder, 'LFP-LFP_syncspctrm*.mat')))
                % load mat file
                fprintf('Reading LFP-LFP cross spectrum from %s\n', ...
                    sitepair_csd_file.name);
                load(fullfile(session_csd_folder, ...
                    sitepair_csd_file.name), 'sitepair_crosspow');
                
                if any(strcmp(lfp_tfa_cfg.analyses, 'sync'))
                    if lfp_tfa_cfg.compute_site_average || ...
                            ~exist(sitepair_sync_folder, 'dir') || ...
                            isempty(dir(fullfile(sitepair_sync_folder, 'sitepair_sync*.mat')))
                        % compute ppc spectrogram between sitepair
                        % get the trial conditions for this session
                        sitepair_sync = lfp_tfa_sitepair_averaged_sync(...
                            sitepair_crosspow, conditions, lfp_tfa_cfg);
                    end
                end
                
                if any(strcmp(lfp_tfa_cfg.analyses, 'syncsp'))
                    if lfp_tfa_cfg.compute_site_average || ...
                            ~exist(sitepair_syncsp_folder, 'dir') || ...
                            isempty(dir(fullfile(sitepair_syncsp_folder, 'LFP-LFP_syncspctrm*.mat')))
                        % compute ppc spectrogram between sitepair
                        % get the trial conditions for this session
                        sitepair_syncsp = lfp_tfa_sitepair_averaged_syncspctrm(...
                            sitepair_crosspow, conditions, lfp_tfa_cfg);
                    end
                end
                
            end
        end
        
    end
    
    % Calculate the session-wise average of LFP-LFP phase sync
    if any(strcmp(lfp_tfa_cfg.analyses, 'sync')) && ...
            any(strcmp(lfp_tfa_cfg.compute_avg_across, 'sessions'))
        session_avg_sync_folder = ...
            fullfile(sessions_info(i).lfp_sync_results_fldr, ...
            'session_avg');
        if lfp_tfa_cfg.compute_site_average || ...
                ~exist(session_avg_sync_folder, 'dir') || ...
                isempty(dir(fullfile(session_avg_sync_folder, 'Avg_LFP_LFP_sync.mat')))
            sessions_info(i).avg_sync_results = lfp_tfa_avg_sitepairs_sync(...
                sessions_info(i), lfp_tfa_cfg);
        else
            sessions_info(i).avg_sync_results = fullfile(...
                session_avg_sync_folder, 'Avg_LFP_LFP_sync.mat');
        end
    end
    
    % Calculate the session-wise average of LFP-LFP phase sync spectrum
    if any(strcmp(lfp_tfa_cfg.analyses, 'syncsp')) && ...
            any(strcmp(lfp_tfa_cfg.compute_avg_across, 'sessions'))
        session_avg_syncsp_folder = ...
            fullfile(sessions_info(i).lfp_syncspctrm_results_fldr, ...
            'session_avg');
        if lfp_tfa_cfg.compute_site_average || ...
                ~exist(session_avg_syncsp_folder, 'dir') || ...
                isempty(dir(fullfile(session_avg_syncsp_folder, 'Avg_LFP_LFP_syncspctrm.mat')))
            sessions_info(i).avg_syncspctrm_results = ...
                lfp_tfa_avg_sitepairs_syncspctrm(sessions_info(i), lfp_tfa_cfg);
        else
            sessions_info(i).avg_syncspctrm_results = fullfile(...
                session_avg_syncsp_folder, 'Avg_LFP_LFP_syncspctrm.mat');
        end
    end
    
end


%% Average across sessions
%% do monkey separation or combination here!?

for m=1:numel(lfp_tfa_cfg.monkeys)
    if isempty(lfp_tfa_cfg.monkeys{m}) %combined
        lfp_tfa_cfg.monkey='';
        m_idx=true(size(sessions_info));
    else
        lfp_tfa_cfg.monkey=[lfp_tfa_cfg.monkeys{m} '_'];
        m_idx=ismember({sessions_info.Monkey},lfp_tfa_cfg.monkeys{m});
    end
    if sum(m_idx) > 1
        % average session averages
        if any(strcmp(lfp_tfa_cfg.compute_avg_across, 'sessions'))
            % Average of session averages of LFP TFR, LFP evoked responce and
            % LFP power
            % LFP TFR
            if any(strcmp(lfp_tfa_cfg.analyses, 'tfs'))
                lfp_tfr.sessions_avg(m) = ...
                    lfp_tfa_avg_tfr_across_sessions(lfp_tfr.session(m_idx), lfp_tfa_cfg);
            end
            % LFP evoked response
            if any(strcmp(lfp_tfa_cfg.analyses, 'evoked'))
                lfp_evoked.sessions_avg(m) = ...
                    lfp_tfa_avg_evoked_LFP_across_sessions(lfp_evoked.session(m_idx), lfp_tfa_cfg);
            end
            % LFP spectral power
            if any(strcmp(lfp_tfa_cfg.analyses, 'pow'))
                lfp_pow.sessions_avg(m) = ...
                    lfp_tfa_avg_pow_across_sessions(lfp_pow.session(m_idx), lfp_tfa_cfg);
            end
            % LFP-LFP phase sync
            if any(strcmp(lfp_tfa_cfg.analyses, 'sync'))
                lfp_tfa_avg_sessions_sync(sessions_info(m_idx), lfp_tfa_cfg);
            end
            % LFP-LFP phase sync
            if any(strcmp(lfp_tfa_cfg.analyses, 'syncsp'))
                lfp_tfa_avg_sessions_syncspctrm(sessions_info(m_idx), lfp_tfa_cfg);
            end
        end
        % average site averages
        if any(strcmp(lfp_tfa_cfg.compute_avg_across, 'sites'))
            % Average of site averages of LFP TFR, LFP evoked response and LFP
            % power spectrum
            if any(strcmp(lfp_tfa_cfg.analyses, 'tfs'))
                lfp_tfr.sites_avg(m) = ...
                    lfp_tfa_avg_tfr_across_sites(lfp_tfr.session(m_idx), lfp_tfa_cfg);
            end
            % LFP evoked response
            if any(strcmp(lfp_tfa_cfg.analyses, 'evoked'))
                lfp_evoked.sites_avg(m) = ...
                    lfp_tfa_avg_evoked_LFP_across_sites(lfp_evoked.session(m_idx), lfp_tfa_cfg);
            end
            % LFP power spectrum
            if any(strcmp(lfp_tfa_cfg.analyses, 'pow'))
                lfp_pow.sites_avg(m) = ...
                    lfp_tfa_avg_pow_across_sites(lfp_pow.session(m_idx), lfp_tfa_cfg);
            end
            % LFP-LFP phase sync
            if any(strcmp(lfp_tfa_cfg.analyses, 'sync'))
                lfp_tfa_avg_sitepairs_sync(sessions_info(m_idx), lfp_tfa_cfg);
            end
            % LFP-LFP phase sync spectrum
            if any(strcmp(lfp_tfa_cfg.analyses, 'syncsp'))
                lfp_tfa_avg_sitepairs_syncspctrm(sessions_info(m_idx), lfp_tfa_cfg);
            end
        end
    end
end

close all;
clear;