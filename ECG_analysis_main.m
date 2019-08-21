%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading processed LFP data, rejection of noise trials
% and task specific analysis using TFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear; 

% file containing settings for LFP analysis

settings_filepath = 'C:\Users\snair\Documents\GitHub\LFP_timefrequency_analysis\settings\lfp_tfa_settings_Cornelius_ecg.m';


% whether the LFP should be processed (true) or not (false)
% if the LFP for the sessions to analyse has already been processed, and no
% settings need to be changed, this flag can be set to false to skip LFP
% being processed again
% If LFP was not previously processed and the flag is set to false,
% analysis won't happen
% TODO: check if LFP is processed, of not, process LFP even if flag is set
% to false
process_LFP = true;

%% INITIALIZATION
close all;


lfp_tfa_cfg = lfp_tfa_define_settings(settings_filepath);

%% Get info about sessions to be analysed

% Read the info about sessions to analyse
sessions_info = lfp_tfa_cfg.session_info;

%% initialize structs to store intermediate results
% struct to store average ECG evoked response for different conditions
ecg_evoked = struct();
tfs_ecg = struct();

try
    % loop through each session to process
    for i = 1:length(sessions_info)
        % name of session = [Monkey name '_' Recording date]
        session_name = [sessions_info(i).Monkey '_' sessions_info(i).Date];
        fprintf('Reading ECG for session %s\n', session_name);
        lfp_tfa_cfg.session = session_name;
        sessions_info(i).sitepair_info = struct();
        % folder to which results of analysis of this session should be
        % stored
        sessions_info(i).proc_results_fldr = ...
            fullfile(lfp_tfa_cfg.proc_lfp_folder, session_name);
        % read LFP data for each site and each trial and calculate the 
        % trial-wise time frequency spectrogram
        session_ecg = ...
            lfp_tfa_process_session_ECG(sessions_info(i), lfp_tfa_cfg);

    
        fprintf('Analysing for session %s\n', session_name);

        % folder to which results of analysis of this session should be
        % stored
        sessions_info(i).analyse_lfp_fldr = ...
            fullfile(lfp_tfa_cfg.analyse_lfp_folder, session_name);
        
        % Calculate and plot the session average ECG, 
        % evoked response for different conditions 
        ecg_evoked.session(i) = lfp_tfa_plot_site_evoked_ECG( session_ecg, ...
            sessions_info(i), lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
        ecg_b2bt.session(i) = lfp_tfa_plot_session_ECG_b2bt( session_ecg, ...
            sessions_info(i), lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
        tfs_ecg.session(i) = lfp_tfa_plot_session_tfs_ECG( session_ecg, ...
            sessions_info(i), lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
        
    end
catch e
    error(e.message());
end