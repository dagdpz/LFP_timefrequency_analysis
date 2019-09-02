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
Rpeak_state_onset = struct();
tfs_ecg = struct();

%try
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
        if isempty(fieldnames(session_ecg))
            continue;
        end
    
        fprintf('Analysing for session %s\n', session_name);

        % folder to which results of analysis of this session should be
        % stored
        sessions_info(i).analyse_lfp_fldr = ...
            fullfile(lfp_tfa_cfg.analyse_lfp_folder, session_name);
        
        % Calculate and plot the session average ECG, 
        % evoked response for different conditions 
        ecg_evoked.session(i) = lfp_tfa_plot_site_evoked_ECG( session_ecg, ...
            sessions_info(i), lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
        Rpeak_state_onset.session(i) = lfp_tfa_compute_Rpeak_evoked_state_onsets...
            ( session_ecg, sessions_info(i), lfp_tfa_cfg.analyse_Rpeak_states, ...
            lfp_tfa_cfg );
%         ecg_b2bt.session(i) = lfp_tfa_plot_session_ECG_b2bt( session_ecg, ...
%             sessions_info(i), lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
%         tfs_ecg.session(i) = lfp_tfa_plot_session_tfs_ECG( session_ecg, ...
%             sessions_info(i), lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
        
    end
%catch e
%    error(e.message());
%end

%% average across sessions
if length(sessions_info) > 1
    % Average task evoked ECG
    ecg_evoked.sessions_avg = lfp_tfa_avg_sessions_ECG_evoked(ecg_evoked, ...
        lfp_tfa_cfg);
    % Average task evoked ECG b2bt
%     ecg_b2bt.sessions_avg = lfp_tfa_avg_sessions_ECG_b2bt(ecg_b2bt, ...
%         lfp_tfa_cfg);
%     % Average task evoked ECG time frequency spectrogram
%     tfs_ecg.sessions_avg = lfp_tfa_avg_sessions_ECG_tfs(tfs_ecg, ...
%         lfp_tfa_cfg);
end