%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading processed LFP data, rejection of noise trials
% and task specific analysis using TFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clear; 

% file containing settings for LFP analysis

settings_filepaths = {'C:\Users\snair\Documents\GitHub\LFP_timefrequency_analysis\settings\lfp_tfa_settings_Magnus_dPul_ECG_LFP.m'};
%    {'C:\Users\snair\Documents\GitHub\LFP_timefrequency_analysis\settings\lfp_tfa_settings_Curius_inactivation.m'};%, ...
%     'C:\Users\snair\Documents\GitHub\LFP_timefrequency_analysis\settings\lfp_tfa_settings_Curius_baseline1.m', ...
%     'C:\Users\snair\Documents\GitHub\LFP_timefrequency_analysis\settings\lfp_tfa_settings_Cornelius_inactivation_ecg_by_block.m', ...
%     'C:\Users\snair\Documents\GitHub\LFP_timefrequency_analysis\settings\lfp_tfa_settings_Cornelius_baseline_ecg_by_block.m'};


% whether the LFP should be processed (true) or not (false)
% if the LFP for the sessions to analyse has already been processed, and no
% settings need to be changed, this flag can be set to false to skip LFP
% being processed again
% If LFP was not previously processed and the flag is set to false,
% analysis won't happen
% TODO: check if LFP is processed, of not, process LFP even if flag is set
% to false
process_ECG = 1;

%% INITIALIZATION
close all;

% loop through settings file
for s = 1:length(settings_filepaths)
    settings_filepath = settings_filepaths{s};

    lfp_tfa_cfg = lfp_tfa_define_settings(settings_filepath);

    %% Get info about sessions to be analysed

    % Read the info about sessions to analyse
    sessions_info = lfp_tfa_cfg.session_info;

    %% initialize structs to store intermediate results
    % struct to store average ECG evoked response for different conditions
    Rpeak_evoked_ECG = struct();
    Rpeak_evoked_event_prob = struct();
    ECG_b2bt_evoked = struct();
    % struct to store average LFP TFR for different conditions
    Rpeak_evoked_lfp_tfr = struct();
    % struct to store average LFP evoked response for different conditions
    Rpeak_evoked_lfp_raw = struct();

    %try
        % loop through each session to process
        for i = 1:length(sessions_info)
            % name of session = [Monkey name '_' Recording date]
            session_name = [sessions_info(i).Monkey '_' sessions_info(i).Date];

            if process_ECG
                fprintf('Reading ECG for session %s\n', session_name);
                lfp_tfa_cfg.session = session_name;
                sessions_info(i).sitepair_info = struct();
                % folder to which results of analysis of this session should be
                % stored
                sessions_info(i).proc_results_fldr = ...
                    fullfile(lfp_tfa_cfg.proc_lfp_folder, session_name);
                % read LFP data for each site and each trial 
                if isfield(sessions_info(i), 'Input_ECG_combined') && ...
                        ~isempty(sessions_info(i).Input_ECG_combined)
                    session_ecg = ...
                        lfp_tfa_read_combined_ECG(sessions_info(i), lfp_tfa_cfg.plottrials);
                elseif isfield(sessions_info(i), 'Input_ECG_preproc') && ...
                        ~isempty(sessions_info(i).Input_ECG_preproc)
                    session_ecg = ...
                        lfp_tfa_read_preproc_ECG(sessions_info(i), lfp_tfa_cfg.plottrials);
                end
                
                % Read LFP data
                if isfield(sessions_info(i), 'Input_LFP') && ...
                        ~isempty(sessions_info(i).Input_LFP)
                    sessions_info(i) = ...
                        lfp_tfa_process_combined_LFP(sessions_info(i), lfp_tfa_cfg);
                end

                if isempty(fieldnames(session_ecg))
                    continue;
                end
            else
                % load session ecg for the session
                session_ecg_filename = fullfile(sessions_info(i).proc_results_fldr, ...
                    ['session_ecg_' session_name '.mat']);
                if exist(session_ecg_filename, 'file')
                    load(session_ecg_filename, 'session_ecg');
                else
                    continue;
                end
            end

            fprintf('Analysing for session %s\n', session_name);

            % folder to which results of analysis of this session should be
            % stored
            sessions_info(i).analyse_lfp_fldr = ...
                fullfile(lfp_tfa_cfg.analyse_lfp_folder, session_name);
            lfp_tfa_cfg.session_results_fldr = ...
            fullfile(lfp_tfa_cfg.analyse_lfp_folder, session_name);  
        
            % Calculate and plot the session average ECG, 
            % evoked response for different conditions 
            
            Rpeak_evoked_ECG.session(i) = lfp_tfa_compute_session_evoked_ECG( session_ecg, ...
                sessions_info(i), lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
            
            Rpeak_evoked_event_prob.session(i) = lfp_tfa_compute_session_Rpeak_evoked_state_onsets...
                ( session_ecg, sessions_info(i), lfp_tfa_cfg.analyse_Rpeak_states, ...
                lfp_tfa_cfg );
            
            ECG_b2bt_evoked.session(i) = lfp_tfa_compute_session_evoked_ECG_b2bt( session_ecg, ...
                    sessions_info(i), lfp_tfa_cfg.event_triggers, lfp_tfa_cfg );
            
            %clear session_ecg;
            
            if any(strcmp(lfp_tfa_cfg.analyses, 'Rpeak_evoked_LFP')) || ...
                    any(strcmp(lfp_tfa_cfg.analyses, 'Rpeak_evoked_TFS'))
                % read the processed lfp mat files for all sites of this session
                sites_lfp_files = dir(fullfile(sessions_info(i).proc_results_fldr, 'site_lfp_*.mat'));
                session_proc_lfp = [];
                for file = {sites_lfp_files.name}
                    %fprintf('Reading processed LFP for site %s\n', file{1});
                    if ~isempty(strfind(file{1}, 'site_lfp_'))
                        fprintf('Reading processed LFP for site %s\n', file{:});
                        load(fullfile(sessions_info(i).proc_results_fldr, file{1}))
                        session_proc_lfp = [session_proc_lfp site_lfp];
                    end            
                end
            end
            
            if any(strcmp(lfp_tfa_cfg.analyses, 'Rpeak_evoked_LFP'))
                Rpeak_evoked_lfp_raw.session(i) = ...
                    lfp_tfa_compute_session_Rpeak_evoked_LFP( session_proc_lfp, ...
                    lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
            end
            if any(strcmp(lfp_tfa_cfg.analyses, 'Rpeak_evoked_TFS'))
                Rpeak_evoked_lfp_tfr.session(i) = ...
                    lfp_tfa_compute_session_Rpeak_evoked_TFS( session_proc_lfp, ...
                    lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
            end
            
    %         tfs_ecg.session(i) = lfp_tfa_plot_session_tfs_ECG( session_ecg, ...
    %             sessions_info(i), lfp_tfa_cfg.event_triggers, lfp_tfa_cfg );       

        end
    %catch e
    %    error(e.message());
    %end

    %% average across sessions
    if length(sessions_info) > 1
        % Average task evoked ECG
        Rpeak_evoked_ECG.sessions_avg = lfp_tfa_avg_sessions_ECG_evoked(Rpeak_evoked_ECG, ...
            lfp_tfa_cfg);
        % Average task evoked ECG b2bt
        ECG_b2bt_evoked.sessions_avg = lfp_tfa_avg_sessions_ECGb2bt_evoked(ECG_b2bt_evoked, ...
            lfp_tfa_cfg);
    % %     % Average task evoked ECG time frequency spectrogram
    %     tfs_ecg.sessions_avg = lfp_tfa_avg_sessions_ECG_tfs(tfs_ecg, ...
    %         lfp_tfa_cfg);
        % Average Rpeak evoked state onset probability
        Rpeak_evoked_event_prob.sessions_avg = lfp_tfa_avg_sessions_Rpeak_evoked_state_onsets( ...
            Rpeak_evoked_event_prob, lfp_tfa_cfg);
        
        if any(strcmp(lfp_tfa_cfg.compute_avg_across, 'sessions'))
            if any(strcmp(lfp_tfa_cfg.analyses, 'Rpeak_evoked_TFS'))
                Rpeak_evoked_lfp_tfr.sessions_avg = ...
                    lfp_tfa_avg_sessions_Rpeak_evoked_TFS(Rpeak_evoked_lfp_tfr, lfp_tfa_cfg);
            end
            if any(strcmp(lfp_tfa_cfg.analyses, 'Rpeak_evoked_LFP'))
                Rpeak_evoked_lfp_raw.sessions_avg = ...
                    lfp_tfa_avg_sessions_Rpeak_evoked_LFP(Rpeak_evoked_lfp_raw, lfp_tfa_cfg);
            end
        end
        
        if any(strcmp(lfp_tfa_cfg.compute_avg_across, 'sites'))
            if any(strcmp(lfp_tfa_cfg.analyses, 'Rpeak_evoked_TFS'))
                Rpeak_evoked_lfp_tfr.sessions_avg = ...
                    lfp_tfa_avg_sites_Rpeak_evoked_TFS(Rpeak_evoked_lfp_tfr, lfp_tfa_cfg);
            end
            if any(strcmp(lfp_tfa_cfg.analyses, 'Rpeak_evoked_LFP'))
                Rpeak_evoked_lfp_raw.sessions_avg = ...
                    lfp_tfa_avg_sites_Rpeak_evoked_LFP(Rpeak_evoked_lfp_raw, lfp_tfa_cfg);
            end
        end
        
    end
    
end