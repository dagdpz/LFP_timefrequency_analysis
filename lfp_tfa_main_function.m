function lfp_tfa_main_function( settings_filepath, process_LFP )
%lfp_tfa_main_function - Function to call other functions for LFP
%processing and time-frequency analysis
%
% USAGE:
%	lfp_tfa_main_function( settings )
%
% INPUTS:
%       settings                  - cell array containing the absolute path
%       to the matlab scripts containing the configurations for running the
%       time-frequency analysis
%
% OUTPUTS:
%		None
%
% REQUIRES:	lfp_tfa_define_settings, lfp_tfa_process_LFP, lfp_tfa_reject_noisy_lfp
% lfp_tfa_plot_site_average_tfr, lfp_tfa_plot_site_evoked_LFP, lfp_tfa_plot_site_powspctrum, 
% lfp_tfa_avg_tfr_across_sessions, lfp_tfa_avg_tfr_across_sites, 
% lfp_tfa_avg_evoked_LFP_across_sessions,
% lfp_tfa_avg_evoked_LFP_across_sites, 
% lfp_tfa_avg_pow_across_sessions, lfp_tfa_avg_pow_across_sites
% 
%
% See also lfp_tfa_define_settings, lfp_tfa_process_LFP, lfp_tfa_reject_noisy_lfp
% lfp_tfa_plot_site_average_tfr, lfp_tfa_plot_site_evoked_LFP, lfp_tfa_plot_site_powspctrum, 
% lfp_tfa_avg_tfr_across_sessions, lfp_tfa_avg_tfr_across_sites, 
% lfp_tfa_avg_evoked_LFP_across_sessions,
% lfp_tfa_avg_evoked_LFP_across_sites, 
% lfp_tfa_avg_pow_across_sessions, lfp_tfa_avg_pow_across_sites 
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
    

    %% INITIALIZATION
    close all;  
    
    %% loop through each settings file        

    %% INITIALIZATION
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

    try
        %% LFP processing
        if process_LFP
            % loop through each session to process
            for i = 1:length(sessions_info)
                % name of session = [Monkey name '_' Recording date]
                session_name = [sessions_info(i).Monkey '_' sessions_info(i).Date];
                fprintf('Processing LFP for session %s\n', session_name);
                lfp_tfa_cfg.session = session_name;
                % folder to which results of analysis of this session should be
                % stored
                sessions_info(i).proc_results_fldr = ...
                    fullfile(lfp_tfa_cfg.proc_lfp_folder, session_name);
                % read LFP data for each site and each trial and calculate the 
                % trial-wise time frequency spectrogram
                sessions_info(i) = ...
                    lfp_tfa_process_LFP(sessions_info(i), lfp_tfa_cfg);

            end
            lfp_tfa_cfg.sessions_info = sessions_info;
        end
        %% loop through each processed session for analysis
        for i = 1:length(sessions_info)
            % clear variables
            clear session_proc_lfp session_lfp;

            session_name = [sessions_info(i).Monkey '_' sessions_info(i).Date];
            fprintf('Analysing LFP for session %s\n', session_name);
            lfp_tfa_cfg.session = session_name;

            % folder to which results of analysis of this session should be
            % stored
            sessions_info(i).analyse_lfp_fldr = ...
                fullfile(lfp_tfa_cfg.analyse_lfp_folder, session_name);
            lfp_tfa_cfg.session_results_fldr = ...
                fullfile(lfp_tfa_cfg.analyse_lfp_folder, session_name);        

            % absolute path of file containing LFP data for this session
            % check if folder exists and there are processed LFP files inside
            if ~exist(sessions_info(i).proc_results_fldr, 'dir') || ...
                    isempty(dir(fullfile(sessions_info(i).proc_results_fldr, '*.mat')))
                warning('No mat files with processed LFP found for session %s\n'...
                    , session_name);
                continue;
            end            
            % read the processed lfp mat files for all sites of this session
            sites_lfp_files = dir(fullfile(sessions_info(i).proc_results_fldr, '*.mat'));
            session_proc_lfp = [];
            for file = {sites_lfp_files.name}
                %fprintf('Reading processed LFP for site %s\n', file{1});
                if ~strcmp(file{1}, 'allsites_lfp.mat')
                    load(fullfile(sessions_info(i).proc_results_fldr, file{1}))
                    session_proc_lfp = [session_proc_lfp site_lfp];
                end            
            end

            lfp_tfa_cfg.data_filepath = lfp_datafiles{i};
            % perturbation blocks
            lfp_tfa_cfg.perturbation_groups = {};        
            % preinjection blocks for this session
            if isfield(sessions_info, 'Preinj_blocks') && ...
                    ~isempty(sessions_info(i).Preinj_blocks)
                lfp_tfa_cfg.perturbation_groups = ...
                    [lfp_tfa_cfg.perturbation_groups sessions_info(i).Preinj_blocks];
            elseif any(lfp_tfa_cfg.compare.perturbations == 0)            
                lfp_tfa_cfg.perturbation_groups = [lfp_tfa_cfg.perturbation_groups 0];            
            end
            % postinjection blocks for this session
            if isfield(sessions_info, 'Postinj_blocks') && ...
                    ~isempty(sessions_info(i).Postinj_blocks)
                lfp_tfa_cfg.perturbation_groups{2} = sessions_info(i).Postinj_blocks;
            elseif any(lfp_tfa_cfg.compare.perturbations == 1)            
                lfp_tfa_cfg.perturbation_groups = [lfp_tfa_cfg.perturbation_groups 'all'];            
            end
            % Calculate and plot the site-wise and session average TFR, 
            % evoked response and power spectral density for 
            % different conditions and hand-space labels
            if any(strcmp(lfp_tfa_cfg.analyses, 'tfs'))
                lfp_tfr.session(i) = ...
                    lfp_tfa_plot_site_average_tfr( session_proc_lfp, ...
                    lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
            end
            if any(strcmp(lfp_tfa_cfg.analyses, 'evoked'))
                lfp_evoked.session(i) = ...
                    lfp_tfa_plot_site_evoked_LFP( session_proc_lfp, ...
                    lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );
            end
            if any(strcmp(lfp_tfa_cfg.analyses, 'pow'))
                lfp_pow.session(i) = ...
                    lfp_tfa_plot_site_powspctrum( session_proc_lfp, ...
                    lfp_tfa_cfg ) ; 
            end

            % Calculate the session-wise average of LFP-LFP phase sync
            if any(strcmp(lfp_tfa_cfg.analyses, 'sync')) && ...
                    any(strcmp(lfp_tfa_cfg.compute_avg_across, 'sessions'))
                sessions_info(i).avg_sync_results = lfp_tfa_avg_sitepairs_sync(sessions_info(i), lfp_tfa_cfg);
            end

            % Calculate the session-wise average of LFP-LFP phase sync spectrum
            if any(strcmp(lfp_tfa_cfg.analyses, 'syncspctrm')) && ...
                    any(strcmp(lfp_tfa_cfg.compute_avg_across, 'sessions'))
                sessions_info(i).avg_syncspctrm_results = lfp_tfa_avg_sitepairs_syncspctrm(sessions_info(i), lfp_tfa_cfg);
            end

        end
    catch e
        error(e.message());
    end

    %% Average across sessions
    if length(sessions_info) > 1
        % average session averages
        if any(strcmp(lfp_tfa_cfg.compute_avg_across, 'sessions'))
            % Average of session averages of LFP TFR, LFP evoked responce and
            % LFP power spectrum response
        if any(strcmp(lfp_tfa_cfg.analyses, 'tfs'))
        lfp_tfr.sessions_avg = ...
            lfp_tfa_avg_tfr_across_sessions(lfp_tfr, lfp_tfa_cfg);
        end
        if any(strcmp(lfp_tfa_cfg.analyses, 'evoked'))
        lfp_evoked.sessions_avg = ...
            lfp_tfa_avg_evoked_LFP_across_sessions(lfp_evoked, lfp_tfa_cfg);
        end
        if any(strcmp(lfp_tfa_cfg.analyses, 'pow'))
        lfp_pow.sessions_avg = ...
            lfp_tfa_avg_pow_across_sessions(lfp_pow, lfp_tfa_cfg);
        end
        % LFP-LFP phase sync
        if any(strcmp(lfp_tfa_cfg.analyses, 'sync'))
            lfp_tfa_avg_sessions_sync(sessions_info, lfp_tfa_cfg);
        end
        % LFP-LFP phase sync
        if any(strcmp(lfp_tfa_cfg.analyses, 'syncspctrm'))
            lfp_tfa_avg_sessions_syncspctrm(sessions_info, lfp_tfa_cfg);
        end
        end
        % average site averages
        if any(strcmp(lfp_tfa_cfg.compute_avg_across, 'sites'))
            % Average of site averages of LFP TFR, LFP evoked response and LFP
            % power spectrum
        % clear memory
        %clear session_proc_lfp;
            if any(strcmp(lfp_tfa_cfg.analyses, 'tfs'))
                lfp_tfr.sites_avg = ...
                    lfp_tfa_avg_tfr_across_sites(lfp_tfr, lfp_tfa_cfg);
            end
            % LFP evoked response
            if any(strcmp(lfp_tfa_cfg.analyses, 'evoked'))
                lfp_evoked.sites_avg = ...
                    lfp_tfa_avg_evoked_LFP_across_sites(lfp_evoked, lfp_tfa_cfg);
            end
            % LFP power spectrum
            if any(strcmp(lfp_tfa_cfg.analyses, 'pow'))
                lfp_pow.sites_avg = ...
                    lfp_tfa_avg_pow_across_sites(lfp_pow, lfp_tfa_cfg);
            end
            % LFP-LFP phase sync
            if any(strcmp(lfp_tfa_cfg.analyses, 'sync'))
                lfp_tfa_avg_sitepairs_sync(sessions_info, lfp_tfa_cfg);
            end
            % LFP-LFP phase sync spectrum
            if any(strcmp(lfp_tfa_cfg.analyses, 'syncspctrm'))
                lfp_tfa_avg_sitepairs_syncspctrm(sessions_info, lfp_tfa_cfg);
            end
        end
    end



end
