function session_ecg = lfp_tfa_read_combined_ECG( session_info, plottrials  )

% lfp_tfa_process_LFP - function to read in the processed lfp and
% compute the time frequency spectrogram for each trial
%
% USAGE:
%	state_lfp = lfp_tfa_process_LFP( session_lfp, lfp_tfa_cfg )
%
% INPUTS:
%       session_lfp         - structure containing raw LFP data for one
%       session
%       lfp_tfa_cfg         - structure containing configurations for
%       reading LFP data and calculating spectrograms
%       Required fields: 
%           datafile_path    	- filename containing the LFP data ( in the
%                               format as the processed LFP from Lukas' pipeline)
%           all_states          - structure containing information about all
%                               states, see lfp_tfa_define_states
%           maxsites            - maximum no:of sites to be analysed, set to inf to
%                               analyse all sites
%           root_results_fldr   - path to save results
%           tfr.method          - method to be used for calculating
%                               spectrogram
%           tfr.width           - width of window in cycles for
%                               multitapering the input data
%           tfr.twin            - length of time windows in seconds
%                               to be used for spectrogram calculations
%
% OUTPUTS:
%		states_lfp      	- structure containing trial data (raw lfp,
%                               timestamps, choice/instructed, block, 
%                               control/inactivation, states info and lfp
%                               time freq spectrogram) for successful 
%                               trials for all sites in a session
%		
%
% See also lfp_tfa_define_states
    
    close all; 
    
    if nargin < 2 
        plottrials = 0;
    end
    
%     % Read in LFP data for the session - check if this is better than the
%     current approach
%     fprintf('Reading processed LFP data \n');
%     session = load(lfp_tfa_cfg.data_filepath);

    % struct to save data for a site
    session_ecg = struct();
        
    if ~exist(session_info.Input_ECG_combined, 'dir')
        fprintf('No file with ECG data found in the specified directory \n%s\n', ...
            session_info.Input_ECG_combined);
        return;
    end
    block_files = dir(fullfile(session_info.Input_ECG_combined, '*.mat'));
    
    % prepare results folder
    results_fldr = fullfile(session_info.proc_results_fldr);
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    if isfield(session_info, 'Input_ECG')
        if ~exist(session_info.Input_ECG, 'file')
            fprintf('No file found \n%s\n', ...
                session_info.Input_ECG_raw);
            return;
        end
        load(session_info.Input_ECG);
        if exist('out', 'var')
            block_Rpeaks = out;
            clear out;
        end        
    end
          
    % for future use
%     usable_sites_table = table;
%     if ~isempty(lfp_tfa_cfg.sites_info)
%        usable_sites_table = lfp_tfa_cfg.sites_info;
%     end
    comp_trial = 0; % iterator for completed trials 
    
    % save data inside struct 
    % first loop through each site
    for b = 1:length(block_files)
        
        % get info about site
        % for future use
            % find if this site's entry is available in usable_sites_table
    %         if isempty(usable_sites_table(strcmp(usable_sites_table.Site_ID, ...
    %                 sites(i).site_ID),:))
    %             continue;
    %         end
            nrblock_combinedFile = str2num(block_files(b).name(end-5:end-4));
            run = str2num(block_files(b).name(end-8:end-7));
            
            % check if ECG peaks exists for this block
            if length(block_Rpeaks) >= b
                block_Rpeak = block_Rpeaks(b);
                if isempty(block_Rpeak) || isempty(block_Rpeak.Rpeak_t)
                    continue;
                end
            else
                continue;
            end
            
            fprintf('=============================================================\n');
            fprintf('Reading ECG for block, %g\n', nrblock_combinedFile);
            
            % block ECG timestamps
            block_ecg_timestamps = [];
            trial_ecg_timestamps = [];           
            
            
            % for future use
            % get 'Set' entry from usable_sites_table
    %         site_lfp.dataset = usable_sites_table(...
    %             strcmp(usable_sites_table.Site_ID, sites(i).site_ID), :).Set(1);
            session_ecg.session = session_info.session;
            
            load(fullfile(session_info.Input_ECG_combined, block_files(b).name));
            
            % add first trial INI to block timestamps
            ts = 1/(trial(1).TDT_ECG1_samplingrate);
             first_trial_INI_timestamps = linspace(0, ...
                 ts*(length(First_trial_INI.ECG1)-1), length(First_trial_INI.ECG1));
             trial_ecg_timestamps = first_trial_INI_timestamps;
%              block_ecg_timestamps = [block_ecg_timestamps, ...
%                  trial_ecg_timestamps];
                                    
            %% get information common to all sites for a session
            
                % now loop through each trial for this site
                for t = 1:length(trial)
                    completed = trial(t).completed;
                    if true
                        type = trial(t).type;
                        effector = trial(t).effector;
                        completed = trial(t).completed;
                        % for future use
                        % check if the block is usable
        %                 if isempty(usable_sites_table(strcmp(usable_sites_table.Site_ID, ...
        %                         sites(i).site_ID) && usable_sites_table.Block == block))
        %                     continue;
        %                 end
                        choice_trial = trial(t).choice;
                        perturbation = nan;
                        
                        if isfield(session_info, 'Preinj_blocks') && ...
                            ~isempty(session_info.Preinj_blocks) && ...
                            ismember(nrblock_combinedFile, session_info.Preinj_blocks)
                                perturbation = 0;
                        elseif exist('ses', 'var') && ...
                                nrblock_combinedFile < ses.first_inj_block
                            perturbation = 0;
                        end
                        if isfield(session_info, 'Postinj_blocks') && ...
                            ~isempty(session_info.Postinj_blocks) && ...
                            ismember(nrblock_combinedFile, session_info.Postinj_blocks)
                                perturbation = 1;
                        elseif exist('ses', 'var') && ...
                                nrblock_combinedFile >= ses.first_inj_block
                            perturbation = 1;
                        end
                                      
                        
                        start_time = 0; % trial start time
                        fs = trial(t).TDT_ECG1_samplingrate; % sample rate
                        ts = (1/fs); % sample time
                        ECG = trial(t).TDT_ECG1; % LFP data
                        nsamples = numel(ECG);
                        end_time = start_time + (ts*(nsamples-1));
                        timestamps = linspace(start_time, end_time, nsamples);  
                        trial_ecg_timestamps = ...
                            timestamps + ts + trial_ecg_timestamps(end);                           
%                         block_ecg_timestamps = [block_ecg_timestamps, ...
%                             trial_ecg_timestamps];
                        % save retrieved data into struct
                        comp_trial = comp_trial + 1;
                        session_ecg.trials(comp_trial).completed = completed;
                        session_ecg.trials(comp_trial).type = type;
                        session_ecg.trials(comp_trial).effector = effector;
                        session_ecg.trials(comp_trial).run = run;
                        session_ecg.trials(comp_trial).block = nrblock_combinedFile;
                        session_ecg.trials(comp_trial).dataset = [];
                        session_ecg.trials(comp_trial).choice_trial = choice_trial;
                        session_ecg.trials(comp_trial).reach_hand = 0;
                        session_ecg.trials(comp_trial).reach_space = 0;
                        session_ecg.trials(comp_trial).hndspc_lbl  = [];
                        session_ecg.trials(comp_trial).time = timestamps;
                        session_ecg.trials(comp_trial).ecg_data = ECG;
                        session_ecg.trials(comp_trial).fsample  = fs;
                        session_ecg.trials(comp_trial).tsample = ts;
                        session_ecg.trials(comp_trial).tstart = start_time;
                        session_ecg.trials(comp_trial).trialperiod = ...
                            [trial_ecg_timestamps(1) trial_ecg_timestamps(end)];
                        session_ecg.trials(comp_trial).perturbation  = perturbation;
                        % flag to mark noisy trials, default False, filled in by
                        % lfp_tfa_reject_noisy_lfp.m
                        session_ecg.trials(comp_trial).noisy = ~completed;

                        % get state onset times and onset samples - test and delete
                        session_ecg.trials(comp_trial).states = struct();
                            st_idx = 0;
                            for st = 1:length(trial(t).TDT_states)
                                % get state ID
                                state_id = trial(t).TDT_states(st);
                                % get state onset time
                                state_onset = trial(t).TDT_state_onsets(trial(t).TDT_states == ...
                                    state_id);
                                % for combined file, first sample
                                % correspond to timestamp = 0
                                if state_onset < 0
                                    continue;
                                end
                                % get sample number of state onset time
                                state_onset_sample = find(abs(timestamps - state_onset(1)) == ...
                                    min(abs(timestamps - state_onset(1))), 1);
                                st_idx = st_idx + 1;
                                % save into struct
                                session_ecg.trials(comp_trial).states(st_idx).id = state_id;
                                session_ecg.trials(comp_trial).states(st_idx).onset_t  = state_onset(1);
                                session_ecg.trials(comp_trial).states(st_idx).onset_s  = state_onset_sample;
                            end
                        %end
                                   

                    end
                end
                
                session_ecg = lfp_tfa_get_block_Rpeak_times( session_ecg, block_Rpeak, nrblock_combinedFile, plottrials );
                      
        
            
        %%% Noise rejection - should this be included within processing check this? %%%
        %state_filt_lfp(i) = lfp_tfa_reject_noisy_lfp( state_lfp(i), lfp_tfa_cfg.noise );

        
    end
    
%     if isfield(session_info, 'Input_ECG')
%         load(session_info.Input_ECG, 'out');
%         if exist('out', 'var')
%             block_Rpeaks = out;
%             clear out;
%         end
%         session_ecg = lfp_tfa_get_ECG_peak_times( session_ecg, block_Rpeaks );
%     end  
    
    % Calculate time frequency spectrogram of ECG
    %session_ecg = lfp_tfa_compute_ECG_spectrogram( session_ecg, lfp_tfa_cfg );
    
    % save allsites_lfp
    results_mat = fullfile(results_fldr, ['session_ecg_' session_info.session '.mat']);
    save(results_mat, 'session_ecg', '-v7.3');
    
    
end

