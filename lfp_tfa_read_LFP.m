function [ states_lfp ] = lfp_tfa_read_LFP( session_filename, all_states, maxsites, root_results_fldr )

% lfp_tfa_read_LFP - function to read in the processed lfp and
% compute the time frequency spectrogram for each trial
%
% USAGE:
%	states_lfp = lfp_tfa_create_ft_datatype( session_filename, all_states, ...
%       maxsites, root_results_fldr )
%
% INPUTS:
%		session_filename 	- filename containing the LFP data ( in the
%                               format as the processed LFP from Lukas' pipeline)
%       all_states          - structure containing information about all
%                               states, see lfp_tfa_define_states
%       maxsites            - maximum no:of sites to be analysed, set to inf to
%                               analyse all sites
%       root_results_fldr   - path to save results
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
    
    % prepare a root data structure for all sites
    states_lfp = struct();

    % Read in LFP data for the session
    fprintf('Reading processed LFP data');
    load(session_filename);
    
    % prepare results folder
    sessionName = sites(1).site_ID(1:12);
    results_fldr = fullfile(root_results_fldr, date, sessionName);
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % save data inside struct 
    % first loop through each site
    for i = 1:min(length(sites), maxsites)
        states_lfp(i).session = sites(1).site_ID(1:12);
        states_lfp(i).site_ID = sites(i).site_ID;
        states_lfp(i).target = sites(i).target;
        states_lfp(i).recorded_hemispace = sites(i).target(end);
        %states_lfp(i).lesioned_hemispace = 
        % now loop through each trial for this site
        comp_trial = 0; % iterator for completed trials
        for t = 1:length(sites(i).trial)
            completed = sites(i).trial(t).completed;
            if completed
                comp_trial = comp_trial + 1;
                block = sites(i).trial(t).block;
                choice_trial = sites(i).trial(t).choice;
                reach_hand = sites(i).trial(t).reach_hand; % 1 = left, 2 = right
                perturbation = sites(i).trial(t).perturbation; % 0 = control
                tar_pos = sites(i).trial(t).tar_pos;
                fix_pos = sites(i).trial(t).fix_pos;
                      
                % reach space         
                if sign(real(tar_pos) - real(fix_pos)) == -1
                    reach_space = 'L'; 
                else
                    reach_space = 'R';
                end
                
                % reach hand
                if reach_hand == 1, reach_hand = 'L'; else reach_hand = 'R'; end               
                
                
                % assign hand-space for the trial
                if reach_hand == 'R' && reach_space == 'R'
                    hs_label = 'RH RS';
                elseif reach_hand == 'R' && reach_space == 'L'
                    hs_label = 'RH LS';
                elseif reach_hand == 'L' && reach_space == 'R'
                    hs_label = 'LH RS';
                elseif reach_hand == 'L' && reach_space == 'L'
                    hs_label = 'LH LS';
                end
%                     if session_lfp(i).ipsilateral == reach_space
%                         hs_label = 'IH IS';
%                     else
%                         hs_label = 'IH CS';
%                     end
%                 else 
%                     if session_lfp(i).ipsilateral == reach_space
%                         hs_label = 'CH IS';
%                     else
%                         hs_label = 'CH CS';
%                     end
%                 end

                %states = trial(i).states;
                start_time = (sites(i).trial(t).TDT_LFPx_tStart); % trial start time
                fs = sites(i).trial(t).TDT_LFPx_SR; % sample rate
                LFP = sites(i).trial(t).LFP; % LFP data
                ts = (1/fs); % sample time
                nsamples = numel(LFP);
                end_time = start_time + (ts*(nsamples-1));
                timestamps = linspace(start_time, end_time, nsamples);
                
                % create an FT_DATATYPE_RAW
                ft_data_lfp = struct();
                ft_data_lfp.trial = {LFP};
                ft_data_lfp.time = {timestamps};
                ft_data_lfp.fsample = fs;
                ft_data_lfp.sampleinfo = [timestamps(1) timestamps(end)];
                ft_data_lfp.label = {strcat(states_lfp(i).site_ID, '_Trial_', num2str(t))};
                
                % calculate the TFR
                cfg              = [];
                cfg.method       = 'wavelet';                
                cfg.width        = 4;
%                 cfg.method       = 'mtmconvol';
%                 cfg.taper        = 'hanning';
                cfg.pad          = 'nextpow2';
                cfg.foi          = 2:2:120;             % analysis 2 to 100 Hz in steps of 2 Hz
                %cfg.t_ftimwin    = ones(length(cfg.foi),1).*500*ts;    % length of time window = 0.2 sec
%                 cfg.t_ftimwin    = 4./cfg.foi;                          % 4 cycles per time window
                %cfg.t_ftimwin    = ones(length(cfg.foi),1).*500*ts;    % length of time window = 0.5 sec
                cfg.toi          = timestamps(1):ts:timestamps(end);% time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
                cfg.channel      = ft_data_lfp.label;
                TFR_hann_fixed   = ft_freqanalysis(cfg, ft_data_lfp);

                % save retrieved data into struct
                states_lfp(i).trials(comp_trial).block = block;
                states_lfp(i).trials(comp_trial).choice_trial = choice_trial;
                states_lfp(i).trials(comp_trial).lfp_data = LFP;
                states_lfp(i).trials(comp_trial).time = timestamps;
                states_lfp(i).trials(comp_trial).fsample  = fs;
                states_lfp(i).trials(comp_trial).tsample = ts;
                states_lfp(i).trials(comp_trial).reach_hand  = reach_hand;
                states_lfp(i).trials(comp_trial).reach_space  = reach_space;
                states_lfp(i).trials(comp_trial).hndspc_lbl  = hs_label;
                states_lfp(i).trials(comp_trial).perturbation  = perturbation;
                states_lfp(i).trials(comp_trial).tfs = TFR_hann_fixed;

    %             % get movement onset time
    %             reach_start = sites(i).trial(t).states_onset(sites(i).trial(t).states == 62);
    %             % give cue_onset zero timestamp 
    %             reach_start_sample = find(abs(timestamps - reach_start) == min(abs(timestamps - reach_start)));

                % get state onset times and onset samples
                states_lfp(i).trials(comp_trial).states = struct();
                for s = 1:length(all_states)
                    % get state onset time
                    state_onset = sites(i).trial(t).states_onset(sites(i).trial(t).states == all_states(s).state_ID);
                    % get sample number of state onset time
                    state_onset_sample = find(abs(timestamps - state_onset) == min(abs(timestamps - state_onset)));
                    % save into struct
                    states_lfp(i).trials(comp_trial).states(s).id = all_states(s).state_ID;
                    states_lfp(i).trials(comp_trial).states(s).name  = all_states(s).state_name;
                    states_lfp(i).trials(comp_trial).states(s).onset_t  = state_onset;
                    states_lfp(i).trials(comp_trial).states(s).start_t  = ...
                        max(timestamps(1), state_onset - all_states(s).tminus_onset);
                    states_lfp(i).trials(comp_trial).states(s).end_t = ...
                        min(timestamps(end), state_onset + all_states(s).tplus_onset);
                    states_lfp(i).trials(comp_trial).states(s).onset_s  = state_onset_sample;
                    states_lfp(i).trials(comp_trial).states(s).start_s  = ...
                        max(1, state_onset_sample - round(all_states(s).tminus_onset / ts));
                    states_lfp(i).trials(comp_trial).states(s).end_s = ...
                        min(length(timestamps), state_onset_sample + ...
                        round(all_states(s).tminus_onset / ts));
                    
                end
                
                trial_start_t = min([states_lfp(i).trials(comp_trial).states.start_t]);
                trial_end_t = max([states_lfp(i).trials(comp_trial).states.end_t]);
                states_lfp(i).trials(comp_trial).trialperiod = [trial_start_t, ...
                    trial_end_t];
                
                % get baseline samples
%                 states_lfp(i).trials(comp_trial).baseline = struct();
%                 states_lfp(i).trials(comp_trial).baseline.ref_t = ...
%                     sites(i).trial(t).states_onset(sites(i).trial(t).states == baseline.ref_state);
%                 states_lfp(i).trials(comp_trial).baseline.start_t = ...
%                     states_lfp(i).trials(comp_trial).baseline.ref_t + baseline.period(1);
%                 states_lfp(i).trials(comp_trial).baseline.end_t = ...
%                     states_lfp(i).trials(comp_trial).baseline.ref_t + baseline.period(2);
%                 states_lfp(i).trials(comp_trial).baseline.start_s = ...
%                     max(1, sample(states_lfp(i).trials(comp_trial).baseline.start_t));
%                 states_lfp(i).trials(comp_trial).baseline.end_s = ...
%                     min(length(timestamps), sample(states_lfp(i).trials(comp_trial).baseline.end_t));
                
            end
        end
        %sites_lfp(i).lfp_data.fsample = fs;
    end

    % save data
    results_mat = fullfile(results_fldr, '\states_lfp.mat');
    save(results_mat, 'states_lfp');

end

