function lfp_tfa_process_LFP( data_filepath, lfp_tfa_cfg )

% lfp_tfa_read_LFP - function to read in the processed lfp and
% compute the time frequency spectrogram for each trial
%
% USAGE:
%	states_lfp = lfp_tfa_create_ft_datatype( session_filename, all_states, ...
%       maxsites, root_results_fldr )
%
% INPUTS:
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
    
    % Read in LFP data for the session
    fprintf('Reading processed LFP data \n');
    session = load(data_filepath);
    
    sites = session.sites;
    
    % prepare results folder
    %sessionName = sites(1).site_ID(1:12);
    results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr);
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % struct to save data
    state_lfp = struct();
    
    % save data inside struct 
    % first loop through each site
    for i = 1:min(length(sites), lfp_tfa_cfg.maxsites)
        fprintf('Processing site, %s\n', sites(i).site_ID);
        %site_lfp = struct();
        state_lfp(i).session = sites(1).site_ID(1:12);
        state_lfp(i).site_ID = sites(i).site_ID;
        state_lfp(i).target = sites(i).target;
        state_lfp(i).recorded_hemispace = sites(i).target(end);
        state_lfp(i).lesioned_hemispace = lfp_tfa_cfg.lesional_hemispace;
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
                if strcmp(state_lfp(i).lesioned_hemispace, reach_space)
                    if strcmp(state_lfp(i).lesioned_hemispace, reach_hand)
                        hs_label = 'IH IS';
                    else
                        hs_label = 'CH IS';
                    end
                else 
                    if strcmp(state_lfp(i).lesioned_hemispace, reach_hand)
                        hs_label = 'IH CS';
                    else
                        hs_label = 'CH CS';
                    end
                end
%                 if reach_hand == 'R' && reach_space == 'R'
%                     hs_label = 'RH RS';
%                 elseif reach_hand == 'R' && reach_space == 'L'
%                     hs_label = 'RH LS';
%                 elseif reach_hand == 'L' && reach_space == 'R'
%                     hs_label = 'LH RS';
%                 elseif reach_hand == 'L' && reach_space == 'L'
%                     hs_label = 'LH LS';
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
                ft_data_lfp.label = {strcat(state_lfp(i).site_ID, '_Trial_', num2str(t))};
                
                % calculate the TFR
                cfg              = [];
                cfg.method       = lfp_tfa_cfg.tfr.method;                
                cfg.width        = lfp_tfa_cfg.tfr.width;
%                 cfg.method       = 'mtmconvol';
                cfg.taper        = lfp_tfa_cfg.tfr.taper;
                cfg.pad          = 'nextpow2';
                cfg.foi          = lfp_tfa_cfg.tfr.foi;   % analysis 2 to 100 Hz in steps of 2 Hz
                %cfg.foi          = 2:2:120;
                %cfg.t_ftimwin    = ones(length(cfg.foi),1).*500*ts;    % length of time window = 0.2 sec
                cfg.t_ftimwin    = lfp_tfa_cfg.tfr.t_ftimwin;           % 4 cycles per time window
                %cfg.t_ftimwin    = ones(length(cfg.foi),1).*500*ts;    % length of time window = 0.5 sec
                cfg.toi          = timestamps(1):lfp_tfa_cfg.tfr.timestep*ts:...
                    timestamps(end);% time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
                cfg.channel      = ft_data_lfp.label;
                TFR_wavelet      = ft_freqanalysis(cfg, ft_data_lfp);

                % save retrieved data into struct
                state_lfp(i).trials(comp_trial).block = block;
                state_lfp(i).trials(comp_trial).choice_trial = choice_trial;
                state_lfp(i).trials(comp_trial).lfp_data = LFP;
                state_lfp(i).trials(comp_trial).time = timestamps;
                state_lfp(i).trials(comp_trial).fsample  = fs;
                state_lfp(i).trials(comp_trial).tsample = ts;
                state_lfp(i).trials(comp_trial).reach_hand  = reach_hand;
                state_lfp(i).trials(comp_trial).reach_space  = reach_space;
                state_lfp(i).trials(comp_trial).hndspc_lbl  = hs_label;
                state_lfp(i).trials(comp_trial).perturbation  = perturbation;
                state_lfp(i).trials(comp_trial).tfs = TFR_wavelet;

    
                % get state onset times and onset samples
                state_lfp(i).trials(comp_trial).states = struct();
                if i > 1 && strcmp(state_lfp(i).session,  state_lfp(i-1).session)
                    state_lfp(i).trials(comp_trial).states = ...
                        state_lfp(i-1).trials(comp_trial).states;
                else
                    for s = 1:length(lfp_tfa_cfg.all_states)
                        % get state onset time
                        state_onset = sites(i).trial(t).states_onset(sites(i).trial(t).states == ...
                            lfp_tfa_cfg.all_states(s).state_ID);
                        % get sample number of state onset time
                        state_onset_sample = find(abs(timestamps - state_onset) == ...
                            min(abs(timestamps - state_onset)));
                        % save into struct
                        state_lfp(i).trials(comp_trial).states(s).id = ...
                            lfp_tfa_cfg.all_states(s).state_ID;
                        state_lfp(i).trials(comp_trial).states(s).name  = ...
                            lfp_tfa_cfg.all_states(s).state_name;
                        state_lfp(i).trials(comp_trial).states(s).onset_t  = state_onset;
                        state_lfp(i).trials(comp_trial).states(s).start_t  = ...
                            max(timestamps(1), state_onset - ...
                            lfp_tfa_cfg.all_states(s).tminus_onset);
                        state_lfp(i).trials(comp_trial).states(s).end_t = ...
                            min(timestamps(end), state_onset + ...
                            lfp_tfa_cfg.all_states(s).tplus_onset);
                        state_lfp(i).trials(comp_trial).states(s).onset_s  = state_onset_sample;
                        state_lfp(i).trials(comp_trial).states(s).start_s  = ...
                            max(1, state_onset_sample - ...
                            round(lfp_tfa_cfg.all_states(s).tminus_onset / ts));
                        state_lfp(i).trials(comp_trial).states(s).end_s = ...
                            min(length(timestamps), state_onset_sample + ...
                            round(lfp_tfa_cfg.all_states(s).tminus_onset / ts));
                    end
                end
                
                % get state onset times and onset samples
                state_lfp(i).trials(comp_trial).epochs = struct();
                if i > 1 && strcmp(state_lfp(i).session,  state_lfp(i-1).session)
                    state_lfp(i).trials(comp_trial).epochs = ...
                        state_lfp(i-1).trials(comp_trial).epochs;
                else
                    for e = 1:length(lfp_tfa_cfg.epochs)
                        % get state onset time
                        epoch_onset = sites(i).trial(t).states_onset(sites(i).trial(t).states == ...
                            lfp_tfa_cfg.epochs(e).ref_state_ID);
                        % save into struct
                        state_lfp(i).trials(comp_trial).epochs(e).name  = ...
                            lfp_tfa_cfg.epochs(e).name;
                        state_lfp(i).trials(comp_trial).epochs(e).ref_state_ID = ...
                            lfp_tfa_cfg.epochs(e).ref_state_ID;
                        state_lfp(i).trials(comp_trial).epochs(e).onset_t  = epoch_onset;
                        state_lfp(i).trials(comp_trial).epochs(e).start_t  = ...
                            max(timestamps(1), epoch_onset + lfp_tfa_cfg.epochs(e).ref_period(1));
                        state_lfp(i).trials(comp_trial).epochs(e).end_t = ...
                            min(timestamps(end), epoch_onset + lfp_tfa_cfg.epochs(e).ref_period(2));
    %                     site_lfp.trials(comp_trial).states(s).onset_s  = state_onset_sample;
    %                     site_lfp.trials(comp_trial).states(s).start_s  = ...
    %                         max(1, state_onset_sample - ...
    %                         round(lfp_tfa_cfg.all_states(s).tminus_onset / ts));
    %                     site_lfp.trials(comp_trial).states(s).end_s = ...
    %                         min(length(timestamps), state_onset_sample + ...
    %                         round(lfp_tfa_cfg.all_states(s).tminus_onset / ts));

                    end
                end
                
                trial_start_t = min([state_lfp(i).trials(comp_trial).states.start_t]);
                trial_end_t = max([state_lfp(i).trials(comp_trial).states.end_t]);
                state_lfp(i).trials(comp_trial).trialperiod = [trial_start_t, ...
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
        
        % call function for rejeting noisy trials for this site
        state_lfp(i) = lfp_tfa_reject_noisy_lfp( site_lfp, cfg_noise );
        
        % call function for calculating the baseline for this site
                
        
        % save data
        results_mat = fullfile(results_fldr, [state_lfp(i).site_ID '.mat']);
        site_lfp = state_lfp(i);
        save(results_mat, 'site_lfp', '-v7.3');
    end  

end

