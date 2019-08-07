function sites_lfp = lfp_coh_process_LFP( session_lfp, lfp_tfa_cfg )

% lfp_coh_process_LFP - function to read in the processed lfp and
% compute the trial-wise power spectrum for each site and cross power
% spectrum between each pair of sites and store them as .mat files
%
% USAGE:
%	sites_lfp = lfp_coh_process_LFP( session_lfp, lfp_tfa_cfg )
%
% INPUTS:
%       session_lfp         - structure containing raw LFP data for one
%       session
%       lfp_tfa_cfg         - structure containing configurations for
%       reading LFP data and calculating spectrograms
%       Required fields: 
%           datafile_path    	- filename containing the LFP data ( in the
%                               format as the processed LFP from Lukas' pipeline)
%           proc_lfp_fldr      - path to save results
%           tfr.method          - method to be used for calculating
%                               spectrogram
%           tfr.width           - width of window in cycles for
%                               multitapering the input data
%           tfr.twin            - length of time windows in seconds
%                               to be used for spectrogram calculations
%
% OUTPUTS:
%		sites_lfp            	- structure containing trial data (raw lfp,
%                               timestamps, choice/instructed, block, 
%                               control/inactivation, states info and lfp
%                               time freq spectrogram) for successful 
%                               trials for all sites in a session
%		
%
% See also lfp_tfa_define_states
    
    close all; 
   
    sites = session_lfp;
    
    % prepare results folder
    %sessionName = sites(1).site_ID(1:12);
    results_fldr = fullfile(lfp_tfa_cfg.proc_lfp_folder);
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % get info from settings
    trial_start_state = lfp_tfa_cfg.trialinfo.start_state; 
    trial_end_state = lfp_tfa_cfg.trialinfo.end_state; 
    trial_start_offset = lfp_tfa_cfg.trialinfo.ref_tstart;
    trial_end_offset = lfp_tfa_cfg.trialinfo.ref_tend;
                    
    % struct to save data
    sites_lfp = struct();
    
    if ~isempty(lfp_tfa_cfg.sites_info)
       usable_sites_table = lfp_tfa_cfg.sites_info;
    end
    
    % initialize data inside struct 
    % information common to all trials (varying across sites)
    sites_lfp.session             = {};
    sites_lfp.site_ID             = {};
    sites_lfp.target              = {};
    sites_lfp.recorded_hemisphere = {}; sites_lfp.ref_hemisphere = {};  
    % information common to all sites (varying across trials)
    sites_lfp.type                = [];
    sites_lfp.effector            = [];
    sites_lfp.block               = [];
    sites_lfp.run                 = []; 
    sites_lfp.perturbation        = []; % 0 - pre, nonzero - post
    sites_lfp.choice              = []; % 0, 1
    sites_lfp.reach_hand          = {}; % 'R', 'L'
    sites_lfp.reach_space         = {}; % 'R', 'L'
    
    sites_lfp.states              = {}; % state id
    sites_lfp.states_onset        = {}; % state onset time
    
    % information varying across trials
    % LFP data and timestamps
    sites_lfp.trial               = {}; % lfp data for each trial
    sites_lfp.time                = {}; % timestamps
    sites_lfp.fsample             = []; % sampling frequency
    sites_lfp.trialperiod         = []; % start and end time of trial  
    sites_lfp.noisy               = []; % flag to indicate if trial is noisy
    
    
    % first loop through each site to get information varying across sites
    for i = 1:min(length(sites), lfp_tfa_cfg.maxsites)
%         if isempty(usable_sites_table(strcmp(usable_sites_table.Site_ID, ...
%                 sites(i).site_ID),:))
%             continue;
%         end
        fprintf('Processing site, %s\n', sites(i).site_ID);

%         sites_lfp(i).dataset = usable_sites_table(...
%             strcmp(usable_sites_table.Site_ID, sites(i).site_ID), :).Set(1);

        % information varying across sites (common to all trials)
        sites_lfp.session               = sites(1).site_ID(1:12);
        sites_lfp.site_ID               = [sites_lfp.site_ID, sites(i).site_ID];
        sites_lfp.target                = [sites_lfp.target, sites(i).target];
        sites_lfp.recorded_hemisphere   = [sites_lfp.recorded_hemisphere, upper(sites(i).target(end))];
        sites_lfp.ref_hemisphere        = [sites_lfp.ref_hemisphere, lfp_tfa_cfg.ref_hemisphere];
        sites_lfp.electrode_depth       = [sites.electrode_depth];
        
        %usable_sites_table(...
        %    strcmp(usable_sites_table.Site_ID, sites(i).site_ID), :).Hemisphere(1);   
        % max length of a trial in this session
        all_trials = {sites(1).trial.LFP};
        trial_maxlen = 0;
        for a = all_trials
            if size(a{:},2) > trial_maxlen
                trial_maxlen = size(a{:},2);
            end
        end
        
    end
        
            
    % concatenate trial LFP data for all sites
    concat_sites_trial = vertcat(sites.trial);
    % loop through each trial to get information varying across trials
    comp_trial = 0; % iterator for completed trials
    for t = 1:length(concat_sites_trial)
        trial = concat_sites_trial(:,t);
        completed = all([trial.completed]);
        if completed
            comp_trial = comp_trial + 1;
            % check if the block is usable
%                 if isempty(usable_sites_table(strcmp(usable_sites_table.Site_ID, ...
%                         sites(i).site_ID) && usable_sites_table.Block == block))
%                     continue;
%                 end

            sites_lfp.type          = [sites_lfp.type, trial(1).type];
            sites_lfp.effector      = [sites_lfp.effector, trial(1).effector];
            sites_lfp.run           = [sites_lfp.run, trial(1).run];
            sites_lfp.block         = [sites_lfp.block, trial(1).block];
            if ~isnan(trial(1).perturbation)
                sites_lfp.perturbation  = [sites_lfp.perturbation, trial(1).perturbation];
            else
                sites_lfp.perturbation  = [sites_lfp.perturbation, 0];
            end
            sites_lfp.choice        = [sites_lfp.choice, trial(1).choice];

            tar_pos = trial.tar_pos; % target coordinates            
            fix_pos = trial.fix_pos; % fixation coordinates
            % reach space         
            if sign(real(tar_pos) - real(fix_pos)) == -1
                reach_space = 'L'; 
            else
                reach_space = 'R';
            end
            % reach hand
            reach_hand = trial(1).reach_hand;
            if reach_hand == 1, reach_hand = 'L'; else reach_hand = 'R'; end

            sites_lfp.reach_hand    = [sites_lfp.reach_hand, reach_hand];
            sites_lfp.reach_space   = [sites_lfp.reach_space, reach_space];
            sites_lfp.states        = [sites_lfp.states, trial(1).states];
            sites_lfp.states_onset  = [sites_lfp.states_onset, trial(1).states_onset];
            % flag to mark noisy trials, default False, filled in by
            % lfp_tfa_reject_noisy_lfp.m
            sites_lfp.noisy         = [sites_lfp.noisy, 0];

            % trial LFP and timestamps, timestamps same for all sites
            start_time = (trial(1).TDT_LFPx_tStart); % trial start time
            fs = trial(1).TDT_LFPx_SR; % sample rate
            ts = (1/fs); % sample time
            LFP = vertcat(trial.LFP); % LFP data for all sites
            nsamples = size(LFP,2);
            end_time = start_time + (ts*(nsamples-1));
            timestamps = linspace(start_time, end_time, nsamples);
            
                        
            % save into struct
            sites_lfp.time          = [sites_lfp.time, timestamps];
            sites_lfp.trial         = [sites_lfp.trial, LFP];
            sites_lfp.fsample       = fs;

            % trial period
            states = sites_lfp.states{comp_trial};
            states_onset = sites_lfp.states_onset{comp_trial};
            trial_start_t = states_onset(states == trial_start_state) + trial_start_offset; 
            % Find a timestamp which is closest to trial start
            trial_start_t = timestamps(abs(timestamps - trial_start_t) == ...
                min(abs(timestamps - trial_start_t)));
            trial_end_t = states_onset(states == trial_end_state) + trial_end_offset; 
            % Find a timestamp which is closest to trial end
            trial_end_t = timestamps(abs(timestamps - trial_end_t) == ...
                min(abs(timestamps - trial_end_t)));
%                 sites_lfp(i).trials(comp_trial).trialperiod = [trial_start_t, ...
%                     trial_end_t];                
            sites_lfp.trialperiod   = [sites_lfp.trialperiod; [trial_start_t, trial_end_t]];          


        end
    end       
    
    % save data
    results_mat = fullfile(results_fldr, ['sites_lfp.mat']);
    %sites_lfp = sites_lfp(i);
    save(results_mat, 'sites_lfp', '-v7.3');
    
    %% % find the power spectrum and cross power spectrum
    sites_pow = struct();
    sites_crosspow = struct();
    sites_crosspow.crsspctrm = {};
    sites_crosspow.time = {};
        
    for t = 1:length(sites_lfp.trial)
    
        % first create an ft_datatype_raw for processing in Fieldtrip
        ft_data_lfp            = struct();
        ft_data_lfp.trial      = sites_lfp.trial(t);
        ft_data_lfp.time       = sites_lfp.time(t);
        ft_data_lfp.fsample    = sites_lfp.fsample;
        ft_data_lfp.label      = sites_lfp.site_ID;
        
        timestamps = ft_data_lfp.time{1};
        
        % Calculate power spectrum and cross power spectrum
        cfg              = [];
        cfg.method       = lfp_tfa_cfg.tfr.method;
        cfg.width        = lfp_tfa_cfg.tfr.width;
        cfg.output       = 'powandcsd';
        cfg.taper        = lfp_tfa_cfg.tfr.taper;
        cfg.pad          = 'nextpow2';
        cfg.foi          = lfp_tfa_cfg.tfr.foi;   % frequencies of interest
        cfg.t_ftimwin    = ones(length(cfg.foi),1).*0.5; %lfp_tfa_cfg.tfr.t_ftimwin;           % number of cycles per time window
        cfg.toi          = timestamps(1):lfp_tfa_cfg.tfr.timestep*(1/ft_data_lfp.fsample):...
            timestamps(end);% time period of interest
        TFR              = ft_freqanalysis(cfg, ft_data_lfp);
        
        for l = 1:length(TFR.label)
            if t == 1
                sites_pow(l).freq = TFR.freq;
                sites_pow(l).dimord = TFR.dimord;
                sites_pow(l).label = TFR.label(l);
                sites_pow(l).cfg = TFR.cfg;
                sites_pow(l).time = {};
                sites_pow(l).powspctrm = {};
            end
            sites_pow(l).powspctrm = [sites_pow(l).powspctrm, TFR.powspctrm(l,:,:)];
            sites_pow(l).time = [sites_pow(l).time, TFR.time];
        end
        
        for l = 1:length(TFR.labelcmb)
            if t == 1
                sites_crosspow(l).freq = TFR.freq;
                sites_crosspow(l).dimord = TFR.dimord;
                sites_crosspow(l).labelcmb = TFR.labelcmb(l,:);
                sites_crosspow(l).cfg = TFR.cfg;
                sites_crosspow(l).time = {};
                sites_crosspow(l).crsspctrm = {};
            end
            sites_crosspow(l).crsspctrm = [sites_crosspow(l).crsspctrm, ...
                TFR.crsspctrm(l, :, :)];
            sites_crosspow(l).time = [sites_crosspow(l).time, TFR.time];
        end                
            
    end
    
    % save results
    for l = 1:length(TFR.label)
        site_pow = sites_pow(l);
        results_filename = fullfile(results_fldr, 'powspctrm', ['sites_pow_' TFR.label{l}]);
        save(results_filename, 'site_pow');
    end
    
    for l = 1:length(TFR.labelcmb)
        sitepair_crosspow = sites_crosspow(l);
        results_filename = fullfile(results_fldr, ['sites_crosspow_' TFR.labelcmb{l,1} '_' TFR.labelcmb{l,2}]);
        save(results_filename, 'sitepair_crosspow');
    end
        
end
