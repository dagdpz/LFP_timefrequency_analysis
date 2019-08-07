function session_info = lfp_tfa_process_LFP( session_info, lfp_tfa_cfg )

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
    
%     % Read in LFP data for the session - check if this is better than the
%     current approach
%     fprintf('Reading processed LFP data \n');
%     session = load(lfp_tfa_cfg.data_filepath);
    
    load(session_info.Input, 'sites');
    
    if isfield(session_info, 'Input_ECG')
        load(session_info.Input_ECG, 'out');
        if exist('out', 'var')
            block_ECG = out;
            clear out;
        end
    end
    
    % prepare results folder
    results_fldr = fullfile(session_info.proc_results_fldr);
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % struct to save data for a site
    site_lfp = struct();
    
    % structure array to store lfp data for all sites 
    % to be used for cross power spectrum calculation
    allsites_lfp = [];    
       
    % for future use
%     usable_sites_table = table;
%     if ~isempty(lfp_tfa_cfg.sites_info)
%        usable_sites_table = lfp_tfa_cfg.sites_info;
%     end
    
    % save data inside struct 
    % first loop through each site
    for i = 1:min(length(sites), lfp_tfa_cfg.maxsites)
        
        
            % for future use
            % find if this site's entry is available in usable_sites_table
    %         if isempty(usable_sites_table(strcmp(usable_sites_table.Site_ID, ...
    %                 sites(i).site_ID),:))
    %             continue;
    %         end
            fprintf('=============================================================\n');
            fprintf('Processing site, %s\n', sites(i).site_ID);
            % for future use
            % get 'Set' entry from usable_sites_table
    %         site_lfp.dataset = usable_sites_table(...
    %             strcmp(usable_sites_table.Site_ID, sites(i).site_ID), :).Set(1);
            site_lfp.site_ID = sites(i).site_ID;
            site_lfp.target = sites(i).target;
            site_lfp.recorded_hemisphere = upper(sites(i).target(end));            
            site_lfp.xpos = sites(i).grid_x;
            site_lfp.ypos = sites(i).grid_y;
            site_lfp.zpos = sites(i).electrode_depth;
            
        %% get information common to all sites for a session
        if i == 1
            site_lfp.session = sites(i).site_ID(1:12);
            site_lfp.ref_hemisphere = lfp_tfa_cfg.ref_hemisphere;
            % now loop through each trial for this site
            comp_trial = 0; % iterator for completed trials
            for t = 1:length(sites(i).trial)
                completed = sites(i).trial(t).completed;
                if true
                    type = sites(i).trial(t).type;
                    effector = sites(i).trial(t).effector;
                    run = sites(i).trial(t).run;
                    block = sites(i).trial(t).block;
                    dataset = sites(i).trial(t).dataset;
                    completed = sites(i).trial(t).completed;
                    % for future use
                    % check if the block is usable
    %                 if isempty(usable_sites_table(strcmp(usable_sites_table.Site_ID, ...
    %                         sites(i).site_ID) && usable_sites_table.Block == block))
    %                     continue;
    %                 end
                    choice_trial = sites(i).trial(t).choice;
                    reach_hand = sites(i).trial(t).reach_hand; % 1 = left, 2 = right
                    perturbation = sites(i).trial(t).perturbation; % 0 = control
                    if isnan(perturbation)
                        perturbation = 0;
                    end
                    tar_pos = sites(i).trial(t).tar_pos;
                    fix_pos = sites(i).trial(t).fix_pos;

                    % reach space         
                    if sign(real(tar_pos) - real(fix_pos)) == -1
                        reach_space = 'L'; 
                    elseif sign(real(tar_pos) - real(fix_pos)) == 1
                        reach_space = 'R';
                    else
                        reach_space = 'N';
                    end

                    % reach hand
                    if reach_hand == 1
                        reach_hand = 'L'; 
                    elseif reach_hand == 2
                        reach_hand = 'R'; 
                    else
                        reach_hand = 'N';  % no hand labeling
                    end               


                    % assign hand-space for the trial
                    if strcmp(site_lfp.ref_hemisphere, reach_space)
                         if strcmp(site_lfp.ref_hemisphere, reach_hand)
                            hs_label = 'IH IS';
                        else
                            hs_label = 'CH IS';
                        end                
                    else 
                        if strcmp(site_lfp.ref_hemisphere, reach_hand)
                            hs_label = 'IH CS';
                        else
                            hs_label = 'CH CS';
                        end
                    end

                    % check if this kind of labeling is required
    %                 if reach_hand == 'R' && reach_space == 'R'
    %                     hs_label = 'RH RS';
    %                 elseif reach_hand == 'R' && reach_space == 'L'
    %                     hs_label = 'RH LS';
    %                 elseif reach_hand == 'L' && reach_space == 'R'
    %                     hs_label = 'LH RS';
    %                 elseif reach_hand == 'L' && reach_space == 'L'
    %                     hs_label = 'LH LS';
    %                 end


                    start_time = (sites(i).trial(t).TDT_LFPx_tStart); % trial start time
                    fs = sites(i).trial(t).TDT_LFPx_SR; % sample rate
                    ts = (1/fs); % sample time
                    LFP = sites(i).trial(t).LFP; % LFP data
                    nsamples = numel(LFP);
                    end_time = start_time + (ts*(nsamples-1));
                    timestamps = linspace(start_time, end_time, nsamples);                    
                    trial_onset_time = sites(i).trial(t).trial_onset_time;

                    % save retrieved data into struct
                    comp_trial = comp_trial + 1;
                    site_lfp.trials(comp_trial).completed = completed;
                    site_lfp.trials(comp_trial).type = type;
                    site_lfp.trials(comp_trial).effector = effector;
                    site_lfp.trials(comp_trial).run = run;
                    site_lfp.trials(comp_trial).block = block;
                    site_lfp.trials(comp_trial).dataset = dataset;
                    site_lfp.trials(comp_trial).choice_trial = choice_trial;
                    site_lfp.trials(comp_trial).time = timestamps;
                    site_lfp.trials(comp_trial).lfp_data = LFP;
                    site_lfp.trials(comp_trial).fsample  = fs;
                    site_lfp.trials(comp_trial).tsample = ts;
                    site_lfp.trials(comp_trial).reach_hand  = reach_hand;
                    site_lfp.trials(comp_trial).reach_space  = reach_space;
                    site_lfp.trials(comp_trial).hndspc_lbl  = hs_label;
                    site_lfp.trials(comp_trial).perturbation  = perturbation;
                    % flag to mark noisy trials, default False, filled in by
                    % lfp_tfa_reject_noisy_lfp.m
                    site_lfp.trials(comp_trial).noisy = ~completed;

                    % get state onset times and onset samples - test and delete
                    site_lfp.trials(comp_trial).states = struct();

                        for s = 1:length(sites(i).trial(t).states)
                            % get state ID
                            state_id = sites(i).trial(t).states(s);
                            % get state onset time
                            state_onset = sites(i).trial(t).states_onset(sites(i).trial(t).states == ...
                                state_id);
                            % get sample number of state onset time
                            state_onset_sample = find(abs(timestamps - state_onset) == ...
                                min(abs(timestamps - state_onset)));
                            % save into struct
                            site_lfp.trials(comp_trial).states(s).id = state_id;
                            site_lfp.trials(comp_trial).states(s).onset_t  = state_onset;
                            site_lfp.trials(comp_trial).states(s).onset_s  = state_onset_sample;
                        end
                    %end
                    if site_lfp.trials(comp_trial).completed
                        trial_start_t = site_lfp.trials(comp_trial).states(...
                            [site_lfp.trials(comp_trial).states.id] == ...
                            lfp_tfa_cfg.trialinfo.start_state).onset_t + ...
                            lfp_tfa_cfg.trialinfo.ref_tstart;
                        trial_end_t = site_lfp.trials(comp_trial).states( ...
                            [site_lfp.trials(comp_trial).states.id] == ...
                            lfp_tfa_cfg.trialinfo.end_state).onset_t + ...
                            lfp_tfa_cfg.trialinfo.ref_tend;
                        site_lfp.trials(comp_trial).trialperiod = [trial_start_t, ...
                            trial_end_t];
                    end                        

                    site_lfp.trials(comp_trial).trial_onset_time = trial_onset_time;                    

                end
            end
            
        % Get ECG spikes
        if exist('block_ECG', 'var')
            site_lfp = lfp_tfa_get_ECG_peaks( site_lfp, block_ECG );
        end    
            
        else
            % loop through each trial for this site to get the LFP data
            comp_trial = 0; % iterator for completed trials
            for t = 1:length(sites(i).trial)
                completed = sites(i).trial(t).completed;
                if true
                    LFP = sites(i).trial(t).LFP; % LFP data
                    comp_trial = comp_trial + 1;
                    % overwrite only LFP data
                    site_lfp.trials(comp_trial).lfp_data = LFP;
                end
            end
        end    
        
        
        %%% Noise rejection - should this be included within processing check this? %%%
        %state_filt_lfp(i) = lfp_tfa_reject_noisy_lfp( state_lfp(i), lfp_tfa_cfg.noise );
        
        %% Time frequency spectrogram calculation
        site_lfp = lfp_tfa_compute_site_tfr( site_lfp, lfp_tfa_cfg );
        
        
               
        % Noise rejection
        if lfp_tfa_cfg.noise.detect
            site_lfp = lfp_tfa_reject_noisy_lfp_trials( site_lfp, lfp_tfa_cfg.noise );
        end
        
        % Baseline power calculation
        site_lfp = lfp_tfa_compute_site_baseline( site_lfp, session_info, lfp_tfa_cfg );
        
        % save data
        results_mat = fullfile(results_fldr, ['site_lfp_pow_' site_lfp.site_ID '.mat']);
        %site_lfp = state_lfp(i);
        save(results_mat, 'site_lfp', '-v7.3');
        
        % accumulate lfp for all sites
        site_lfp.trials = rmfield(site_lfp.trials, 'tfs');
        allsites_lfp = [allsites_lfp, site_lfp];
        
    end
    
    % save allsites_lfp
    results_mat = fullfile(results_fldr, 'allsites_lfp.mat');
    save(results_mat, 'allsites_lfp', '-v7.3');
    
    %% calculate cross power spectrum between sites and sync measure spectrogram
    p = 0; % sitepair index
    if any(strcmp(lfp_tfa_cfg.analyses, 'sync'))
        session_info.sitepair_info = struct();
        % prepare results folder
        results_fldr = fullfile(session_info.proc_results_fldr, 'crossspectrum');
        if ~exist(results_fldr, 'dir')
            mkdir(results_fldr);
        end
        % loop through each site
        for i = 1:length(allsites_lfp)-1
            site1_lfp = allsites_lfp(i);
            % pair a site
            for j = i+1:length(allsites_lfp)
                site2_lfp = allsites_lfp(j);
                fprintf('Computing cross power spectrum for site pair %s - %s\n', ...
                    site1_lfp.site_ID, site2_lfp.site_ID);
                sitepair_crosspow = lfp_tfa_compute_sitepair_crosspow(site1_lfp, site2_lfp, lfp_tfa_cfg);
                % save data
                results_mat = fullfile(results_fldr, ['sites_crosspow_', sitepair_crosspow.sites{1} '-' sitepair_crosspow.sites{2} '.mat']);
                save(results_mat, 'sitepair_crosspow', '-v7.3');
                
                % store info about sitepair
                p = p + 1;
                session_info.sitepair_info(p).sites = sitepair_crosspow.sites;
                session_info.sitepair_info(p).targets = sitepair_crosspow.targets;
                session_info.sitepair_info(p).session = sitepair_crosspow.session;
                session_info.sitepair_info(p).csd_results_file = results_mat;

                % compute ppc spectrogram between sitepair
                % get the trial conditions for this session
                conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg, ...
                    {session_info.Preinj_blocks, session_info.Postinj_blocks});
                sitepair_sync_file = lfp_tfa_sitepair_averaged_sync(sitepair_crosspow, conditions, lfp_tfa_cfg);
                session_info.sitepair_info(p).sync_results_file = sitepair_sync_file;
                
                % compute ppc spectrum between sitepair
                if any(strcmp(lfp_tfa_cfg.analyses, 'syncspctrm'))
                    % get the trial conditions for this session
                    conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg, ...
                        {session_info.Preinj_blocks, session_info.Postinj_blocks});
                    sitepair_syncspctrm_file = lfp_tfa_sitepair_averaged_syncspctrm(sitepair_crosspow, conditions, lfp_tfa_cfg);  
                    session_info.sitepair_info(p).syncspctrm_results_file = sitepair_syncspctrm_file;
                end               
                
            end
        end         
    end
    
    %% calculate sync measure spectrum
    if false%any(strcmp(lfp_tfa_cfg.analyses, 'syncspctrm'))
        % loop through each site
        for i = 1:length(allsites_lfp)-1
            site1_lfp = allsites_lfp(i);
            % pair a site
            for j = i+1:length(allsites_lfp)
                site2_lfp = allsites_lfp(j);
                fprintf('Computing sync spectrum for site pair %s - %s\n', ...
                    site1_lfp.site_ID, site2_lfp.site_ID);
                % compute ppc spectrum between sitepair
                % get the trial conditions for this session
                conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg, ...
                    {session_info.Preinj_blocks, session_info.Postinj_blocks});
                sitepair_syncspctrm = lfp_tfa_sitepair_avg_syncspctrum(site1_lfp, site2_lfp, conditions, lfp_tfa_cfg);                
            end
        end  
    end

end

