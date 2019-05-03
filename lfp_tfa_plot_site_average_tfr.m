function [session_tfs] = lfp_tfa_plot_site_average_tfr( states_lfp, analyse_states, lfp_tfa_cfg ) 

% lfp_tfa_compute_plot_tfr  - plots lfp time freq response for
% different hand-space tuning conditions for each site and across sites
%
% USAGE:
%	[ session_tfs ] = lfp_tfa_plot_average_tfr( sites_lfp_folder, analyse_states, lfp_tfa_cfg )
%
% INPUTS:
%		states_lfp  	- folder containing lfp data, output from
%		lfp_tfa_compute_baseline
%       analyse_states  - cell array containing ids of states to be
%       analysed
%       lfp_tfa_cfg     - struct containing configuration for TFR 
%           Required fields:
%           trial_condition.blocks              - blocks to be analysed, 
%           leave empty to analyse each block separately
%           trial_condition.baseline_method     - method to be used for 
%           baseline normalization ('zscore', 'subtraction', 'division',
%           'relchange')
%           results_folder                      - folder to save results
% OUTPUTS:
%		session_tfs    	- output structure which saves the average tfs for  
%                         trials of a given condition for different handspace 
%                         tunings and periods around the states analysed
%                         same datastructure as input ft_data_sites, but
%                         with additional field to store condition-wise lfp
%                         time freq response		
%
%
%
% See also lfp_tfa_compute_baseline, lfp_tfa_define_states
    
    % suppress warning for xticklabel
    warning ('off', 'MATLAB:hg:willberemoved');

    % make a folder to save figures
    %sessionName = states_lfp(1).session;
    results_folder_tfr = fullfile(lfp_tfa_cfg.session_results_fldr, 'Condition_based_TFS');
    if ~exist(results_folder_tfr, 'dir')
        mkdir(results_folder_tfr);
    end
        
    recorded_hemisphere = unique([states_lfp.recorded_hemisphere]);
    choice = unique([states_lfp(1).trials.choice_trial]);
    perturbation = unique([states_lfp(1).trials.perturbation]);
    blocks = unique([states_lfp(1).trials.block]);
    
    % create conditions
    site_conditions = struct();
       
    i = 0;
    for rec_hem = recorded_hemisphere        
        for c = choice
            for b = blocks
                i = i + 1;
                site_conditions(i).recorded_hemisphere = rec_hem;
                site_conditions(i).choice = c;
                site_conditions(i).block = b;
                site_conditions(i).perturbation = ...
                    unique([states_lfp(1).trials([states_lfp(1).trials.block] == b).perturbation]);
                cond_label = [];
                if site_conditions(i).recorded_hemisphere == 'L'
                    cond_label = [cond_label 'Left_hemisphere_'];
                else
                    cond_label = [cond_label 'Right_hemisphere_'];
                end
                if site_conditions(i).choice == 0
                    cond_label = [cond_label 'Instructed_'];
                else
                    cond_label = [cond_label 'Choice_'];
                end
                if site_conditions(i).perturbation == 0
                    cond_label = [cond_label 'Control_'];
                else
                    cond_label = [cond_label 'Inactivation_'];
                end
                cond_label = [cond_label, 'Block_', num2str(site_conditions(i).block)];
                site_conditions(i).label = cond_label;
                
                % create a folder for storing results for this condition
%                 cfg_conditions(i).results_folder = fullfile(results_folder_tfr); %, cfg_conditions(i).label
%                 if ~exist(cfg_conditions(i).results_folder, 'dir')
%                     mkdir(cfg_conditions(i).results_folder)
%                 end
                                
            end
        end
    end
    
    
    if isfield(lfp_tfa_cfg, 'add_conditions')
        for c = 1:length(lfp_tfa_cfg.add_conditions)
            if ~isempty(lfp_tfa_cfg.add_conditions(c))
                if ~isempty(lfp_tfa_cfg.add_conditions(c).blocks)
                    for rec_hem = recorded_hemisphere        
                        for ch = choice
                            i = i + 1;
                            if strcmp(lfp_tfa_cfg.add_conditions(c).blocks, 'inactivation')
                                site_conditions(i).block = blocks(perturbation ~= 0);
                                site_conditions(i).perturbation = 1;
                            else
                                site_conditions(i).block = lfp_tfa_cfg.add_conditions(c).blocks;
                                site_conditions(i).perturbation = perturbation(blocks == lfp_tfa_cfg.add_conditions(c).blocks(1));
                                
                            end                    
                            site_conditions(i).choice = ch;
                            if isfield(lfp_tfa_cfg.add_conditions(c), 'perturbation')
                                site_conditions(i).perturbation = lfp_tfa_cfg.add_conditions(c).perturbation;
                            end
                            site_conditions(i).recorded_hemisphere = rec_hem;
                            cond_label = [];
                            if site_conditions(i).recorded_hemisphere == 'L'
                                cond_label = [cond_label 'Left_hemispace_'];
                            else
                                cond_label = [cond_label 'Right_hemispace_'];
                            end
                            if site_conditions(i).choice == 0
                                cond_label = [cond_label 'Instructed_'];
                            else
                                cond_label = [cond_label 'Choice_'];
                            end
                            if site_conditions(i).perturbation == 0
                                cond_label = [cond_label 'Control_'];
                            else
                                cond_label = [cond_label 'Inactivation_'];
                            end
                            cond_label = [cond_label, 'Block_', num2str(site_conditions(i).block)];
                            site_conditions(i).label = cond_label;
                            
                            % create a folder for storing results for this condition
%                             cfg_conditions(i).results_folder = fullfile(results_folder_tfr); %, cfg_conditions(i).label
%                             if ~exist(cfg_conditions(i).results_folder, 'dir')
%                                 mkdir(cfg_conditions(i).results_folder)
%                             end
                        end
                    end
                end
            end
        end
    end
    
    %lfp_tfa_cfg.conditions = cfg_conditions;
    
    %usable_sites_table = lfp_tfa_cfg.sites_info;
       
    % condition based TFS
    sites_tfr = struct();
    session_tfs = struct();
    session_tfs.session = states_lfp(1).session;
    for i = 1:length(states_lfp)
        % check if site is usable
%         if isempty(usable_sites_table(strcmp(usable_sites_table.Site_ID, ...
%                 states_lfp(i).site_ID)))
%             continue;
%         end
        
        % target for this site
        target = unique({states_lfp(i).target});
        site_conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg, target);
        
        % folder to save sitewise results
        site_results_folder = fullfile(results_folder_tfr, states_lfp(i).site_ID);
        if ~exist(site_results_folder, 'dir')
            mkdir(site_results_folder);
        end
        sites_tfr(i).condition = struct();
        sites_tfr(i).site_ID = states_lfp(i).site_ID;
        sites_tfr(i).session = states_lfp(i).session;
        sites_tfr(i).target = states_lfp(i).target;
        % flag to indicate if this site should be used for
        % averaging based on minimum no:of trials per condition
        sites_tfr(i).use_for_avg = 1;
        for cn = 1:length(site_conditions)

            % hand-space tuning of LFP
            hs_labels = site_conditions(cn).hs_labels;

            % current trial condition analysed
            cfg_condition = site_conditions(cn);


            % cell array to store time frequency average across sites
            TFR_avg = cell(length(analyse_states),length(hs_labels));

            % loop through each site
            nsites = length(states_lfp);

            %cond_based_tfs(cn).sites = struct();

        % loop through each site
        %for cn = 1:length(cfg_conditions)            
            
            %cond_based_tfs(cn).sites(i) = struct();
            % consider site based on recorded hemispace
            sites_tfr(i).condition(cn).label = site_conditions(cn).label;
            sites_tfr(i).condition(cn).cfg_condition = site_conditions(cn);
            %if strcmp(states_lfp(i).target, cfg_conditions(cn).target)             
            % make a struct for concatenating TFR for all states
            sites_tfr(i).condition(cn).hs_tuned_tfs = struct(); 
            %cond_based_tfs(i).tfs_avg_site.powspctrm = cell(length(analyse_states),length(hs_labels));
            sites_tfr(i).condition(cn).ntrials = zeros(1,length(hs_labels));                

            for hs = 1:length(hs_labels)

                cond_trials = lfp_tfa_get_condition_trials(states_lfp(i), site_conditions(cn));

                cond_trials = cond_trials & ...
                    strcmp({states_lfp(i).trials.hndspc_lbl}, hs_labels(hs));
                sites_tfr(i).condition(cn).ntrials(hs) = sum(cond_trials);

                fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
                fprintf('Total number of trials %g\n', sum(cond_trials));

                sites_tfr(i).condition(cn).noisytrials(hs) = ...
                    sum(cond_trials & [states_lfp(i).trials.noisy]); 

                % consider only non noisy trials
                fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                    & [states_lfp(i).trials.noisy]));
                cond_trials = cond_trials & ~[states_lfp(i).trials.noisy];

                % check if the site contains a specified minimum number
                % of trials for all conditions
                if sum(cond_trials) < lfp_tfa_cfg.mintrials_percondition
                    sites_tfr(i).use_for_avg = 0;
                end
                % loop through trials 

                for st = 1:size(analyse_states, 1)
                    state_id = analyse_states{st, 1};
                    state_name = analyse_states{st, 2};
                    state_ref_tstart = analyse_states{st, 3};
                    state_ref_tend = analyse_states{st, 4};
                    
                    state_tfs.powspctrm = {}; % power spectrogram
                    state_tfs.time = {}; % timebins fo spectrogram
                    state_tfs.freq = {}; % freq bins

                    for t = find(cond_trials)

                        states          = states_lfp(i).trials(t).states;
                        state_onset_t   = states([states(:).id] == ...
                            state_id).onset_t;
                        state_start_t   = states([states(:).id] == ...
                            state_id).onset_t + state_ref_tstart;
                        state_end_t     = states([states(:).id] == ...
                            state_id).onset_t + state_ref_tend;
                        % sampling frequency
                        fs = states_lfp(i).trials(t).fsample;

                        % crop the tfs for this state
                        state_tfs.powspctrm = [state_tfs.powspctrm, ...
                            states_lfp(i).trials(t).tfs.powspctrm(1, :, ...
                            (states_lfp(i).trials(t).tfs.time >= state_start_t & ...
                            states_lfp(i).trials(t).tfs.time <= state_end_t))];
    %                     if sum (states_lfp(i).trials(t).tfs.time >= state_start_t & ...
    %                         states_lfp(i).trials(t).tfs.time <= state_end_t) < ntimebins
                        state_tfs.time = [state_tfs.time, ...
                            states_lfp(i).trials(t).tfs.time(1, ...
                            (states_lfp(i).trials(t).tfs.time >= state_start_t & ...
                            states_lfp(i).trials(t).tfs.time <= state_end_t)) - state_onset_t];
                        state_tfs.freq = [state_tfs.freq states_lfp(i).trials(t).tfs.freq]; 
                        state_tfs.cfg = states_lfp(i).trials(t).tfs.cfg;
    %                     end

                    end

                    % find number of time bins in power
                    % spectrogram
                    ntimebins = min(cellfun('length', state_tfs.time));
                    nfreqbins = numel(state_tfs.freq);
                    % crop each tfs to the ntimebins
                    for k = 1:length(state_tfs.powspctrm)
                        state_tfs.powspctrm{k} = state_tfs.powspctrm{k}(1,:,1:ntimebins);
                        state_tfs.time{k} = state_tfs.time{k}(1:ntimebins);
                    end

                    % average power spectrum for each state
                    arr_state_pow = zeros(1, nfreqbins, ntimebins);

                    if ~isempty(state_tfs.powspctrm)

                        % find the average TFS for each state
                        arr_state_pow = cat(1, state_tfs.powspctrm{:});
                        state_tfs.powspctrm_rawmean = nanmean(arr_state_pow, 1);

                        % baseline normalization
                        cfg_baseline.method = lfp_tfa_cfg.baseline_method;
                        cfg_baseline.mean = states_lfp(i).baseline_mean;
                        cfg_baseline.std = states_lfp(i).baseline_std;
                        state_tfs.powspctrm_normmean = lfp_tfa_baseline_normalization(...
                            state_tfs.powspctrm_rawmean, cfg_baseline); 

                        % save average tfs
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm = state_tfs.powspctrm_normmean;
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).time = state_tfs.time{1};
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq = state_tfs.freq{1}; 
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).cfg = state_tfs.cfg;
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).hs_label = hs_labels(hs);
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).state = state_id;
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).state_name = state_name;
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).trials = find(cond_trials);
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).ntrials = length(find(cond_trials));
                    end

                end

            end
%         else
%             continue;                
%         end
            
            % TFR
            if ~isempty(sites_tfr(i).condition(cn).hs_tuned_tfs)
                if site_conditions(cn).perturbation == 0
                    injection = 'Pre';
                else
                    injection = 'Post';
                end
                plottitle = ['LFP TFR (' injection '): Site ' sites_tfr(i).site_ID ', Target ' sites_tfr(i).target ', '  ...
                '(Perturbation ' num2str(site_conditions(cn).perturbation_group{1}) ')'];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                else
                    plottitle = [plottitle 'Choice trials'];
                end
                
                result_file = fullfile(site_results_folder, ...
                    ['LFP_TFR_' sites_tfr(i).site_ID '_' site_conditions(cn).label '.png']);
                lfp_tfa_plot_hs_tuned_tfr(sites_tfr(i).condition(cn).hs_tuned_tfs, ...
                    lfp_tfa_cfg, plottitle, result_file);
            end

        end
        
        % difference between pre- and post-injection
        if sum(lfp_tfa_cfg.compare.perturbations == [0, 1]) > 1
            sites_tfr(i).difference = lfp_tfa_compute_diff_tfr(sites_tfr(i), lfp_tfa_cfg);
            % Plot TFR difference
            for dcn = 1:length(sites_tfr(i).difference)
                if ~isempty(sites_tfr(i).difference(dcn).hs_tuned_tfs)
                    plottitle = ['LFP Diff TFR(Post - Pre): Target ' sites_tfr(i).target, ', Site ', sites_tfr(i).site_ID ];% ', '  ...
                        %'(block ' num2str(cfg_conditions(cn).block) '), '];
                    if site_conditions(dcn).choice == 0
                        plottitle = [plottitle ', Instructed trials'];
                    else
                        plottitle = [plottitle ', Choice trials'];
                    end

                    result_file = fullfile(site_results_folder, ...
                        ['LFP_DiffTFR_' sites_tfr(i).site_ID '_' sites_tfr(i).difference(dcn).label '.png']);
                    lfp_tfa_plot_hs_tuned_tfr(sites_tfr(i).difference(dcn).hs_tuned_tfs, ...
                        lfp_tfa_cfg, plottitle, result_file, 'bluewhitered');
                end
            end
        end
        
        % save mat file for each site
        site_tfr = sites_tfr(i);
        save(fullfile(site_results_folder, ...
            ['LFP_TFR_' site_tfr.site_ID '.mat']), 'site_tfr');
        % save into a mother struct
        session_tfs.sites(i) = site_tfr;       
        
    end
       
    
    % Calculate average TFR across all sites
    session_avg = struct();
    % targets for this session
    targets = unique({states_lfp.target});
    for t = 1:length(targets)
        session_avg(t).target = targets{t};
        for cn = 1:length(site_conditions) 
            session_avg(t).condition(cn).hs_tuned_tfs = [];
            isite = 0;
            for i = 1:length(states_lfp)
                if strcmp(states_lfp(i).target, targets{t})
                    % check if this site should be used for averaging
                    if sites_tfr(i).use_for_avg
                        % calculate the average across sites for this condition 
                        if ~isempty(sites_tfr(i).condition(cn).hs_tuned_tfs) && ... 
                            isfield(sites_tfr(i).condition(cn).hs_tuned_tfs, 'powspctrm')
                            isite = isite + 1;                                
                            % struct to store average tfr across sites
                            %cond_based_tfs(cn).tfs_avg_session = struct();
                            %cell(length(analyse_states),length(hs_labels));

                            %for i = 1:length(cond_based_tfs(cn).sites)

                            for hs = 1:length(hs_labels)
                                for st = 1:size(lfp_tfa_cfg.analyse_states, 1)                        
                                    if ~isempty(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm)
                                        if isite == 1%~isfield(cond_based_tfs(cn).tfs_avg_session(st, hs), 'powspctrm')
                                            session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm = ...
                                                sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm ;
                                            session_avg(t).condition(cn).hs_tuned_tfs(st, hs).time = ...
                                                sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).time;
                                        else
                                            ntimebins = size(session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm, 3);
                                            % average same number of time bins
                                            if ntimebins > length(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).time)
                                                ntimebins = length(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).time);
                                                session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm = ...
                                                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins) + ...
                                                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins) ;
                                            else
                                                if ~isempty(session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm)
                                                    session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm = ...
                                                            session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm + ...
                                                            sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins) ;
                                                else
                                                    session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm = ...
                                                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins) ;
                                                end
                                            end
                                            session_avg(t).condition(cn).hs_tuned_tfs(st, hs).time = ...
                                                sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).time(1:ntimebins);
                                        end

                                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq;
                                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).hs_label = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).hs_label;
                                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).state = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).state;
                                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).state_name = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).state_name;
                                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).cfg = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).cfg;
                                        %session_tfs(t).condition(cn).hs_tuned_tfs(st, hs).nsites = nsites;
                                        session_avg(t).condition(cn).cfg_condition = site_conditions(cn);
                                        session_avg(t).condition(cn).label = site_conditions(cn).label;
                                        session_avg(t).condition(cn).session = states_lfp(i).session;
                                        session_avg(t).condition(cn).target = states_lfp(i).target;
                                    end
                                end
                            end
                        end
                    end
                else
                    continue;
                end            
            end
            % average TFR across sites for a session
            if isfield(session_avg(t).condition(cn).hs_tuned_tfs, 'powspctrm') 
                for hs = 1:length(hs_labels)
                    for st = 1:size(analyse_states,1)
                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm = ...
                            session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm / isite;
                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsites = isite;
                    end
                end
            end
            
            % plot average TFR for this condition and target
            if ~isempty(session_avg(t).condition(cn).hs_tuned_tfs)
                if site_conditions(cn).perturbation == 0
                    injection = 'Pre';
                else
                    injection = 'Post';
                end
                plottitle = ['LFP TFR (' injection '): Target = ' session_avg(t).condition(cn).target ', '  ...
                'Session ', session_avg(t).condition(cn).session, 'Perturbation ' num2str(site_conditions(cn).perturbation_group{1}) ', '];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                else
                    plottitle = [plottitle 'Choice trials'];
                end
                result_file = fullfile(results_folder_tfr, ...
                                ['LFP_TFR_' session_avg(t).condition(cn).target '_'...
                                session_avg(t).condition(cn).session '_' site_conditions(cn).label '.png']);
                lfp_tfa_plot_hs_tuned_tfr(session_avg(t).condition(cn).hs_tuned_tfs, ...
                            lfp_tfa_cfg, plottitle, result_file);

            end
        
        end
        
        % Difference TFR for session
        if sum(lfp_tfa_cfg.compare.perturbations == [0, 1]) > 1
            session_avg(t).difference = lfp_tfa_compute_diff_tfr(session_avg(t), lfp_tfa_cfg);
            % plot average TFR difference across sites for this session
            for dcn = 1:length(session_avg(t).difference)
                if ~isempty(session_avg(t).difference(dcn).hs_tuned_tfs)
                    plottitle = ['LFP Diff TFR (Post - Pre): Target ' session_avg(t).difference(dcn).target ', '  ...
                    'Session ', session_avg(t).difference(dcn).session];
                    if site_conditions(dcn).choice == 0
                        plottitle = [plottitle 'Instructed trials'];
                    else
                        plottitle = [plottitle 'Choice trials'];
                    end
                    result_file = fullfile(results_folder_tfr, ...
                                    ['LFP_DiffTFR_' session_avg(t).difference(dcn).target '_'...
                                    session_avg(t).difference(dcn).session '_' session_avg(t).difference(dcn).label '.png']);
                    lfp_tfa_plot_hs_tuned_tfr(session_avg(t).difference(dcn).hs_tuned_tfs, ...
                                lfp_tfa_cfg, plottitle, result_file, 'bluewhitered');

                end
            end
        end
        
    end
    
    session_tfs.session_avg = session_avg;
    
    
    % save session average tfs
    save(fullfile(results_folder_tfr, ['LFP_TFR_' session_tfs.session '.mat']), 'session_tfs');
    % save settings file
    save(fullfile(results_folder_tfr, 'lfp_tfa_settings.mat'), 'lfp_tfa_cfg');

        % plots

        % site averages
        
    

%         % plot average TFR across sites for this session
%         plottitle = ['Session: ', cond_based_tfs(cn).sites(1).session ', Target = ' cond_based_tfs(cn).sites(1).target ', '  ...
%         'Block ' num2str(cfg_conditions(cn).block) ', '];
%         if cfg_conditions(cn).choice == 0
%             plottitle = [plottitle 'Instructed trials'];
%         else
%             plottitle = [plottitle 'Choice trials'];
%         end
%         results_file = fullfile(cfg_conditions(cn).results_folder, ...
%             ['AvgSpectrogram_' cond_based_tfs(cn).sites(1).session '_' cfg_conditions(cn).label '.png']);
%         lfp_tfa_plot_hs_tuned_tfr(cond_based_tfs(cn).tfs_avg_session, ...
%                     lfp_tfa_cfg, plottitle, results_file);
      
        
    %end
    
    close all;
    % save mat files
%     save(fullfile(results_folder_tfr, 'cfg_conditions.mat'), 'cfg_conditions');
%     save(fullfile(results_folder_tfr, 'cond_based_tfs.mat'), 'cond_based_tfs');
%     save(fullfile(results_folder_tfr, 'lfp_tfa_cfg.mat'), 'lfp_tfa_cfg');
end
        
%         % read baseline tfs from all trials
%         baseline_pow = cell(1, length(states_lfp(i).trials));
%         for t = 1:length(states_lfp(i).trials)
%             trial = states_lfp(i).trials(t);
%             % trials to be considered for baseline calculation
%             consider_trial = trial.choice_trial == choice & trial.ipsilateral == hemisphere & ...
%                 trial.perturbation == 0 & trial.noisy == 0;
%             if consider_trial
% %             if t == 1
% %                 baseline_pow{t} = trial.tfs.powspctrm;
% %             else
%                 baseline_pow{t} = trial.tfs.powspctrm(1,:,...
%                 trial.tfs.time >= trial.baseline.start_t & ...
%                 trial.tfs.time <= trial.baseline.end_t);
%                 %baseline_pow{t} = tfs_baseline.powspctrm + baseline_pow;
%             end
%         end
%         
%         % find mean and std of baseline power for each frequency bin
%         concat_baseline_pow = cat(3, baseline_pow{:});
%         baseline_mean = nanmean(concat_baseline_pow, 3);
%         baseline_std  = nanstd(concat_baseline_pow, 0, 3);
%         
%         % handspace tuned LFP TF analysis
%         tfs_sites = struct();
%         % initialize 
%         for h = 1:length(hs_label)
%             tfs_sites(i).handspace(h) = struct();
%             tfs_sites(i).handspace(h).name = hs_label(h);
%             for s = 1:length(analyse_states)
%                 tfs_sites(i).handspace(h).state(s).name = analyse_states(s);
%                 tfs_sites(i).handspace(h).state(s).tfr_avg = struct();
%                 %tfs_sites(i).handspace(h).state(s).pow_norm = [];
%             end
%         end
%         
%         powspctrm = cell(length(hs_label), length(analyse_states));
%         %powspctrm_norm = length(states_lfp(i).trials);
%         for t = 1:length(states_lfp(i).trials)
%             trial = states_lfp(i).trials(t);
%             % trial to be considered for analysis
%             consider_trial = trial.choice_trial == choice & trial.ipsilateral == hemisphere & ...
%                 isany(trial.block == blocks) & trial.noisy == 0;
%             if consider_trial
%                 % get stored info about all states
%                 states = trial.states;
%                 % handspace tuning for the trial
%                 hndspc = trial.hndspc_lbl;
%                 hs_idx = find(hs_label == hndspc);
%                 for s = 1:length(analyse_states)
%                     st_idx = find(analyse_states(s) == states);
%                     % get the tfs for state and handspace tuning
%                     tfs_sites(i).handspace(hs_idx).state(s).tfr_avg = trial.tfs;
%                     powspctrm{t} = trial.tfs.powspctrm(1,:,...
%                         trial.tfs.time >= trial.states(st_idx).start_t & ...
%                         trial.tfs.time <= trial.states(st_idx).end_t);
%                     tfs_sites(i).handspace(hs_idx).state(s).tfr_avg.time = trial.tfs.time(1, ...
%                         trial.tfs.time >= trial.states(st_idx).start_t & ...
%                         trial.tfs.time <= trial.states(st_idx).end_t);
%                     tfs_sites(i).handspace(hs_idx).state(s).tfr_avg.time = ...
%                         tfs_sites(i).handspace(hs_idx).state(s).tfr_avg.time - ...
%                         trial.states(st_idx).onset_t;
%                     
%                     % baseline nomrlaization
%                     tfs_sites(i).handspace(hs_idx).state(s).tfr_avg.powspctrm = baselineNormalization(...
%                         tfs_sites.handspace(hs_idx).state(s).pow, baseline_mean, baseline_std, 'division');
%                     % combined TFR for all states
%                     tfs_sites(i).handspace(hs_idx).tfr_avg = ...
%                         tfs_sites(i).handspace(hs_idx).state(s).tfr_avg;
%                     
%                 end
%             end
%             
%         end
%         
%         % save tfs for the site
%         states_lfp(i).tfs = tfs_sites(i);
%         
%         % plot TFR for each site and calculate average across sites
%         figure; nsubplot = 0;
%         for h = 1:length(tfs_sites(i).handspace)
%             for s = 1:length(tfs_sites(i).handspace(h).state)
%                 tfr_avg = tfs_sites(i).handspace(h).state(s).tfr_avg;
%                 if ~isempty(tfr_avg)
%                     % configuration for ft_singleplotTFR
%                     cfg = [];
%                     cfg.baseline     = 'no'; %baseline_shift;                 % -400ms to -100ms before the onset of first state
%                     cfg.baselinetype = 'absolute';  
%                     cfg.maskstyle    = 'saturation';
%                     %cfg.interactive  = 'no';
%                     cfg.zlim         = [-2 2];
%                     %cfg.masknans     = 'yes';
%                     cfg.channel  = concat_states_TFR.label;
%                     nsubplot = nsubplot + 1;
%                     subplot(2,4,nsubplot)
%                     ft_singleplotTFR(cfg, concat_states_TFR);
%                     line([0 0], ylim);
%                     title(tfs_sites(i).handspace(h).name + ': ' + ...
%                         tfs_sites(i).handspace(h).state(s).name);
%                 end
%             end
%         end
%     end
%     % calculate the average TFR across sites
%     tfs_avg_sites = struct();
%     % initialize
%     tfs_avg_sites = tfs_sites(1);
%     tfs_powspctrm = cat(3, tfs_sites(:).handspace(h).state(s).tfr_avg.powspctrm);
%     for i = 1:length(tfs_sites)
%         for h = 1:length(tfs_sites(i).handspace)
%             for s = 1:length(tfs_sites(i).handspace(h).state)
%                 tfr_avg = tfs_sites(i).handspace(h).state(s).tfr_avg;
% end
%         
%         
%         
%         
%                     
%         
% %         % baseline for this site
% %         data_lfp_baseline = ft_data_sites(i).baseline;
% %         % get baseline power for normalization
% %         if ~isempty(data_lfp_baseline.trial)
% %             nsample = length(data_lfp_baseline.time{1});
% %             fs = data_lfp_baseline.fsample;
% %             ts = 1/fs;
% %             tstart = data_lfp_baseline.time{1}(1);
% %             tend = data_lfp_baseline.time{1}(end);
% %             % configuration for TFR of baseline
% %             cfg              = [];
% %             cfg.output       = 'pow';
% %             cfg.method       = 'mtmconvol';
% %             cfg.taper        = 'hanning';
% %             cfg.pad          = 'nextpow2';
% %             cfg.foi          = fs/500:fs/500:fs/10;             % analysis 2 to 100 Hz in steps of 2 Hz
% %             cfg.t_ftimwin    = ones(length(cfg.foi),1).*200*ts;    % length of time window = 0.2 sec
% %             %cfg.t_ftimwin    = 4./cfg.foi;    % 7 cycles per time window
% %             cfg.toi          = linspace(tstart, tend, nsample);% time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
% %             cfg.channel      = data_lfp_baseline.label;
% %             disp(['Processing ' cfg.channel]);
% %             TFR_baseline     = ft_freqanalysis(cfg, data_lfp_baseline);
% %             % now save the baseline TFR to corresponding hand space tuning
% %             ft_data_sites(i).TFR_baseline = TFR_baseline;
% %             % get mean value of baseline power for each frequency
% %             baseline_mean = nanmean(TFR_baseline.powspctrm, 3);
% %             % get std value of baseline power for each frequency
% %             baseline_std = nanstd(TFR_baseline.powspctrm, 0, 3);
% %         else
% %             ft_data_sites(i).TFR_baseline = [];
% %         end
% %         
% %         % loop through each hand-space tuning
% %         for hs = 1:length(hs_label)            
% %             % loop through each state
% %             for st = 1:length(ft_data_sites(i).states)
% %                 % now this is a FT_DATATYPE_RAW on which FT time frequency
% %                 % analzsis can be done
% %                 data_lfp_st_hs = ft_data_sites(i).states(st).hndspc(hs);
% %                 nsample = length(ft_data_sites(i).states(st).hndspc(hs).time{1});
% %                 if ~isempty(data_lfp_st_hs.trial)
% %                     fs = data_lfp_st_hs.fsample;
% %                     ts = 1/fs;
% %                     tstart = data_lfp_st_hs.time{1}(1);
% %                     tend = data_lfp_st_hs.time{1}(end);
% %                     % configuration for TFR
% %                     cfg              = [];
% %                     cfg.output       = 'pow';
% %                     cfg.method       = 'mtmconvol';
% %                     cfg.taper        = 'hanning';
% %                     cfg.pad          = 'nextpow2';
% %                     cfg.foi          = fs/500:fs/500:fs/10;             % analysis 2 to 100 Hz in steps of 2 Hz
% %                     cfg.t_ftimwin    = ones(length(cfg.foi),1).*200*ts;    % length of time window = 0.2 sec
% %                     % 7 cycles per time window
% %                     cfg.toi          = linspace(tstart, tend, nsample);% time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
% %                     cfg.channel      = data_lfp_st_hs.label;
% %                     disp(['Processing ' cfg.channel]);
% %                     TFRhann_fixed = ft_freqanalysis(cfg, data_lfp_st_hs); 
% %                     % now save the TFR to corresponding hand space tuning
% %                     ft_data_sites(i).TFR{st, hs} = TFRhann_fixed;
% %                     % baseline normalization
% % %                     TFRhann_fixed_norm = baselinePowerNormalize(TFRhann_fixed, ...
% % %                         TFR_baseline, 'subtraction');
% %                     % subtract baseline mean for each frequency
% %                     TFRhann_fixed_norm = TFRhann_fixed;
% %                     TFRhann_fixed_norm.powspctrm = (TFRhann_fixed.powspctrm - ...
% %                         repmat(baseline_mean, [1 1 size(TFRhann_fixed.powspctrm, 3)])) ./ ...
% %                         repmat(baseline_std, [1 1 size(TFRhann_fixed.powspctrm, 3)]);
% %                     % concatenate TFRs for each state
% %                     if st == 1
% %                         concat_states_TFR = TFRhann_fixed_norm;
% %                     else
% %                         % concatenate times
% %                         concat_states_TFR.time = [concat_states_TFR.time, TFRhann_fixed_norm.time];
% %                         concat_states_TFR.powspctrm = cat(3, concat_states_TFR.powspctrm, TFRhann_fixed_norm.powspctrm);
% %                     end
% % 
% %                     % accumulate TFR per site for averaging across sites
% %                     if i == 1
% %                         TFR_avg{st, hs} = TFRhann_fixed_norm;
% %                         TFR_avg{st, hs}.powspctrm = (1/nsites) * TFRhann_fixed_norm.powspctrm;
% %                         TFR_avg{st, hs}.ntrials = length(data_lfp_st_hs.trial);
% %                     else
% %                         TFR_avg{st, hs}.powspctrm = TFR_avg{st, hs}.powspctrm + (1/nsites) * TFRhann_fixed_norm.powspctrm;
% %                         TFR_avg{st, hs}.ntrials = TFR_avg{st, hs}.ntrials + length(data_lfp_st_hs.trial);
% %                     end
% %                     TFR_avg{st, hs}.label = strcat(analyse_states(st),hs_label(hs));
% %                 else
% %                     ft_data_sites(i).TFR{st, hs} = [];
% %                 end
% %             
% %             end 
% % 
% %             % Visualization
% %             if ~isempty(concat_states_TFR) > 0
% %                 % normalize power spectrogram by z-score
% % %                 concat_states_TFR.powspctrm = (concat_states_TFR.powspctrm - ...
% % %                     nanmean(concat_states_TFR.powspctrm(:))) / nanstd(concat_states_TFR.powspctrm(:));
% %                 time_to_event = concat_states_TFR.time;
% %                 tshift = time_to_event(1);
% %                 % first rearrange the time axis of TFR
% %                 bin_t = concat_states_TFR.time(2) - concat_states_TFR.time(1);
% %                 concat_states_TFR.time = 0:bin_t:bin_t*(length(concat_states_TFR.time)-1);
% %                 % now find the state onsets in the new time axis
% %                 rearr_state_onset_t = concat_states_TFR.time(time_to_event == 0);
% %                 %600*ts%0:1:length(concat_states_TFR.time)-1;
% %                 % rearrange the baseline
% %                 baseline = [-0.4 -0.1];
% %                 baseline_shift = baseline - tshift;
% % 
% %                 cfg = [];
% %                 cfg.baseline     = 'no'; %baseline_shift;                 % -400ms to -100ms before the onset of first state
% %                 cfg.baselinetype = 'absolute';  
% %                 cfg.maskstyle    = 'saturation';
% %                 %cfg.interactive  = 'no';
% %                 cfg.zlim         = [-2 2];
% %                 %cfg.masknans     = 'yes';
% %                 cfg.channel  = concat_states_TFR.label;
% %                 subplot(2,2,hs)
% %                 ft_singleplotTFR(cfg, concat_states_TFR);
% %                 % mark state onsets
% %                 set(gca,'xtick',rearr_state_onset_t)
% %                 set(gca,'xticklabels',strrep(analyse_states, '_', '\_'))
% %                 for t = rearr_state_onset_t
% %                     line([t t], ylim);
% %                 end
% %                 xlabel('Time');
% %                 ylabel('Frequency');
% %                 title(sprintf('%s (ntrials = %g)', cfg.channel{1}, length(data_lfp_st_hs.trial)));
% %                 line([0 0], ylim, 'color', 'k');
% %             end
% %         end
% %         % plot title
% %     %     sgtitle(sprintf('Session: %s, Target = %s, Site: %s (nblocks = %g)', ft_data_sites(i).session, ...
% %     %         ft_data_sites(i).target, ft_data_sites(i).siteID, ft_data_sites(i).nblocks));
% %         plottitle = ['Session: ' ft_data_sites(i).session ', Target = ' ft_data_sites(i).target ', Site ID: ' ft_data_sites(i).siteID ...
% %             '(nblocks = ' num2str(ft_data_sites(i).nblocks) ')'];
% %         ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
% %             , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% %         saveas(gca, fullfile(fig_folder_tfr, [ft_data_sites(i).siteID '.png']));
% %     end
% % 
% %     % save variables
% %     save ft_data_sites ft_data_sites;
% %     %save TFR_avg TFR_avg;
% % 
% %     %%% Visualiyation of time frequency average across site %%%
% %     figure; colormap jet;
% %     for hs = 1:length(hs_label)
% %         for st = 1:length(analyse_states)
% %             % first concatenate TFR for different states for one hs tuning for
% %             % visualization
% %             if st == 1
% %                 concat_TFR_avg = TFR_avg{st,hs};
% %                 concat_TFR_avg.label = hs_label{hs};
% %             else
% %                 concat_TFR_avg.time = [concat_TFR_avg.time TFR_avg{st,hs}.time];
% %                 concat_TFR_avg.powspctrm = cat(3, concat_TFR_avg.powspctrm, TFR_avg{st,hs}.powspctrm);
% %             end
% % 
% %         end
% % 
% %         % Visualization
% % 
% %         time_to_event = concat_TFR_avg.time;
% %         tshift = time_to_event(1);
% %         % first rearrange the time axis of TFR
% %         bin_t = concat_TFR_avg.time(2) - concat_TFR_avg.time(1);
% %         concat_TFR_avg.time = 0:bin_t:bin_t*(length(concat_TFR_avg.time)-1);
% %         % now find the state onsets in the new time axis
% %         rearr_state_onset_t = concat_TFR_avg.time(time_to_event == 0);
% %         % rearrange the baseline
% %         baseline = [-0.4 -0.1];
% %         baseline_shift = baseline - tshift;
% % 
% %         cfg = [];
% %         cfg.baseline     = 'no';%baseline_shift;                 % -400ms to -100ms before the state onset
% %         cfg.baselinetype = 'absolute';  
% %         cfg.maskstyle    = 'saturation';
% %         cfg.zlim         = [-2 2];
% %         cfg.masknans     = 'yes';
% %         %concat_TFR_avg.label = TFR_avg{st,hs}.label{1};
% %         cfg.channel  = concat_TFR_avg.label;
% %         subplot(2,2,hs)
% %         ft_singleplotTFR(cfg, concat_TFR_avg);
% %         % mark state onsets
% %         set(gca,'xtick',rearr_state_onset_t)
% %         set(gca,'xticklabels',strrep(analyse_states, '_', '\_'))
% %         for t = rearr_state_onset_t
% %             line([t t], ylim);
% %         end
% %         xlabel('Time');
% %         ylabel('Frequency');
% %         title(sprintf('%s', cfg.channel{1}));
% %         line([0 0], ylim, 'color', 'k');
% %     end
% % 
% %     % plot title
% %     plottitle = ['Session: ' ft_data_sites(i).session '(nsites = ', ...
% %         num2str(length(ft_data_sites)), 'nblocks = ' num2str(ft_data_sites(i).nblocks) ')'];
% %     ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
% %         , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% %     % sgtitle(sprintf('Session: %s, Target = %s, (nsites = %g)', ft_data_sites(1).session, ...
% %     %     ft_data_sites(1).target(1:end-2), length(ft_data_sites)));
% %     saveas(gca, fullfile(fig_folder_tfr, [ft_data_sites(1).session '.png']));
% % 
% % end
% % 
