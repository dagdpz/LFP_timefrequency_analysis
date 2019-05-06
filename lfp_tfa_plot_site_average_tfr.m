function [session_tfs] = lfp_tfa_plot_site_average_tfr( states_lfp, analyse_states, lfp_tfa_cfg ) 
% lfp_tfa_compute_plot_tfr  - compute and plot average lfp time freq response for
% different hand-space tuning conditions for each site and across all sites
% of a session
%
% USAGE:
%	[ session_tfs ] = lfp_tfa_plot_average_tfr( sites_lfp_folder, analyse_states, lfp_tfa_cfg )
%
% INPUTS:
%		states_lfp  	- folder containing lfp data for all sites of  
%       a session, output from lfp_tfa_compute_baseline
%       analyse_states  - cell array containing states to be
%       analysed, see lfp_tfa_settings
%       lfp_tfa_cfg     - struct containing configuration for TFR 
%           Required fields:
%               session_results_fldr            - folder to which the
%               results of the session should be saved
%               mintrials_percondition          - minimum number of trials
%               required per condition for considering the site for
%               averaging
%               analyse_states                  - states to analyse 
%               baseline_method                 - method used for baseline
%               normalization ('zscore', 'relchange', 'subtraction',
%               'division')
%               compare.perturbations           - perturbations to compare
%               (0 = pre-injection, 1 = post-injection)
%
% OUTPUTS:
%		session_tfs    	- output structure which saves the average tfs for  
%                         trials of a given condition for different handspace 
%                         tunings and periods around the states analysed
%
% REQUIRES:	lfp_tfa_compare_conditions, lfp_tfa_plot_hs_tuned_tfr,
% lfp_tfa_compute_diff_tfr, bluewhitered
%
% See also lfp_tfa_process_lfp, lfp_tfa_compute_baseline, 
% lfp_tfa_compare_conditions, lfp_tfa_compute_diff_tfr, 
% lfp_tfa_plot_hs_tuned_tfr, bluewhitered
    
    % suppress warning for xticklabel
    warning ('off', 'MATLAB:hg:willberemoved');

    % make a folder to save figures
    %sessionName = states_lfp(1).session;
    results_folder_tfr = fullfile(lfp_tfa_cfg.session_results_fldr, 'Condition_based_TFS');
    if ~exist(results_folder_tfr, 'dir')
        mkdir(results_folder_tfr);
    end
       
    % condition based TFS
    sites_tfr = struct();
    session_tfs = struct();
    session_tfs.session = states_lfp(1).session;
    for i = 1:length(states_lfp)
        
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

end
