function [session] = lfp_tfa_compute_session_average( analysis, sites_lfp, site_conditions, lfp_tfa_cfg ) 
% lfp_tfa_plot_site_average_tfr  - compute and plot average lfp time freq
% response for different hand-space tuning and trial conditions for each site and
% across all sites of a session
%
% USAGE:
%	[ session_tfs ] = lfp_tfa_plot_site_average_tfr( states_lfp, ...
%       analyse_states, lfp_tfa_cfg )
%
% INPUTS:
%		states_lfp  	- 1xM struct containing lfp data for all sites of  
%       a session (M = number of sites), see lfp_tfa_process_LFP
%       site_conditions - struct containing the conditions to analyse. 
%       i.e., the output of lfp_tfa_compare_conditions (A condition is a 
%       combination of type-effector, choice, and perturbation)       
%       lfp_tfa_cfg     - struct containing configuration for TFR, see 
%       settings/lfp_tfa_settings_example  
%           Required fields:
%               session_results_fldr            - folder to which the
%               results of the session should be saved
%               random_seed                     - for reproducibility of
%               random numbers, see rng
%               mintrials_percondition          - minimum number of trials
%               required per condition for considering the site for
%               averaging (Condition is a combination of perturbation,
%               choice, type-effector, and hand-space tuning)
%               analyse_states  - cell array containing info about the 
%               states to be analysed
%               ref_hemisphere                  - reference hemispehere for
%               ipsi- and contra-labeling
%               diff_condition                  - trial conditions between
%               which the average LFP TFR difference should be computed
%               baseline_method                 - method used for baseline
%               normalization ('zscore', 'relchange', 'subtraction',
%               'division')
%               compare.perturbations           - perturbations to compare
%               (0 = pre-injection, 1 = post-injection)
%
% OUTPUTS:
%		session_tfs    	- output structure which saves the condition-wise
%                       average tfs for the specified periods around the 
%                       states analysed
%           Fields:
%           sites       - 1xN struct containing condition-wise average tfs 
%                       across trials from a single site
%           session_avg - 1xT struct containing condition-wise average tfs 
%                       averaged across muliple sites in a target area in
%                       one session (T = number of target areas)
%                       
%
% REQUIRES:	lfp_tfa_get_condition_trials, 
% lfp_tfa_get_combined_tfs, lfp_tfa_get_state_tfs, 
% lfp_tfa_compute_difference_condition_tfr
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_define_settings, 
% lfp_tfa_process_lfp, lfp_tfa_compare_conditions, 
% lfp_tfa_compute_difference_condition_tfr, lfp_tfa_plot_site_evoked_LFP, 
% lfp_tfa_plot_site_powspctrum, lfp_tfa_plot_hs_tuned_tfr_multiple_img
    
    % suppress warning for xticklabel
    warning ('off', 'MATLAB:hg:willberemoved');

    % make a folder to save figures
    results_folder = [];
    if strcmp(analysis, 'tfs')
        results_folder = fullfile(lfp_tfa_cfg.session_results_fldr, 'Condition_based_TFS');
    elseif strcmp(analysis, 'pow')
        results_folder = fullfile(lfp_tfa_cfg.session_results_fldr, 'Condition_based_LFP_Power');
    elseif strcmp(analysis, 'evoked')
        results_folder = fullfile(lfp_tfa_cfg.session_results_fldr, 'Condition_based_LFP_Evoked');
    end
    if ~exist(results_folder, 'dir')
        mkdir(results_folder);
    end
       
    % condition based TFS
    % struct to store TFR for each site
    sites_avg = struct();
    % struct to accumulate TFR for each site and store session average
    session = struct();
    % session name
    session.session = sites_lfp(1).session;
        
    % loop through each site
    for i = 1:length(sites_lfp)
        
        %rng(lfp_tfa_cfg.random_seed); % set random seed for reproducibility        
        
        % folder to save sitewise results
        site_results_folder = fullfile(results_folder, 'sites');
        if ~exist(site_results_folder, 'dir')
            mkdir(site_results_folder);
        end
        
        % structure to store condition-wise tfs
        sites_avg(i).condition = struct();
        % info about session and site
        sites_avg(i).site_ID = sites_lfp(i).site_ID;
        sites_avg(i).session = sites_lfp(i).session;
        sites_avg(i).target = sites_lfp(i).target;
        % flag to indicate if this site should be used for
        % averaging based on minimum no:of trials per condition
        sites_avg(i).use_for_avg = 1;
        
        % loop through each condition
        for cn = 1:length(site_conditions)

            % hand-space tuning of LFP
            hs_labels = site_conditions(cn).hs_labels;
                         
            % store details of condition analysed
            sites_avg(i).condition(cn).label = site_conditions(cn).label;
            sites_avg(i).condition(cn).cfg_condition = site_conditions(cn);
            sites_avg(i).condition(cn).hs_tuned_tfs = struct(); 
            sites_avg(i).condition(cn).ntrials = zeros(1,length(hs_labels));                

            % loop through hand space labels
            for hs = 1:length(hs_labels)
                % get the trial indices which satisfy the given condition
                cond_trials = lfp_tfa_get_condition_trials(sites_lfp(i), site_conditions(cn));
                                
                % get the trial indices which satisfy the given hand-space
                % label for the given condition
                if ~strcmp(site_conditions(cn).reach_hands{hs}, 'any')
                    cond_trials = cond_trials & ...
                        strcmp({sites_lfp(i).trials.reach_hand}, ...
                        site_conditions(cn).reach_hands{hs});
                end
                if ~strcmp(site_conditions(cn).reach_spaces{hs}, 'any')
                    cond_trials = cond_trials & ...
                        strcmp({sites_lfp(i).trials.reach_space}, ...
                        site_conditions(cn).reach_spaces{hs});
                end
                
                sites_avg(i).condition(cn).ntrials(hs) = sum(cond_trials);

                fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
                fprintf('Total number of trials %g\n', sum(cond_trials));
                
                sites_avg(i).condition(cn).noisytrials(hs) = ...
                    sum(cond_trials & [sites_lfp(i).trials.noisy]); 

                % consider only non noisy trials
                fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                    & [sites_lfp(i).trials.noisy]));
                cond_trials = cond_trials & ~[sites_lfp(i).trials.noisy];

                % check if the site contains a specified minimum number
                % of trials for all conditions
                if sum(cond_trials) < lfp_tfa_cfg.mintrials_percondition
                    sites_avg(i).use_for_avg = 0;
                end
                
                if strcmp(analysis, 'tfs')
                    % loop through states to analyse 
                    for st = 1:size(lfp_tfa_cfg.analyse_states, 1)

                        if strcmp(lfp_tfa_cfg.analyse_states{st, 1}, 'combined')
                            state_tfs = lfp_tfa_get_combined_tfs(sites_lfp(i), ...
                                cond_trials, lfp_tfa_cfg.analyse_states(st, :), lfp_tfa_cfg);
                        end

                        if strcmp(lfp_tfa_cfg.analyse_states{st, 1}, 'single')
                            state_tfs = lfp_tfa_get_state_tfs(sites_lfp(i), ...
                                cond_trials, lfp_tfa_cfg.analyse_states(st, :), lfp_tfa_cfg);
                        end                                       

                        if ~isempty(state_tfs.powspctrm)

                            % save average tfs for this condition, hand-space
                            % label, and state
                            sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = state_tfs.powspctrm;
    %                         sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.norm_mean = state_tfs.powspctrm_normmean;
    %                         sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.raw_mean = state_tfs.powspctrm_rawmean;
                            sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.time = state_tfs.time;
                            sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.freq = state_tfs.freq; 
                            sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.cfg = state_tfs.cfg;
                            sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.baseline = state_tfs.baseline;
                            if isfield(state_tfs, 'state_id') && isfield(state_tfs, 'state_name')
                                sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).state = state_tfs.state_id;
                                sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).state_name = state_tfs.state_name;
                            end
                            sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).hs_label = hs_labels(hs);
                            sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).trials = find(cond_trials);
                            sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).ntrials = length(find(cond_trials));
                        end

                    end
                elseif strcmp(analysis, 'evoked')
                    % loop through time windows around the states to analyse
                    for st = 1:size(lfp_tfa_cfg.analyse_states, 1)

                        cond_trials_lfp = sites_lfp(i).trials(cond_trials);

                        if strcmp(lfp_tfa_cfg.analyse_states{st, 1}, 'single')
                            state_evoked = lfp_tfa_get_state_evoked_lfp(cond_trials_lfp, ...
                                lfp_tfa_cfg.analyse_states(st, :));
                        elseif strcmp(lfp_tfa_cfg.analyse_states{st, 1}, 'combined')
                            state_evoked = lfp_tfa_get_combined_evoked_lfp(cond_trials_lfp, ...
                                lfp_tfa_cfg.analyse_states(st, :));
                        end                        


                        if ~isempty(state_evoked.lfp)

                            % save evoked LFP
                            sites_avg(i).condition(cn).hs_tuned_evoked(st, hs).lfp = state_evoked.lfp;
                            [sites_avg(i).condition(cn).hs_tuned_evoked(st, hs).mean, ...
                                sites_avg(i).condition(cn).hs_tuned_evoked(st, hs).error] = ...
                                lfp_tfa_compute_statistics(state_evoked.lfp, lfp_tfa_cfg.error_measure);
                            %sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).error = state_tfs.error; 
                            sites_avg(i).condition(cn).hs_tuned_evoked(st, hs).time = state_evoked.lfp_time;
                            sites_avg(i).condition(cn).hs_tuned_evoked(st, hs).dimord = 'ntrials_time';
                            sites_avg(i).condition(cn).hs_tuned_evoked(st, hs).trials = find(cond_trials);
                            sites_avg(i).condition(cn).hs_tuned_evoked(st, hs).hs_label = hs_labels(hs);
                            if isfield(state_evoked, 'state_id') && isfield(state_evoked, 'state_name')
                                sites_avg(i).condition(cn).hs_tuned_evoked(st, hs).state = state_evoked.state_id;
                                sites_avg(i).condition(cn).hs_tuned_evoked(st, hs).state_name = state_evoked.state_name;
                            end

                        end

                    end
                elseif strcmp(analysis, 'pow')
                    % loop through epochs to analyse
                    for ep = 1:size(lfp_tfa_cfg.analyse_epochs, 1)
                        epoch_refstate   = lfp_tfa_cfg.analyse_epochs{ep, 1};
                        epoch_name       = lfp_tfa_cfg.analyse_epochs{ep, 2};
                        epoch_reftstart  = lfp_tfa_cfg.analyse_epochs{ep, 3};
                        epoch_reftend    = lfp_tfa_cfg.analyse_epochs{ep, 4};

                        epoch_pow.psd = {}; % power spectrum
                        epoch_pow.psd_f = {}; % power spectrum freq

                        for t = find(cond_trials)

                            % get timing information of epoch
                            states          = states_lfp(i).trials(t).states;
                            state_onset_t   = states([states(:).id] == ...
                                epoch_refstate).onset_t;
                            epoch_start_t   = states([states(:).id] == ...
                                epoch_refstate).onset_t + epoch_reftstart;
                            epoch_end_t     = states([states(:).id] == ...
                                epoch_refstate).onset_t + epoch_reftend;
                            % sampling frequency
                            fs = states_lfp(i).trials(t).fsample;

                            % LFP power spectrum
                            epoch_powspctrm = sites_lfp(i).trials(t).tfs.powspctrm(1, :, ...
                                (states_lfp(i).trials(t).tfs.time >= epoch_start_t & ...
                                states_lfp(i).trials(t).tfs.time <= epoch_end_t));  
                            epoch_pow.psd = [epoch_pow.psd sum(epoch_powspctrm, 3)];
                            epoch_pow.psd_f = sites_lfp(i).trials(t).tfs.freq;

                        end

                        if ~isempty(epoch_pow.psd)
                            % power spectrum average
                            arr_state_psd = vertcat(epoch_pow.psd{:});
                            %epoch_psd_mean = nanmean(arr_state_psd, 1);

                            % save LFP power spectrum
                            sites_avg(i).condition(cn).hs_tuned_power(ep, hs).psd = arr_state_psd;
                            sites_avg(i).condition(cn).hs_tuned_power(ep, hs).dimord = 'ntrials_freq';
                            [sites_avg(i).condition(cn).hs_tuned_power(ep, hs).mean, ...
                                sites_avg(i).condition(cn).hs_tuned_power(ep, hs).error] = ...
                                lfp_tfa_compute_statistics(arr_state_psd, lfp_tfa_cfg.error_measure);
                            sites_avg(i).condition(cn).hs_tuned_power(ep, hs).freq = epoch_pow.psd_f;
                            sites_avg(i).condition(cn).hs_tuned_power(ep, hs).trials = find(cond_trials);
                            sites_avg(i).condition(cn).hs_tuned_power(ep, hs).hs_label = hs_labels(hs);
                            sites_avg(i).condition(cn).hs_tuned_power(ep, hs).epoch_name = epoch_name;

                        end
                    end
                end

            end

            
            % plot TFR
            if strcmp(analysis, 'tfs') && ~isempty(fieldnames(sites_avg(i).condition(cn).hs_tuned_tfs))
                if site_conditions(cn).perturbation == 0
                    injection = 'Pre';
                elseif site_conditions(cn).perturbation == 1
                    injection = 'Post';
                end
                plottitle = ['LFP TFR (' injection '): Site ' sites_avg(i).site_ID ...
                    ', Target ' sites_avg(i).target '(ref_' lfp_tfa_cfg.ref_hemisphere '), '  ...
                    site_conditions(cn).label];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                elseif site_conditions(cn).choice == 1
                    plottitle = [plottitle 'Choice trials'];
                end
                
                result_file = fullfile(site_results_folder, ...
                    ['LFP_TFR_' sites_avg(i).site_ID '_' site_conditions(cn).label ]);
                lfp_tfa_plot_hs_tuned_tfr_multiple_img(sites_avg(i).condition(cn).hs_tuned_tfs, ...
                    lfp_tfa_cfg, plottitle, result_file);
            end
            % Evoked LFP
            if strcmp(analysis, 'evoked') && ~isempty(fieldnames(sites_avg(i).condition(cn).hs_tuned_evoked))
                plottitle = ['Site ID: ', sites_avg(i).site_ID ', Target = ' ...
                    sites_avg(i).target '(ref_' lfp_tfa_cfg.ref_hemisphere '), '  ...
                    site_conditions(cn).label '), '];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                elseif site_conditions(cn).choice == 1
                    plottitle = [plottitle 'Choice trials'];
                end
                result_file = fullfile(site_results_folder, ...
                    ['LFP_Evoked_' sites_avg(i).site_ID '_' site_conditions(cn).label ]);

                lfp_tfa_plot_evoked_lfp (sites_avg(i).condition(cn).hs_tuned_evoked, lfp_tfa_cfg, ...
                    plottitle, result_file);
            end
            % plot site average power spectrum
            if strcmp(analysis, 'pow') && ~isempty(fieldnames(sites_avg(i).condition(cn).hs_tuned_power))

                plottitle = ['Site ID: ', sites_avg(i).site_ID ...
                    ', Target = ' sites_avg(i).target ' (ref_' lfp_tfa_cfg.ref_hemisphere '), ' ...
                'Perturb ' num2str(site_conditions(cn).perturbation_group{1}) ', '];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                else
                    plottitle = [plottitle 'Choice trials'];
                end
                result_file = fullfile(site_results_folder, ...
                    ['LFP_Power_' sites_avg(i).site_ID '_' site_conditions(cn).label]);
                lfp_tfa_plot_hs_tuned_psd_2(sites_avg(i).condition(cn).hs_tuned_power, ...
                    lfp_tfa_cfg, plottitle, result_file);
            end

        end
        
        if strcmp(analysis, 'tfs')
            sites_avg(i).difference = [];
            % difference between conditions
            for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
                diff_condition = lfp_tfa_cfg.diff_condition{diff};
                sites_avg(i).difference = [sites_avg(i).difference, ...
                    lfp_tfa_compute_difference_condition_tfr(sites_avg(i).condition, diff_condition)];
            end
            % Plot TFR difference
            for dcn = 1:length(sites_avg(i).difference)
                if ~isempty(fieldnames(sites_avg(i).difference(dcn).hs_tuned_tfs))
                    plottitle = [' Target ' ...
                        sites_avg(i).target, ' (ref_', lfp_tfa_cfg.ref_hemisphere, ...
                        '), Site ', sites_avg(i).site_ID ...
                        sites_avg(i).difference(dcn).label ];
    %                     if sites_tfr(i).difference(dcn).cfg_condition.choice == 0
    %                         plottitle = [plottitle ', Instructed trials'];
    %                     else
    %                         plottitle = [plottitle ', Choice trials'];
    %                     end

                    result_file = fullfile(site_results_folder, ...
                        ['LFP_DiffTFR_' sites_avg(i).site_ID '_' ...
                        'diff_condition' num2str(dcn) ]);%sites_tfr(i).difference(dcn).label '.png']);
                    lfp_tfa_plot_hs_tuned_tfr_multiple_img(sites_avg(i).difference(dcn).hs_tuned_tfs, ...
                        lfp_tfa_cfg, plottitle, result_file, 'bluewhitered');
                end
            end
            %end
        end
        
        % save mat file for each site
        %site_tfr = sites_avg(i);
        % save into a mother struct
        session.sites = sites_avg(i);
        
        close all;
        
    end
    
    % calculate average across all sites of the session
    temp = struct();
    temp.session.sites = sites_avg(i);        
    if strcmp(analysis, 'tfr')
        session.session_avg = lfp_tfa_avg_tfr_across_sites(temp, lfp_tfa_cfg, results_folder);
    elseif strcmp(analysis, 'evoked')
        session.session_avg = lfp_tfa_avg_evoked_LFP_across_sites(temp, lfp_tfa_cfg, results_folder);
    elseif strcmp(analysis, 'pow')
        session.session_avg = lfp_tfa_avg_pow_across_sites(temp, lfp_tfa_cfg, results_folder);
    end
    
%     % Calculate average TFR across all sites
%     session_avg = struct();
%     % targets for this session
%     targets = unique({session_lfp.target});
%     % average each target separately
%     for t = 1:length(targets)
%         session_avg(t).target = targets{t};
%         session_avg(t).session = session_lfp(1).session;
%         % loop through conditions
%         for cn = 1:length(site_conditions) 
%             % condition-wise session average tfs
%             session_avg(t).condition(cn).hs_tuned_tfs = struct();
%             session_avg(t).condition(cn).cfg_condition = site_conditions(cn);
%             session_avg(t).condition(cn).label = site_conditions(cn).label;
%             session_avg(t).condition(cn).session = session_lfp(i).session;
%             session_avg(t).condition(cn).target = session_lfp(i).target;
%             % initialize number of site pairs for each handspace
%             % label
%             for st = 1:size(sites_avg(1).condition(cn).hs_tuned_tfs, 1)
%                 for hs = 1:size(sites_avg(1).condition(cn).hs_tuned_tfs, 2)
%                     session_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsites = 0;
%                     session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = [];
%                 end
%             end
%             isite = 0;
%             
%             for i = 1:length(session_lfp)
%                 if strcmp(session_lfp(i).target, targets{t})
%                     % check if this site should be used for averaging
%                     if sites_avg(i).use_for_avg
%                         % calculate the average across sites for this condition 
%                         if ~isempty(sites_avg(i).condition(cn).hs_tuned_tfs) && ... 
%                             isfield(sites_avg(i).condition(cn).hs_tuned_tfs, 'freq')% && ...
%                             %isite = isite + 1;                                
%                             
%                             for hs = 1:size(sites_avg(i).condition(cn).hs_tuned_tfs, 2)
%                                 for st = 1:size(sites_avg(i).condition(cn).hs_tuned_tfs, 1)                        
%                                     if isfield(sites_avg(i).condition(cn).hs_tuned_tfs(st, hs), 'freq') ...
%                                         && ~isempty(sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq) ...
%                                         && isfield(sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq, 'powspctrm') ...
%                                         && ~isempty(sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm)
%                                         
%                                         % increment number of site pairs
%                                         session_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsites = ...
%                                             session_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsites + 1;
%                                         
%                                         if session_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsites == 1
%                                             session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
%                                                 nanmean(sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm, 1);
%                                             session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.time = ...
%                                                 sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.time;
%                                             session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.freq = ...
%                                                 sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.freq;
%                                             session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.cfg = ...
%                                                 sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.cfg;
%                                             
%                                             session_avg(t).condition(cn).hs_tuned_tfs(st, hs).hs_label = sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).hs_label;
%                                             if isfield(sites_avg(i).condition(cn).hs_tuned_tfs(st, hs), 'state') && ...
%                                                     isfield(sites_avg(i).condition(cn).hs_tuned_tfs(st, hs), 'state_name')
%                                                 session_avg(t).condition(cn).hs_tuned_tfs(st, hs).state = sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).state;
%                                                 session_avg(t).condition(cn).hs_tuned_tfs(st, hs).state_name = sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).state_name;
%                                             end
%                                             
%                                         else
%                                             ntimebins = size(session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm, 3);
%                                             % average same number of time bins
%                                             if ntimebins > length(sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.time)
%                                                 ntimebins = length(sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.time);
%                                                 session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
%                                                         cat(1, session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins), ...
%                                                         nanmean(sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins), 1)) ;
%                                             else
%                                                 if ~isempty(session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm)
%                                                     session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
%                                                             cat(1, session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm, ...
%                                                             nanmean(sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins), 1)) ;
%                                                 else
%                                                     session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
%                                                         sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins) ;
%                                                 end
%                                             end
%                                             session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.time = ...
%                                                 sites_avg(i).condition(cn).hs_tuned_tfs(st, hs).freq.time(1:ntimebins);
%                                         end
%                                         % store session tfs
%                                         
%                                         %session_avg(t).condition(cn).hs_tuned_tfs(st, hs).cfg = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).cfg;
%                                         
%                                     end
%                                 end
%                             end
%                         end
%                     end
%                 else
%                     continue;
%                 end            
%             end
%             % average TFR across sites for a session
% %             if isfield(session_avg(t).condition(cn).hs_tuned_tfs, 'powspctrm') 
% %                 for hs = 1:size(session_avg(t).condition(cn).hs_tuned_tfs, 2)
% %                     for st = 1:size(session_avg(t).condition(cn).hs_tuned_tfs, 1)
% %                         session_avg(t).condition(cn).hs_tuned_tfs(st, hs).norm_mean = ...
% %                             nanmean(session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm, 1);
% %                         session_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsites = nsites;
% %                     end
% %                 end
% %             end
%             
%             % plot average TFR for this condition and target
%             if ~isempty(session_avg(t).condition(cn).hs_tuned_tfs)
%                 if isfield(session_avg(t).condition(cn).hs_tuned_tfs, 'freq') 
%                     plottitle = ['LFP TFR: Target = ' session_avg(t).target ...
%                         ', (ref_', lfp_tfa_cfg.ref_hemisphere, ') ',  ...
%                         'Session ', session_avg(t).condition(cn).session, ...
%                         ' ', site_conditions(cn).label];
%                     if site_conditions(cn).choice == 0
%                         plottitle = [plottitle 'Instructed trials'];
%                     elseif site_conditions(cn).choice == 1
%                         plottitle = [plottitle 'Choice trials'];
%                     end
%                     result_file = fullfile(results_folder, ...
%                                     ['LFP_TFR_' session_avg(t).target '_'...
%                                     session_avg(t).condition(cn).session '_' site_conditions(cn).label ]);
%                     lfp_tfa_plot_hs_tuned_tfr_multiple_img(session_avg(t).condition(cn).hs_tuned_tfs, ...
%                                 lfp_tfa_cfg, plottitle, result_file);
%                 end
%             end
%         
%         end
%         
%         % Difference TFR for session
%         % check if both pre- and post- injection blocks exist
%         %if sum(lfp_tfa_cfg.compare.perturbations == [0, 1]) > 1
%         %session_avg(t).difference = lfp_tfa_compute_diff_tfr(session_avg(t), lfp_tfa_cfg);
%         session_avg(t).difference = [];
%         % difference between conditions
%         for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
%             diff_condition = lfp_tfa_cfg.diff_condition{diff};
%             session_avg(t).difference = [session_avg(t).difference, ...
%                 lfp_tfa_compute_difference_condition_tfr(session_avg(t).condition, diff_condition)];
%         end
% 
%         % plot average TFR difference across sites for this session
%         for dcn = 1:length(session_avg(t).difference)
%             if ~isempty(fieldnames(session_avg(t).difference(dcn).hs_tuned_tfs))
%                 plottitle = ['LFP Diff TFR: Target ' session_avg(t).target ...
%                     '(ref_' lfp_tfa_cfg.ref_hemisphere '), '  ...
%                     'Session ', session_avg(t).difference(dcn).session, ...
%                     ' ', session_avg(t).difference(dcn).label];
% %                     if session_avg(t).difference(dcn).cfg_condition.choice == 0
% %                         plottitle = [plottitle 'Instructed trials'];
% %                     else
% %                         plottitle = [plottitle 'Choice trials'];
% %                     end
%                 result_file = fullfile(results_folder, ...
%                                 ['LFP_DiffTFR_' session_avg(t).target '_' ...
%                                 session_avg(t).difference(dcn).session '_' ...
%                                 'diff_condition' num2str(dcn) ]); 
%                                 %session_avg(t).difference(dcn).label '.png']);
%                 lfp_tfa_plot_hs_tuned_tfr_multiple_img(session_avg(t).difference(dcn).hs_tuned_tfs, ...
%                             lfp_tfa_cfg, plottitle, result_file, 'bluewhitered');
% 
%             end
%         end
%         %end
%         
%     end
%     
%     session_avg.session_avg = session_avg;
    
    % close figures
    close all;    
    
    % save session average tfs
    if strcmp(analysis, 'tfs')
        session_tfs = session;
        save(fullfile(results_folder, ['LFP_TFR_' session_tfs.session '.mat']), 'session_tfs');
    elseif strcmp(analysis, 'evoked')
        session_evoked = session;
        save(fullfile(results_folder, ['LFP_evoked_' session_evoked.session '.mat']), 'session_evoked');
    elseif strcmp(analysis, 'pow')
        session_pow = session;
        save(fullfile(results_folder, ['LFP_Power_' session_pow.session '.mat']), 'session_pow');
    end
    % save settings file
    %save(fullfile(results_folder, 'lfp_tfa_settings.mat'), 'lfp_tfa_cfg');

end
