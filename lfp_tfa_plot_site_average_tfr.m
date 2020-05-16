function [session_tfs] = lfp_tfa_plot_site_average_tfr( states_lfp, site_conditions, lfp_tfa_cfg ) 
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
    results_folder_tfr = fullfile(lfp_tfa_cfg.session_results_fldr, 'Condition_based_TFS');
    if ~exist(results_folder_tfr, 'dir')
        mkdir(results_folder_tfr);
    end
       
    % condition based TFS
    % struct to store TFR for each site
    sites_tfr = struct();
    % struct to accumulate TFR for each site and store session average
    session_tfs = struct();
    % session name
    session_tfs.session = states_lfp(1).session;
        
    % loop through each site
    for i = 1:length(states_lfp)
        
        %rng(lfp_tfa_cfg.random_seed); % set random seed for reproducibility        
        
        % folder to save sitewise results
        site_results_folder = fullfile(results_folder_tfr, 'sites');
        if ~exist(site_results_folder, 'dir')
            mkdir(site_results_folder);
        end
        
        % structure to store condition-wise tfs
        sites_tfr(i).condition = struct();
        % info about session and site
        sites_tfr(i).site_ID = states_lfp(i).site_ID;
        sites_tfr(i).session = states_lfp(i).session;
        sites_tfr(i).target = states_lfp(i).target;
        % flag to indicate if this site should be used for
        % averaging based on minimum no:of trials per condition
        sites_tfr(i).use_for_avg = 1;
        
        % loop through each condition
        for cn = 1:length(site_conditions)

            % hand-space tuning of LFP
            hs_labels = site_conditions(cn).hs_labels;
                         
            % store details of condition analysed
            sites_tfr(i).condition(cn).label = site_conditions(cn).label;
            sites_tfr(i).condition(cn).cfg_condition = site_conditions(cn);
            sites_tfr(i).condition(cn).hs_tuned_tfs = struct(); 
            sites_tfr(i).condition(cn).ntrials = zeros(1,length(hs_labels));                

            % loop through hand space labels
            for hs = 1:length(hs_labels)
                % get the trial indices which satisfy the given condition
                cond_trials = lfp_tfa_get_condition_trials(states_lfp(i), site_conditions(cn));
                                
                % get the trial indices which satisfy the given hand-space
                % label for the given condition
                if ~strcmp(site_conditions(cn).reach_hands{hs}, 'any')
                    cond_trials = cond_trials & ...
                        strcmp({states_lfp(i).trials.reach_hand}, ...
                        site_conditions(cn).reach_hands{hs});
                end
                if ~strcmp(site_conditions(cn).reach_spaces{hs}, 'any')
                    cond_trials = cond_trials & ...
                        strcmp({states_lfp(i).trials.reach_space}, ...
                        site_conditions(cn).reach_spaces{hs});
                end
                
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
                
                % Remove trials which do not have timing for both windows
                % (e.g no saccades initiation detected
                
                for st = 1:size(lfp_tfa_cfg.analyse_states, 1)
                    
                end
                
                % loop through states to analyse 

                for st = 1:size(lfp_tfa_cfg.analyse_states, 1)
                    
                    if strcmp(lfp_tfa_cfg.analyse_states{st, 1}, 'combined')
                        state_tfs = lfp_tfa_get_combined_tfs(states_lfp(i), ...
                            cond_trials, lfp_tfa_cfg.analyse_states(st, :), lfp_tfa_cfg);
                    end
                    
                    if strcmp(lfp_tfa_cfg.analyse_states{st, 1}, 'single')
                        state_tfs = lfp_tfa_get_state_tfs(states_lfp(i), ...
                            cond_trials, lfp_tfa_cfg.analyse_states(st, :), lfp_tfa_cfg);
                    end                                       

                    if ~isempty(state_tfs.powspctrm)

                        % save average tfs for this condition, hand-space
                        % label, and state
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = state_tfs.powspctrm;
%                         sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.norm_mean = state_tfs.powspctrm_normmean;
%                         sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.raw_mean = state_tfs.powspctrm_rawmean;
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.time = state_tfs.time;
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.freq = state_tfs.freq; 
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.cfg = state_tfs.cfg;
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.baseline = state_tfs.baseline;
                        if isfield(state_tfs, 'state_id') && isfield(state_tfs, 'state_name')
                            sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).state = state_tfs.state_id;
                            sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).state_name = state_tfs.state_name;
                        end
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).hs_label = hs_labels(hs);
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).trials = find(cond_trials);
                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).ntrials = length(find(cond_trials));
                    end

                end

            end

            
            % plot TFR
            if ~isempty(fieldnames(sites_tfr(i).condition(cn).hs_tuned_tfs))
                if site_conditions(cn).perturbation == 0
                    injection = 'Pre';
                elseif site_conditions(cn).perturbation == 1
                    injection = 'Post';
                end
                plottitle = ['LFP TFR (' injection '): Site ' sites_tfr(i).site_ID ...
                    ', Target ' sites_tfr(i).target '(ref_' lfp_tfa_cfg.ref_hemisphere '), '  ...
                    site_conditions(cn).label];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                elseif site_conditions(cn).choice == 1
                    plottitle = [plottitle 'Choice trials'];
                end
                
                result_file = fullfile(site_results_folder, ...
                    ['LFP_TFR_' sites_tfr(i).site_ID '_' site_conditions(cn).label ]);
                lfp_tfa_plot_hs_tuned_tfr_multiple_img(sites_tfr(i).condition(cn).hs_tuned_tfs, ...
                    lfp_tfa_cfg, plottitle, result_file);
            end

        end
        
        sites_tfr(i).difference = [];
        % difference between conditions
        for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
            diff_condition = lfp_tfa_cfg.diff_condition{diff};
            sites_tfr(i).difference = [sites_tfr(i).difference, ...
                lfp_tfa_compute_difference_condition_tfr(sites_tfr(i).condition, diff_condition)];
        end
        % Plot TFR difference
        for dcn = 1:length(sites_tfr(i).difference)
            if ~isempty(fieldnames(sites_tfr(i).difference(dcn).hs_tuned_tfs))
                plottitle = [' Target ' ...
                    sites_tfr(i).target, ' (ref_', lfp_tfa_cfg.ref_hemisphere, ...
                    '), Site ', sites_tfr(i).site_ID ...
                    sites_tfr(i).difference(dcn).label ];
%                     if sites_tfr(i).difference(dcn).cfg_condition.choice == 0
%                         plottitle = [plottitle ', Instructed trials'];
%                     else
%                         plottitle = [plottitle ', Choice trials'];
%                     end

                result_file = fullfile(site_results_folder, ...
                    ['LFP_DiffTFR_' sites_tfr(i).site_ID '_' ...
                    'diff_condition' num2str(dcn) ]);%sites_tfr(i).difference(dcn).label '.png']);
                lfp_tfa_plot_hs_tuned_tfr_multiple_img(sites_tfr(i).difference(dcn).hs_tuned_tfs, ...
                    lfp_tfa_cfg, plottitle, result_file, 'bluewhitered');
            end
        end
        %end
        
        % save mat file for each site
        site_tfr = sites_tfr(i);
        save(fullfile(site_results_folder, ...
            ['LFP_TFR_' site_tfr.site_ID '.mat']), 'site_tfr');
        % save into a mother struct
        session_tfs.sites(i) = site_tfr;
        
        close all;
        
    end
       
    
    % Calculate average TFR across all sites
    session_avg = struct();
    % targets for this session
    targets = unique({states_lfp.target});
    % average each target separately
    for t = 1:length(targets)
        session_avg(t).target = targets{t};
        session_avg(t).session = states_lfp(1).session;
        % loop through conditions
        for cn = 1:length(site_conditions) 
            % condition-wise session average tfs
            session_avg(t).condition(cn).hs_tuned_tfs = struct();
            session_avg(t).condition(cn).cfg_condition = site_conditions(cn);
            session_avg(t).condition(cn).label = site_conditions(cn).label;
            session_avg(t).condition(cn).session = states_lfp(i).session;
            session_avg(t).condition(cn).target = states_lfp(i).target;
            % initialize number of site pairs for each handspace
            % label
            for st = 1:size(sites_tfr(1).condition(cn).hs_tuned_tfs, 1)
                for hs = 1:size(sites_tfr(1).condition(cn).hs_tuned_tfs, 2)
                    session_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsites = 0;
                    session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = [];
                end
            end
            isite = 0;
            
            for i = 1:length(states_lfp)
                if strcmp(states_lfp(i).target, targets{t})
                    % check if this site should be used for averaging
                    if sites_tfr(i).use_for_avg
                        % calculate the average across sites for this condition 
                        if ~isempty(sites_tfr(i).condition(cn).hs_tuned_tfs) && ... 
                            isfield(sites_tfr(i).condition(cn).hs_tuned_tfs, 'freq')% && ...
                            %isite = isite + 1;                                
                            
                            for hs = 1:size(sites_tfr(i).condition(cn).hs_tuned_tfs, 2)
                                for st = 1:size(sites_tfr(i).condition(cn).hs_tuned_tfs, 1)                        
                                    if isfield(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs), 'freq') ...
                                        && ~isempty(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq) ...
                                        && isfield(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq, 'powspctrm') ...
                                        && ~isempty(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm)
                                        
                                        % increment number of site pairs
                                        session_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsites = ...
                                            session_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsites + 1;
                                        
                                        if session_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsites == 1
                                            session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
                                                nanmean(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm, 1);
                                            session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.time = ...
                                                sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.time;
                                            session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.freq = ...
                                                sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.freq;
                                            session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.cfg = ...
                                                sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.cfg;
                                            
                                            session_avg(t).condition(cn).hs_tuned_tfs(st, hs).hs_label = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).hs_label;
                                            if isfield(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs), 'state') && ...
                                                    isfield(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs), 'state_name')
                                                session_avg(t).condition(cn).hs_tuned_tfs(st, hs).state = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).state;
                                                session_avg(t).condition(cn).hs_tuned_tfs(st, hs).state_name = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).state_name;
                                            end
                                            
                                        else
                                            ntimebins = size(session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm, 3);
                                            % average same number of time bins
                                            if ntimebins > length(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.time)
                                                ntimebins = length(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.time);
                                                session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
                                                        cat(1, session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins), ...
                                                        nanmean(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins), 1)) ;
                                            else
                                                if ~isempty(session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm)
                                                    session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
                                                            cat(1, session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm, ...
                                                            nanmean(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins), 1)) ;
                                                else
                                                    session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = ...
                                                        sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins) ;
                                                end
                                            end
                                            session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.time = ...
                                                sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq.time(1:ntimebins);
                                        end
                                        % store session tfs
                                        
                                        %session_avg(t).condition(cn).hs_tuned_tfs(st, hs).cfg = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).cfg;
                                        
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
%             if isfield(session_avg(t).condition(cn).hs_tuned_tfs, 'powspctrm') 
%                 for hs = 1:size(session_avg(t).condition(cn).hs_tuned_tfs, 2)
%                     for st = 1:size(session_avg(t).condition(cn).hs_tuned_tfs, 1)
%                         session_avg(t).condition(cn).hs_tuned_tfs(st, hs).norm_mean = ...
%                             nanmean(session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm, 1);
%                         session_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsites = nsites;
%                     end
%                 end
%             end
            
            % plot average TFR for this condition and target
            if ~isempty(session_avg(t).condition(cn).hs_tuned_tfs)
                if isfield(session_avg(t).condition(cn).hs_tuned_tfs, 'freq') 
                    plottitle = ['LFP TFR: Target = ' session_avg(t).target ...
                        ', (ref_', lfp_tfa_cfg.ref_hemisphere, ') ',  ...
                        'Session ', session_avg(t).condition(cn).session, ...
                        ' ', site_conditions(cn).label];
                    if site_conditions(cn).choice == 0
                        plottitle = [plottitle 'Instructed trials'];
                    elseif site_conditions(cn).choice == 1
                        plottitle = [plottitle 'Choice trials'];
                    end
                    result_file = fullfile(results_folder_tfr, ...
                                    ['LFP_TFR_' session_avg(t).target '_'...
                                    session_avg(t).condition(cn).session '_' site_conditions(cn).label ]);
                    lfp_tfa_plot_hs_tuned_tfr_multiple_img(session_avg(t).condition(cn).hs_tuned_tfs, ...
                                lfp_tfa_cfg, plottitle, result_file);
                end
            end
        
        end
        
        % Difference TFR for session
        % check if both pre- and post- injection blocks exist
        %if sum(lfp_tfa_cfg.compare.perturbations == [0, 1]) > 1
        %session_avg(t).difference = lfp_tfa_compute_diff_tfr(session_avg(t), lfp_tfa_cfg);
        session_avg(t).difference = [];
        % difference between conditions
        for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
            diff_condition = lfp_tfa_cfg.diff_condition{diff};
            session_avg(t).difference = [session_avg(t).difference, ...
                lfp_tfa_compute_difference_condition_tfr(session_avg(t).condition, diff_condition)];
        end

        % plot average TFR difference across sites for this session
        for dcn = 1:length(session_avg(t).difference)
            if ~isempty(fieldnames(session_avg(t).difference(dcn).hs_tuned_tfs))
                plottitle = ['LFP Diff TFR: Target ' session_avg(t).target ...
                    '(ref_' lfp_tfa_cfg.ref_hemisphere '), '  ...
                    'Session ', session_avg(t).difference(dcn).session, ...
                    ' ', session_avg(t).difference(dcn).label];
%                     if session_avg(t).difference(dcn).cfg_condition.choice == 0
%                         plottitle = [plottitle 'Instructed trials'];
%                     else
%                         plottitle = [plottitle 'Choice trials'];
%                     end
                result_file = fullfile(results_folder_tfr, ...
                                ['LFP_DiffTFR_' session_avg(t).target '_' ...
                                session_avg(t).difference(dcn).session '_' ...
                                'diff_condition' num2str(dcn) ]); 
                                %session_avg(t).difference(dcn).label '.png']);
                lfp_tfa_plot_hs_tuned_tfr_multiple_img(session_avg(t).difference(dcn).hs_tuned_tfs, ...
                            lfp_tfa_cfg, plottitle, result_file, 'bluewhitered');

            end
        end
        %end
        
    end
    
    session_tfs.session_avg = session_avg;
    
    % close figures
    close all;    
    
    % save session average tfs
    save(fullfile(results_folder_tfr, ['LFP_TFR_' session_tfs.session '.mat']), 'session_tfs');
    % save settings file
    save(fullfile(results_folder_tfr, 'lfp_tfa_settings.mat'), 'lfp_tfa_cfg');

end
