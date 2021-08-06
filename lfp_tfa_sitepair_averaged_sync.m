function [sitepair_sync] = lfp_tfa_sitepair_averaged_sync( sitepair_crosspow, site_conditions, lfp_tfa_cfg ) 
% lfp_tfa_sitepair_averaged_sync  - compute and plot the condition-wise average
% LFP-LFP phase synchronization spectrogram during a specifed time window
% of interest averaged across different trials which satisfy given 
% conditions (A condition is a combination of
% perturbation/choice/type-effector/hand-space tuning). This function calls
% the ft_connectivityanalysis funtion of FieldTrip toolbox to calculate the
% LFP-LFP phase synchronization spectrogram from the LFP-LFP cross power 
% spectra
%
% USAGE:
%	[sitepair_sync] = lfp_tfa_sitepair_averaged_sync( ... 
%       sitepair_crosspow, site_conditions, lfp_tfa_cfg )
%
% INPUTS:
%		sitepair_crosspow  	- struct containing LFP-LFP cross power
%                           spectrum for a pair of sites. i.e., the output
%                           of lfp_tfa_compute_sitepair_csd
%       
%		site_conditions     - strcuture containing the trial conditions to 
%                           be analysed (Condition is a combination of
%                           perturbation/choice/type-effector/hand-space
%                           tuning)
%       
%       lfp_tfa_cfg         - struct containing configuration for sync
%                           calculation
%           Required fields:
%               analyse_lfp_folder      - root folder under which the
%               results of LFP-LFP sync should be saved. The results of the
%               session will be saved in a sub-folder
%               [sitepair_crosspow.session, '/Condition_based_Sync'] under
%               this root folder
%               mintrials_percondition  - minimum number of trials
%               required per condition for considering the site for
%               averaging
%               analyse_states          - cell array containing states and 
%               time windows to analyse 
%               ref_hemisphere          - reference hemisphere ('R'/'L')
%               for ipsi- and contra- hand-space labeling
%               sync.measure            - measure to be used for LFP-LFP
%               synchronization (Can only be 'ppc' for now). This setting 
%               is used by the ft_connectivityanalysis routine of the 
%               FieldTrip toolbox to calculate the LFP-LFP phase synchronization
%               diff_condition          - conditions between which
%               difference in LFP-LFP sync should be computed and plotted
%
% OUTPUTS:
%		sitepair_sync  	- condition-wise LFP-LFP sync average for the 
%                         specified time window around the specified states
%                         (events) across trials for a given site pair
%
% REQUIRES:	lfp_tfa_get_condition_trials, ft_connectivityanalysis,
% lfp_tfa_compute_diff_condition_tfsync, lfp_tfa_plot_hs_tuned_sync
%
% See also settings/lfp_tfa_settings_example, ft_connectivityanalysis, 
% lfp_tfa_compute_sitepair_csd, lfp_tfa_compute_diff_condition_tfsync, 
% lfp_tfa_plot_hs_tuned_sync, lfp_tfa_sitepair_averaged_syncspctrm
    
    % suppress warning for xticklabel
    warning ('off', 'MATLAB:hg:willberemoved');

    % make a folder to save figures
    results_folder_session = fullfile(lfp_tfa_cfg.analyse_lfp_folder, sitepair_crosspow.session, 'Condition_based_Sync');
    if ~exist(results_folder_session, 'dir')
        mkdir(results_folder_session);
    end
    
    % folder to save sitewise results
    sitepair_results_folder = fullfile(results_folder_session, [sitepair_crosspow.sites{1}, ...
        '-', sitepair_crosspow.sites{2}]);
    if ~exist(sitepair_results_folder, 'dir')
        mkdir(sitepair_results_folder);
    end
       
    % condition based sync
    % struct to accumulate sync for each site and store session average
    sitepair_sync = struct();
    % info about session and site pair
    sitepair_sync.sites = sitepair_crosspow.sites;
    sitepair_sync.session = sitepair_crosspow.session;
    sitepair_sync.targets = sitepair_crosspow.targets;
    
    
    % loop through each site
    %for i = 1:length(sitepair_crosspow)
        
        
        
        
        % structure to store condition-wise tfs
        sitepair_sync.condition = struct();
        
        % flag to indicate if this site should be used for
        % averaging based on minimum no:of trials per condition
        sitepair_sync.use_for_avg = 1;
        
        % loop through each condition
        for cn = 1:length(site_conditions)

            % hand-space tuning of LFP
            hs_labels = site_conditions(cn).hs_labels;
            
            % store details of condition analysed
            sitepair_sync.condition(cn).label = site_conditions(cn).label;
            sitepair_sync.condition(cn).cfg_condition = site_conditions(cn);
            sitepair_sync.condition(cn).hs_tuned_sync = struct(); 
            sitepair_sync.condition(cn).ntrials = zeros(1,length(hs_labels));                

            % loop through hand space labels
            for hs = 1:length(hs_labels)
                % get the trial indices which satisfy the given condition
                cond_trials = lfp_tfa_get_condition_trials(sitepair_crosspow, site_conditions(cn));
                % get the trial indices which satisfy the given hand-space
                % label for the given condition
                if ~strcmp(site_conditions(cn).reach_hands{hs}, 'any')
                    cond_trials = cond_trials & ...
                        strcmp({sitepair_crosspow.trials.reach_hand}, ...
                        site_conditions(cn).reach_hands{hs});
                end
                if ~strcmp(site_conditions(cn).reach_spaces{hs}, 'any')
                    cond_trials = cond_trials & ...
                        strcmp({sitepair_crosspow.trials.reach_space}, ...
                        site_conditions(cn).reach_spaces{hs});
                end
                
                sitepair_sync.condition(cn).ntrials(hs) = sum(cond_trials);

                fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
                fprintf('Total number of trials %g\n', sum(cond_trials));
                
                sitepair_sync.condition(cn).noisytrials(hs) = ...
                    sum(cond_trials & [sitepair_crosspow.trials.noisy]); 

                % consider only non noisy trials
                fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                    & [sitepair_crosspow.trials.noisy]));
                cond_trials = cond_trials & ~[sitepair_crosspow.trials.noisy];

                % check if the site contains a specified minimum number
                % of trials for all conditions
                if sum(cond_trials) < lfp_tfa_cfg.mintrials_percondition
                    sitepair_sync.use_for_avg = 0;
                end
                % loop through states to analyse 

                for st = 1:size(lfp_tfa_cfg.analyse_states, 1)
                    
                    if strcmp(lfp_tfa_cfg.analyse_states{st, 1}, 'combined')
                        continue;
                    end
                    
                    state_id = lfp_tfa_cfg.analyse_states{st, 2};
                    state_name = lfp_tfa_cfg.analyse_states{st, 3};
                    state_ref_tstart = lfp_tfa_cfg.analyse_states{st, 4};
                    state_ref_tend = lfp_tfa_cfg.analyse_states{st, 5};
                    
                    state_freq = sitepair_crosspow.trials(1).csd;                    
                    state_freq.powspctrm = {}; % power spectrogram
                    state_freq.crsspctrm = {}; % power spectrogram
                    state_freq.time = {}; % timebins fo spectrogram
                    % loop through trials
                    for t = find(cond_trials)
                        % get the state information for this trial
                        states          = sitepair_crosspow.trials(t).states;
                        state_onset_t   = states([states(:).id] == ...
                            state_id).onset_t;
                        state_start_t   = states([states(:).id] == ...
                            state_id).onset_t + state_ref_tstart;
                        state_end_t     = states([states(:).id] == ...
                            state_id).onset_t + state_ref_tend;
                        % sampling frequency
                        fs = sitepair_crosspow.trials(t).fsample;
                        
                        % crop the power spectrum for this state
                        state_freq.powspctrm = [state_freq.powspctrm, ...
                            sitepair_crosspow.trials(t).csd.powspctrm(:, :, ...
                            (sitepair_crosspow.trials(t).csd.time >= state_start_t & ...
                            sitepair_crosspow.trials(t).csd.time <= state_end_t))];
                        % crop the cross spectrum
                        state_freq.crsspctrm = [state_freq.crsspctrm, ...
                            sitepair_crosspow.trials(t).csd.crsspctrm(:, :, ...
                            (sitepair_crosspow.trials(t).csd.time >= state_start_t & ...
                            sitepair_crosspow.trials(t).csd.time <= state_end_t))];
                        % time bins
                        state_freq.time = [state_freq.time, ...
                            sitepair_crosspow.trials(t).csd.time(1, ...
                            (sitepair_crosspow.trials(t).csd.time >= state_start_t & ...
                            sitepair_crosspow.trials(t).csd.time <= state_end_t)) - state_onset_t];                       

                    end
                    
                    % freq bins
                    state_freq.freq = sitepair_crosspow.trials(t).csd.freq; 
                    
                    % find number of time bins in power
                    % spectrogram
                    ntimebins = min(cellfun('length', state_freq.time));
                    % crop each tfs to the ntimebins
                    for k = 1:length(state_freq.powspctrm)
                        state_freq.powspctrm{k} = state_freq.powspctrm{k}(:,:,1:ntimebins);
                        state_freq.crsspctrm{k} = state_freq.crsspctrm{k}(:,:,1:ntimebins);
                        state_freq.time{k} = state_freq.time{k}(1:ntimebins);
                    end
                    
                    % now concatenate the trials and form a
                    % ft_datatype_freq
                    state_freq.powspctrm = permute(...
                        cat(4, state_freq.powspctrm{:}), [4, 1, 2,3]);
                    state_freq.crsspctrm = permute(...
                        cat(4, state_freq.crsspctrm{:}), [4, 1, 2,3]);
                    state_freq.time = state_freq.time{1};
                    state_freq.label = sitepair_crosspow.sites;
                    state_freq.dimord = 'rpt_chan_freq_time';
                    
                    % calculate the LFP-LFP phase sync
                    cfg = [];
                    cfg.method     = lfp_tfa_cfg.sync.measure;
                    state_sync     = ft_connectivityanalysis(cfg, state_freq);
                    
                    if ~isempty(state_sync.ppcspctrm)

                        % save average sync for this condition, hand-space
                        % label, and state
                        sitepair_sync.condition(cn).hs_tuned_sync(st, hs).ppc = state_sync;
                        sitepair_sync.condition(cn).hs_tuned_sync(st, hs).ppc.ppcspctrm = single(state_sync.ppcspctrm);
                        sitepair_sync.condition(cn).hs_tuned_sync(st, hs).time = single(state_sync.time);
                        sitepair_sync.condition(cn).hs_tuned_sync(st, hs).freq = single(state_sync.freq); 
                        sitepair_sync.condition(cn).hs_tuned_sync(st, hs).cfg = state_sync.cfg;
                        sitepair_sync.condition(cn).hs_tuned_sync(st, hs).hs_label = hs_labels(hs);
                        sitepair_sync.condition(cn).hs_tuned_sync(st, hs).state = state_id;
                        sitepair_sync.condition(cn).hs_tuned_sync(st, hs).state_name = state_name;
                        sitepair_sync.condition(cn).hs_tuned_sync(st, hs).trials = find(cond_trials);
                        sitepair_sync.condition(cn).hs_tuned_sync(st, hs).ntrials = length(find(cond_trials));
                    end

                end

            end

            
            % plot TFR
            if ~isempty(fieldnames(sitepair_sync.condition(cn).hs_tuned_sync))
                plottitle = sprintf('LFP-LFP Sync Session: %s, Targets %s-%s (ref: %s), %s', sitepair_sync.session, sitepair_sync.targets{:}, ...
                    lfp_tfa_cfg.ref_hemisphere, site_conditions(cn).label);
                result_file = fullfile(sitepair_results_folder, ...
                    ['LFP-LFP_Sync_' sitepair_sync.sites{1} '-' sitepair_sync.sites{2} '_condition' num2str(cn) ]); %site_conditions(cn).label
%                 lfp_tfa_plot_hs_tuned_sync(sitepair_sync.condition(cn).hs_tuned_sync, ...
%                     lfp_tfa_cfg, plottitle, result_file, 'imscale', [0, 1]);
            end

        end
        
        sitepair_sync.difference = [];
        % difference between conditions
        for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
            diff_condition = lfp_tfa_cfg.diff_condition{diff};
            sitepair_sync.difference = [sitepair_sync.difference, ...
                lfp_tfa_compute_diff_condition_tfsync(sitepair_sync.condition, diff_condition)];
        end
        % Plot TFR difference
        for dcn = 1:length(sitepair_sync.difference)
            if ~isempty(fieldnames(sitepair_sync.difference(dcn).hs_tuned_sync))
                plottitle = sprintf('LFP Diff Sync Session: %s, Targets %s-%s (ref: %s), %s', sitepair_sync.session, sitepair_sync.targets{:}, ...
                    lfp_tfa_cfg.ref_hemisphere, sitepair_sync.difference(dcn).label );
%                     if sites_tfr(i).difference(dcn).cfg_condition.choice == 0
%                         plottitle = [plottitle ', Instructed trials'];
%                     else
%                         plottitle = [plottitle ', Choice trials'];
%                     end

                result_file = fullfile(sitepair_results_folder, ...
                    ['LFP_DiffSync_' sitepair_sync.sites{1} '-' sitepair_sync.sites{2} '_' ...
                    'diff_condition' num2str(dcn)]);%sites_tfr(i).difference(dcn).label '.png']);
                lfp_tfa_plot_hs_tuned_sync(sitepair_sync.difference(dcn).hs_tuned_sync, ...
                    lfp_tfa_cfg, plottitle, result_file, 'cmap', 'bluewhitered', 'imscale', [-0.3, 0.3]);
            end
        end
        %end
        
        % save mat file for each site
        save(fullfile(sitepair_results_folder, ...
            ['sitepair_sync_' sitepair_sync.sites{1} '-' sitepair_sync.sites{2} '.mat']), 'sitepair_sync');
        % save into a mother struct
        % sitepair_sync.sites(i) = site_tfr;    
        
        close all;
        
    %end
       
    
%     % Calculate average TFR across all sites
%     session_avg = struct();
%     % targets for this session
%     targets = unique({sitepair_crosspow.target});
%     % average each target separately
%     for t = 1:length(targets)
%         session_avg(t).target = targets{t};
%         session_avg(t).session = sitepair_crosspow(1).session;
%         % loop through conditions
%         for cn = 1:length(site_conditions) 
%             % condition-wise session average tfs
%             session_avg(t).condition(cn).hs_tuned_tfs = [];
%             isite = 0;
%             for i = 1:length(sitepair_crosspow)
%                 if strcmp(sitepair_crosspow(i).target, targets{t})
%                     % check if this site should be used for averaging
%                     if sites_tfr(i).use_for_avg
%                         % calculate the average across sites for this condition 
%                         if ~isempty(sites_tfr(i).condition(cn).hs_tuned_tfs) && ... 
%                             isfield(sites_tfr(i).condition(cn).hs_tuned_tfs, 'powspctrm')
%                             isite = isite + 1;                                
%                             
%                             for hs = 1:length(hs_labels)
%                                 for st = 1:size(lfp_tfa_cfg.analyse_states, 1)                        
%                                     if ~isempty(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm)
%                                         if isite == 1
%                                             session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm = ...
%                                                 sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm ;
%                                             session_avg(t).condition(cn).hs_tuned_tfs(st, hs).time = ...
%                                                 sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).time;
%                                         else
%                                             ntimebins = size(session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm, 3);
%                                             % average same number of time bins
%                                             if ntimebins > length(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).time)
%                                                 ntimebins = length(sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).time);
%                                                 session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm = ...
%                                                         session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins) + ...
%                                                         sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins) ;
%                                             else
%                                                 if ~isempty(session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm)
%                                                     session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm = ...
%                                                             session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm + ...
%                                                             sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins) ;
%                                                 else
%                                                     session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm = ...
%                                                         sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins) ;
%                                                 end
%                                             end
%                                             session_avg(t).condition(cn).hs_tuned_tfs(st, hs).time = ...
%                                                 sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).time(1:ntimebins);
%                                         end
%                                         % store session tfs
%                                         session_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).freq;
%                                         session_avg(t).condition(cn).hs_tuned_tfs(st, hs).hs_label = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).hs_label;
%                                         session_avg(t).condition(cn).hs_tuned_tfs(st, hs).state = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).state;
%                                         session_avg(t).condition(cn).hs_tuned_tfs(st, hs).state_name = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).state_name;
%                                         session_avg(t).condition(cn).hs_tuned_tfs(st, hs).cfg = sites_tfr(i).condition(cn).hs_tuned_tfs(st, hs).cfg;
%                                         session_avg(t).condition(cn).cfg_condition = site_conditions(cn);
%                                         session_avg(t).condition(cn).label = site_conditions(cn).label;
%                                         session_avg(t).condition(cn).session = sitepair_crosspow(i).session;
%                                         session_avg(t).condition(cn).target = sitepair_crosspow(i).target;
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
%             if isfield(session_avg(t).condition(cn).hs_tuned_tfs, 'powspctrm') 
%                 for hs = 1:length(hs_labels)
%                     for st = 1:size(analyse_states,1)
%                         session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm = ...
%                             session_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm / isite;
%                         session_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsites = isite;
%                     end
%                 end
%             end
%             
%             % plot average TFR for this condition and target
%             if ~isempty(session_avg(t).condition(cn).hs_tuned_tfs)
%                 if site_conditions(cn).perturbation == 0
%                     injection = 'Pre';
%                 else
%                     injection = 'Post';
%                 end
%                 plottitle = ['LFP TFR (' injection '): Target = ' session_avg(t).target ...
%                     ', (ref_', lfp_tfa_cfg.ref_hemisphere, ') ',  ...
%                     'Session ', session_avg(t).condition(cn).session, ...
%                     ' ', site_conditions(cn).label];
%                 if site_conditions(cn).choice == 0
%                     plottitle = [plottitle 'Instructed trials'];
%                 else
%                     plottitle = [plottitle 'Choice trials'];
%                 end
%                 result_file = fullfile(results_folder_sync, ...
%                                 ['LFP_TFR_' session_avg(t).target '_'...
%                                 session_avg(t).condition(cn).session '_' site_conditions(cn).label '.png']);
%                 lfp_tfa_plot_hs_tuned_tfr_multiple_img(session_avg(t).condition(cn).hs_tuned_tfs, ...
%                             lfp_tfa_cfg, plottitle, result_file);
% 
%             end
%         
%         end
%         
%         % Difference TFR for session
%         % check if both pre- and post- injection blocks exist
%         if sum(lfp_tfa_cfg.compare.perturbations == [0, 1]) > 1
%             %session_avg(t).difference = lfp_tfa_compute_diff_tfr(session_avg(t), lfp_tfa_cfg);
%             session_avg(t).difference = [];
%             % difference between conditions
%             for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
%                 diff_condition = lfp_tfa_cfg.diff_condition{diff};
%                 session_avg(t).difference = [session_avg(t).difference, ...
%                     lfp_tfa_compute_difference_condition_tfr(session_avg(t).condition, diff_condition)];
%             end
%             
%             % plot average TFR difference across sites for this session
%             for dcn = 1:length(session_avg(t).difference)
%                 if ~isempty(session_avg(t).difference(dcn).hs_tuned_tfs)
%                     plottitle = ['LFP Diff TFR: Target ' session_avg(t).target ...
%                         '(ref_' lfp_tfa_cfg.ref_hemisphere '), '  ...
%                         'Session ', session_avg(t).difference(dcn).session, ...
%                         ' ', session_avg(t).difference(dcn).label];
% %                     if session_avg(t).difference(dcn).cfg_condition.choice == 0
% %                         plottitle = [plottitle 'Instructed trials'];
% %                     else
% %                         plottitle = [plottitle 'Choice trials'];
% %                     end
%                     result_file = fullfile(results_folder_sync, ...
%                                     ['LFP_DiffTFR_' session_avg(t).target '_' ...
%                                     session_avg(t).difference(dcn).session '_' ...
%                                     'diff_condition' num2str(dcn) '.png']); 
%                                     %session_avg(t).difference(dcn).label '.png']);
%                     lfp_tfa_plot_hs_tuned_tfr_multiple_img(session_avg(t).difference(dcn).hs_tuned_tfs, ...
%                                 lfp_tfa_cfg, plottitle, result_file, 'bluewhitered');
% 
%                 end
%             end
%         end
%         
%     end
%     
%     sitepair_sync.session_avg = session_avg;
%     
%     % close figures
%     close all;    
%     
%     % save session average tfs
%     save(fullfile(results_folder_sync, ['LFP_TFR_' sitepair_sync.session '.mat']), 'session_tfs');
%     % save settings file
%     save(fullfile(results_folder_sync, 'lfp_tfa_settings.mat'), 'lfp_tfa_cfg');

end
