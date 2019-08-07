function [results_file, sitepair_syncspctrm] = lfp_tfa_sitepair_averaged_syncspctrm( sitepair_crosspow, site_conditions, lfp_tfa_cfg ) 
% lfp_tfa_compute_plot_tfr  - compute and plot average lfp time freq
% response for different hand-space tuning conditions for each site and
% across all sites of a session
%
% USAGE:
%	[ results_file, sitepair_syncspctrm ] = ...
%       lfp_tfa_sitepair_averaged_syncspctrm( sitepair_crosspow, site_conditions, lfp_tfa_cfg )
%
% INPUTS:
%		sitepair_crosspow  	- data structure containing the cross power
%                             spectrogram for a given site pair. i.e., the output of 
%                             lfp_tfa_compute_sitepair_crosspow
%       site_conditions     - structure containing the different trial
%                             conditions to be analysed. i.e., the output of 
%                             lfp_tfa_compare_conditions ( a condition is a combination of
%                             type, effector, choice/instructed, control/inactivation)
%       lfp_tfa_cfg         - struct containing required configurations
%           Required fields:
%               session_info            - structure containing info about
%               the session. This should contain a field
%               'lfp_syncspctrm_results_fldr' where the results of sync
%               spectrum analysis for this sitepair will be stored
%               mintrials_percondition  - minimum number of trials
%               required per condition for considering the site for
%               averaging
%               analyse_epochs          - the epochs at which the sync
%               spectrum should be calculated
%               ref_hemisphere          - the reference hemisphere for
%               labeling ipsi and contra hand and space
%
% OUTPUTS:
%		session_tfs    	- output structure which saves the average tfs for  
%                         trials of a given condition for different handspace 
%                         tunings and periods around the states analysed
%
% REQUIRES:	lfp_tfa_compare_conditions, lfp_tfa_plot_hs_tuned_tfr,
% lfp_tfa_compute_diff_tfr, bluewhitered
%
% See also lfp_tfa_process_lfp, lfp_tfa_compute_baseline_power, 
% lfp_tfa_compare_conditions, lfp_tfa_compute_diff_tfr, 
% lfp_tfa_plot_hs_tuned_tfr, bluewhitered
    
    % suppress warning for xticklabel
    warning ('off', 'MATLAB:hg:willberemoved');

    % make a folder to save figures
    results_folder_session = lfp_tfa_cfg.session_info...
        (strcmp({lfp_tfa_cfg.session_info.session}, sitepair_crosspow.session)).lfp_syncspctrm_results_fldr;
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
    sitepair_syncspctrm = struct();
    % info about session and site pair
    sitepair_syncspctrm.sites = sitepair_crosspow.sites;
    sitepair_syncspctrm.session = sitepair_crosspow.session;
    sitepair_syncspctrm.targets = sitepair_crosspow.targets;
    
    
    % loop through each site
    %for i = 1:length(sitepair_crosspow)
        
        
        
        
        % structure to store condition-wise tfs
        sitepair_syncspctrm.condition = struct();
        
        % flag to indicate if this site should be used for
        % averaging based on minimum no:of trials per condition
        sitepair_syncspctrm.use_for_avg = 1;
        
        % loop through each condition
        for cn = 1:length(site_conditions)

            % hand-space tuning of LFP
            hs_labels = site_conditions(cn).hs_labels;
            
            % store details of condition analysed
            sitepair_syncspctrm.condition(cn).label = site_conditions(cn).label;
            sitepair_syncspctrm.condition(cn).cfg_condition = site_conditions(cn);
            sitepair_syncspctrm.condition(cn).hs_tuned_sync = struct(); 
            sitepair_syncspctrm.condition(cn).ntrials = zeros(1,length(hs_labels));                

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
                
                sitepair_syncspctrm.condition(cn).ntrials(hs) = sum(cond_trials);

                fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
                fprintf('Total number of trials %g\n', sum(cond_trials));
                
                sitepair_syncspctrm.condition(cn).noisytrials(hs) = ...
                    sum(cond_trials & [sitepair_crosspow.trials.noisy]); 

                % consider only non noisy trials
                fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                    & [sitepair_crosspow.trials.noisy]));
                cond_trials = cond_trials & ~[sitepair_crosspow.trials.noisy];

                % check if the site contains a specified minimum number
                % of trials for all conditions
                if sum(cond_trials) < lfp_tfa_cfg.mintrials_percondition
                    sitepair_syncspctrm.use_for_avg = 0;
                end
                % loop through states to analyse 

                for ep = 1:size(lfp_tfa_cfg.analyse_epochs, 1)
                    
                    epoch_state_id = lfp_tfa_cfg.analyse_epochs{ep, 1};
                    epoch_name = lfp_tfa_cfg.analyse_epochs{ep, 2};
                    epoch_ref_tstart = lfp_tfa_cfg.analyse_epochs{ep, 3};
                    epoch_ref_tend = lfp_tfa_cfg.analyse_epochs{ep, 4};
                    
                    epoch_freq = sitepair_crosspow.trials(1).csd;                    
                    epoch_freq.powspctrm = {}; % power spectrogram
                    epoch_freq.crsspctrm = {}; % power spectrogram
                    epoch_freq.time = {}; % timebins fo spectrogram
                    % loop through trials
                    for t = find(cond_trials)
                        % get the state information for this trial
                        states          = sitepair_crosspow.trials(t).states;
                        epoch_onset_t   = states([states(:).id] == ...
                            epoch_state_id).onset_t;
                        epoch_start_t   = states([states(:).id] == ...
                            epoch_state_id).onset_t + epoch_ref_tstart;
                        epoch_end_t     = states([states(:).id] == ...
                            epoch_state_id).onset_t + epoch_ref_tend;
                        % sampling frequency
                        fs = sitepair_crosspow.trials(t).fsample;
                        
                        % sum the power spectrogram for this epoch
                        epoch_freq.powspctrm = [epoch_freq.powspctrm, ...
                            sum(sitepair_crosspow.trials(t).csd.powspctrm(:, :, ...
                            (sitepair_crosspow.trials(t).csd.time >= epoch_start_t & ...
                            sitepair_crosspow.trials(t).csd.time <= epoch_end_t)), 3)];
                        % sum the cross spectrum for this epoch
                        epoch_freq.crsspctrm = [epoch_freq.crsspctrm, ...
                            sum(sitepair_crosspow.trials(t).csd.crsspctrm(:, :, ...
                            (sitepair_crosspow.trials(t).csd.time >= epoch_start_t & ...
                            sitepair_crosspow.trials(t).csd.time <= epoch_end_t)), 3)];                                          

                    end
                    
                    % freq bins
                    epoch_freq.freq = sitepair_crosspow.trials(t).csd.freq; 
                    
                    % now concatenate the trials and form a
                    % ft_datatype_freq
                    epoch_freq.powspctrm = permute(...
                        cat(3, epoch_freq.powspctrm{:}), [3,1,2]);
                    epoch_freq.crsspctrm = permute(...
                        cat(3, epoch_freq.crsspctrm{:}), [3,1,2]);
                    epoch_freq.label = sitepair_crosspow.sites;
                    epoch_freq.dimord = 'rpt_chancmb_freq';      
                    
                    
                    % calculate the LFP-LFP phase sync
                    cfg = [];
                    cfg.method     = lfp_tfa_cfg.sync.measure;
                    epoch_sync     = ft_connectivityanalysis(cfg, epoch_freq);
                    
                    % save csd for this condition and hand-space label
                    if ~isempty(epoch_freq.crsspctrm)
                        epoch_freq = rmfield(epoch_freq, 'powspctrm');
                        sitepair_syncspctrm.condition(cn).hs_tuned_sync(ep, hs).csd = epoch_freq;
                        sitepair_syncspctrm.condition(cn).hs_tuned_sync(ep, hs).csd.crsspctrm_abs_mean = ...
                            squeeze(nanmean(abs(epoch_freq.crsspctrm), 1));
                        sitepair_syncspctrm.condition(cn).hs_tuned_sync(ep, hs).csd.crsspctrm_abs_std = ...
                            squeeze(nanstd(abs(epoch_freq.crsspctrm), 0, 1));
                    end
                    
                    if ~isempty(epoch_sync.ppcspctrm)

                        % save average sync for this condition, hand-space
                        % label, and state
                        sitepair_syncspctrm.condition(cn).hs_tuned_sync(ep, hs).ppc = epoch_sync;
                        sitepair_syncspctrm.condition(cn).hs_tuned_sync(ep, hs).ppc.ppcspctrm = single(epoch_sync.ppcspctrm);
                        sitepair_syncspctrm.condition(cn).hs_tuned_sync(ep, hs).ppc.freq = single(epoch_sync.freq); 
                        sitepair_syncspctrm.condition(cn).hs_tuned_sync(ep, hs).hs_label = hs_labels(hs);
                        sitepair_syncspctrm.condition(cn).hs_tuned_sync(ep, hs).state = epoch_state_id;
                        sitepair_syncspctrm.condition(cn).hs_tuned_sync(ep, hs).state_name = epoch_name;
                        sitepair_syncspctrm.condition(cn).hs_tuned_sync(ep, hs).trials = find(cond_trials);
                        sitepair_syncspctrm.condition(cn).hs_tuned_sync(ep, hs).ntrials = length(find(cond_trials));
                    end

                end

            end

            
            % plot cross spectrum and sync spectrum
            if ~isempty(fieldnames(sitepair_syncspctrm.condition(cn).hs_tuned_sync))
                plottitle = sprintf('LFP-LFP Sync Session: %s, Targets %s-%s (ref: %s), %s', sitepair_syncspctrm.session, sitepair_syncspctrm.targets{:}, ...
                    lfp_tfa_cfg.ref_hemisphere, site_conditions(cn).label);
                result_file = fullfile(sitepair_results_folder, ...
                    ['LFP-LFP_Sync_' sitepair_syncspctrm.sites{1} '-' sitepair_syncspctrm.sites{2} '_condition' num2str(cn) '.png']); %site_conditions(cn).label
                lfp_tfa_plot_hs_tuned_syncsp(sitepair_syncspctrm.condition(cn).hs_tuned_sync, ...
                    lfp_tfa_cfg, plottitle, result_file);
                % cross spectrum
                plottitle = sprintf('LFP-LFP cross-spectrum Session: %s, Targets %s-%s (ref: %s), %s', sitepair_syncspctrm.session, sitepair_syncspctrm.targets{:}, ...
                    lfp_tfa_cfg.ref_hemisphere, site_conditions(cn).label);
                result_file_crs = fullfile(sitepair_results_folder, ...
                    ['LFP-LFP_crsspctrm_' sitepair_syncspctrm.sites{1} '-' sitepair_syncspctrm.sites{2} '_condition' num2str(cn) '.png']); %site_conditions(cn).label
                lfp_tfa_plot_hs_tuned_csd(sitepair_syncspctrm.condition(cn).hs_tuned_sync, ...
                    lfp_tfa_cfg, plottitle, result_file_crs);
            end

        end
        
                
        % save mat file for each site
        results_file = fullfile(sitepair_results_folder, ...
            ['sitepair_syncspctrm_' sitepair_syncspctrm.sites{1} '-' sitepair_syncspctrm.sites{2} '.mat']);
        save(results_file, 'sitepair_syncspctrm');
        % save into a mother struct
        % sitepair_sync.sites(i) = site_tfr;       
        
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
