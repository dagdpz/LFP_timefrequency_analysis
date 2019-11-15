function [ sitepair_syncspctrm ] = lfp_tfa_sitepair_averaged_syncspctrm( sitepair_crosspow, site_conditions, lfp_tfa_cfg ) 

% lfp_tfa_sitepair_averaged_syncspctrm  - compute and plot the 
% condition-wise average LFP-LFP phase synchronization spectra during a 
% specifed epochs averaged across different trials which 
% satisfy given conditions (A condition is a combination of
% perturbation/choice/type-effector/hand-space tuning). This function calls
% the ft_connectivityanalysis funtion of FieldTrip toolbox to calculate the
% LFP-LFP phase synchronization spectrogram (sync vs. time vs. freq) from 
% the LFP-LFP cross power spectra for the given epoch. Further the sync
% spectrogram value is averaged across the time bins in the epoch to obtain
% a sync spctra (sync vs. freq) for the epoch. 
% PPC(f) = 1/T\sigma_{t=1..T}PPC(t, f)
%
% USAGE:
%	[sitepair_syncspctrm] = lfp_tfa_sitepair_averaged_syncspctrm( ... 
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
%               time windows to analyse, see settings/lfp_tfa_settings_example 
%               ref_hemisphere          - reference hemisphere ('R'/'L')
%               for ipsi- and contra- hand-space labeling
%               sync.measure            - measure to be used for LFP-LFP
%               synchronization (Can only be 'ppc' for now). This setting 
%               is used by the ft_connectivityanalysis routine of the 
%               FieldTrip toolbox to calculate the LFP-LFP phase synchronization
%
% OUTPUTS:
%		sitepair_syncspctrm - condition-wise LFP-LFP sync spectral average 
%                           for the specified time window around the 
%                           specified states (events) across trials for a 
%                           given site pair
%
% REQUIRES:	lfp_tfa_get_condition_trials, ft_connectivityanalysis,
% lfp_tfa_plot_hs_tuned_syncsp
%
% See also settings/lfp_tfa_settings_example, ft_connectivityanalysis, 
% lfp_tfa_compute_sitepair_csd, lfp_tfa_plot_hs_tuned_syncsp, 
% lfp_tfa_sitepair_averaged_sync
    
    % suppress warning for xticklabel
    warning ('off', 'MATLAB:hg:willberemoved');

    % make a folder to save figures
    results_folder_session = lfp_tfa_cfg.session_info...
        (strcmp({lfp_tfa_cfg.session_info.session}, sitepair_crosspow.session)).lfp_syncspctrm_results_fldr;
    if ~exist(results_folder_session, 'dir')
        mkdir(results_folder_session);
    end
       
    % condition based TFS
    sitepair_syncspctrm = struct();
           
        
    % folder to save sitewise results
    sitepair_results_folder = fullfile(results_folder_session, [sitepair_crosspow.sites{1}, ...
        '-', sitepair_crosspow.sites{2}]);
    if ~exist(sitepair_results_folder, 'dir')
        mkdir(sitepair_results_folder);
    end
    % struct to store condition-wise LFP power spectra average
    sitepair_syncspctrm.condition = struct();
    sitepair_syncspctrm.sites = sitepair_crosspow.sites;
    sitepair_syncspctrm.session = sitepair_crosspow.session;
    sitepair_syncspctrm.targets = sitepair_crosspow.targets;
    % flag to indicate if this site should be used for
    % averaging based on minimum no:of trials per condition
    sitepair_syncspctrm.use_for_avg = 1;
        
    % loop through conditions
    for cn = 1:length(site_conditions)

        % hand-space tuning of LFP
        hs_labels = site_conditions(cn).hs_labels;

        % store details of condition
        sitepair_syncspctrm.condition(cn).label = site_conditions(cn).label;
        sitepair_syncspctrm.condition(cn).cfg_condition = site_conditions(cn);
        sitepair_syncspctrm.condition(cn).hs_tuned_syncsp = struct(); 
        sitepair_syncspctrm.condition(cn).ntrials = zeros(1,length(hs_labels));           

        for hs = 1:length(hs_labels)
            % get trial indices for this condition
            % since both site1_lfp and site2_lfp are from same session,
            % trial conditions are same for both sites
            cond_trials = lfp_tfa_get_condition_trials(sitepair_crosspow, site_conditions(cn));
            % get trial indices for this condition and hand-space label
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
            
            % reject noisy trials for both sites
            sitepair_syncspctrm.condition(cn).noisytrials(hs) = ...
                sum(cond_trials & ([sitepair_crosspow.trials.noisy])); 

            % consider only non noisy trials
            fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                & ([sitepair_crosspow.trials.noisy])));
            cond_trials = cond_trials & ~([sitepair_crosspow.trials.noisy]);
            % check if the site contains a specified minimum number
            % of trials for all conditions
            if sum(cond_trials) < lfp_tfa_cfg.mintrials_percondition
                sitepair_syncspctrm.use_for_avg = 0;
            end



            % loop through epochs to analyse
            for ep = 1:size(lfp_tfa_cfg.analyse_epochs, 1)
                epoch_refstate   = lfp_tfa_cfg.analyse_epochs{ep, 1};
                epoch_name       = lfp_tfa_cfg.analyse_epochs{ep, 2};
                epoch_reftstart  = lfp_tfa_cfg.analyse_epochs{ep, 3};
                epoch_reftend    = lfp_tfa_cfg.analyse_epochs{ep, 4};

                epoch_freq = sitepair_crosspow.trials(1).csd;                    
                epoch_freq.powspctrm = {}; % power spectrogram
                epoch_freq.crsspctrm = {}; % power spectrogram
                epoch_freq.time = {}; % timebins fo spectrogram

                for t = find(cond_trials)

                    % get timing information of epoch, same for both sites
                    states          = sitepair_crosspow.trials(t).states;
                    state_onset_t   = states([states(:).id] == ...
                        epoch_refstate).onset_t;
                    epoch_start_t   = states([states(:).id] == ...
                        epoch_refstate).onset_t + epoch_reftstart;
                    epoch_end_t     = states([states(:).id] == ...
                        epoch_refstate).onset_t + epoch_reftend;
                    % sampling frequency
                    fs = sitepair_crosspow.trials(t).fsample; 

                    % epoch power spectra
                    % crop the power spectrum for this epoch
                    epoch_freq.powspctrm = [epoch_freq.powspctrm, ...
                        sitepair_crosspow.trials(t).csd.powspctrm(:, :, ...
                        (sitepair_crosspow.trials(t).csd.time >= epoch_start_t & ...
                        sitepair_crosspow.trials(t).csd.time <= epoch_end_t))];
                    % crop the cross spectrum
                    epoch_freq.crsspctrm = [epoch_freq.crsspctrm, ...
                        sitepair_crosspow.trials(t).csd.crsspctrm(:, :, ...
                        (sitepair_crosspow.trials(t).csd.time >= epoch_start_t & ...
                        sitepair_crosspow.trials(t).csd.time <= epoch_end_t))];
                    % time bins
                    epoch_freq.time = [epoch_freq.time, ...
                        sitepair_crosspow.trials(t).csd.time(1, ...
                        (sitepair_crosspow.trials(t).csd.time >= epoch_start_t & ...
                        sitepair_crosspow.trials(t).csd.time <= epoch_end_t)) - state_onset_t];                     
                    
                end
                
                 % freq bins
                 epoch_freq.freq = sitepair_crosspow.trials(t).csd.freq; 

                % find number of time bins in power
                % spectrogram
                ntimebins = min(cellfun('length', epoch_freq.time));
                % crop each tfs to the ntimebins
                for k = 1:length(epoch_freq.powspctrm)
                    epoch_freq.powspctrm{k} = epoch_freq.powspctrm{k}(:,:,1:ntimebins);
                    epoch_freq.crsspctrm{k} = epoch_freq.crsspctrm{k}(:,:,1:ntimebins);
                    epoch_freq.time{k} = epoch_freq.time{k}(1:ntimebins);
                end
                    
                
                % now concatenate the trials and form a
                % ft_datatype_freq
                epoch_freq.powspctrm = permute(...
                    cat(4, epoch_freq.powspctrm{:}), [4, 1, 2,3]);
                epoch_freq.crsspctrm = permute(...
                    cat(4, epoch_freq.crsspctrm{:}), [4, 1, 2,3]);
                epoch_freq.time = epoch_freq.time{1};
                epoch_freq.label = sitepair_crosspow.sites;
                epoch_freq.dimord = 'rpt_chan_freq_time';

                % calculate the LFP-LFP phase sync
                cfg = [];
                cfg.method     = lfp_tfa_cfg.sync.measure;
                epoch_sync     = ft_connectivityanalysis(cfg, epoch_freq);
                
                % average sync across time bins to get sync spectrum
                epoch_sync = rmfield(epoch_sync, 'time');
                epoch_sync.dimord = 'chancmb_freq';
                epoch_sync.ppcspctrm = nanmean(...
                    epoch_sync.ppcspctrm, 3);

                % save LFP power spectrum
                sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).ppc = epoch_sync;
                sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).trials = find(cond_trials);
                sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).hs_label = hs_labels(hs);
                sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).epoch_name = epoch_name;


            end

        end

        % plot site average power spectrum
        if ~isempty(fieldnames(sitepair_syncspctrm.condition(cn).hs_tuned_syncsp))

            plottitle = sprintf('LFP-LFP Sync spectrum Session: %s, Targets %s-%s (ref: %s), %s', sitepair_syncspctrm.session, sitepair_syncspctrm.targets{:}, ...
                lfp_tfa_cfg.ref_hemisphere, site_conditions(cn).label);
            result_file = fullfile(sitepair_results_folder, ...
                ['LFP-LFP_syncspctrm_' sitepair_syncspctrm.sites{1} '-' sitepair_syncspctrm.sites{2} '_condition' num2str(cn) '.png']); %site_conditions(cn).label
            lfp_tfa_plot_hs_tuned_syncsp(sitepair_syncspctrm.condition(cn).hs_tuned_syncsp, ...
                lfp_tfa_cfg, plottitle, result_file);
        end

    end
    % save mat file for site
    save(fullfile(sitepair_results_folder, ...
        ['LFP-LFP_syncspctrm_' sitepair_syncspctrm.sites{1} '-' sitepair_syncspctrm.sites{2} '.mat']), 'sitepair_syncspctrm');
    %end
    
    close all;        
    
    
%     % Calculate average power spectrum across all sites
%     session_avg = struct();
%     % targets for this session
%     targets = unique({states_lfp.target});
%     for t = 1:length(targets)
%         session_avg(t).target = targets{t};
%         for cn = 1:length(site_conditions)
%             session_avg(t).condition(cn).hs_tuned_power = [];
%             % variable to store no:of sites with trials satisfying this
%             % condition
%             isite = 0;            
%             for i = 1:length(states_lfp)
%                 % if site's target is the target beong considered
%                 if ~strcmp(states_lfp(i).target, targets{t})
%                     continue;
%                 end
%                 if sitepair_syncspctrm(i).use_for_avg                    
%                 % calculate the average LFP power spectrum across sites for this condition 
%                     if ~isempty(sitepair_syncspctrm(i).condition(cn).hs_tuned_power) && ...
%                             isfield(sitepair_syncspctrm(i).condition(cn).hs_tuned_power, 'mean')
%                         isite = isite + 1;
% 
% 
%                         % hand-space labels
%                         for hs = 1:length(hs_labels)
%                             % epochs
%                             for ep = 1:size(lfp_tfa_cfg.analyse_epochs, 1)
%                                 if ~isempty(sitepair_syncspctrm(i).condition(cn).hs_tuned_power(ep, hs).mean)
% 
%                                     if isite == 1
%                                         session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean = ...
%                                             sitepair_syncspctrm(i).condition(cn).hs_tuned_power(ep, hs).mean ;
%                                         session_avg(t).condition(cn).hs_tuned_power(ep, hs).freq = ...
%                                             sitepair_syncspctrm(i).condition(cn).hs_tuned_power(ep, hs).freq;
%                                     else
%                                         if ~isempty(session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean)
%                                             nfreqs = length(session_avg(t).condition(cn).hs_tuned_power(ep, hs).freq);
%                                             % average same number fo frequency
%                                             % bins
%                                             if nfreqs > length(sitepair_syncspctrm(i).condition(cn).hs_tuned_power(ep, hs).freq)
%                                                 nfreqs = length(sitepair_syncspctrm(i).condition(cn).hs_tuned_power(ep, hs).freq);
%                                             end                               
%                                             session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean = ...
%                                                 session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean(1:nfreqs) + ...
%                                                 sitepair_syncspctrm(i).condition(cn).hs_tuned_power(ep, hs).mean(1:nfreqs) ;
%                                         else
%                                             session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean = ...
%                                                 sitepair_syncspctrm(i).condition(cn).hs_tuned_power(ep, hs).mean ;
%                                         end
%                                     end
%                                     % store lfp power spectra average for
%                                     % session
%                                     session_avg(t).condition(cn).hs_tuned_power(ep, hs).hs_label = ...
%                                         sitepair_syncspctrm(i).condition(cn).hs_tuned_power(ep, hs).hs_label;
%                                     session_avg(t).condition(cn).hs_tuned_power(ep, hs).epoch_name = ...
%                                         sitepair_syncspctrm(i).condition(cn).hs_tuned_power(ep, hs).epoch_name;
%                                     session_avg(t).condition(cn).condition = site_conditions(cn);
%                                     session_avg(t).condition(cn).session = states_lfp(i).session;
%                                     session_avg(t).condition(cn).target = states_lfp(i).target;
%                                     session_avg(t).condition(cn).condition = site_conditions(cn);
%                                     session_avg(t).condition(cn).label = site_conditions(cn).label;
%                                 end
%                             end
%                         end
%                     end
%                 end
%             end
%             
%             % average TFR across sites for a session
%             if isfield(session_avg(t).condition(cn).hs_tuned_power, 'mean')
%                 for hs = 1:length(hs_labels)
%                     for ep = 1:size(lfp_tfa_cfg.analyse_epochs, 1)
%                         session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean = ...
%                             session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean / isite;
%                         session_avg(t).condition(cn).hs_tuned_tfs(ep, hs).nsites = isite;
%                     end
%                 end
%             end 
%             
%             % plot average power spectrum across sites for this session
%             if ~isempty(session_avg(t).condition(cn).hs_tuned_power)
%                 plottitle = ['Session: ', session_avg(t).condition(cn).session ...
%                     ', Target = ' session_avg(t).condition(cn).target '(ref_' lfp_tfa_cfg.ref_hemisphere '), '  ...
%                     'Perturb ' num2str(site_conditions(cn).perturbation_group{1}) ', '];
%                 if site_conditions(cn).choice == 0
%                     plottitle = [plottitle 'Instructed trials'];
%                 else
%                     plottitle = [plottitle 'Choice trials'];
%                 end
%                 results_file = fullfile(results_folder_psd, ...
%                     ['LFP_Power_' session_avg(t).condition(cn).session '_' site_conditions(cn).label '.png']);
%                 lfp_tfa_plot_hs_tuned_psd(session_avg(t).condition(cn).hs_tuned_power, ...
%                             lfp_tfa_cfg, plottitle, results_file);
%             end
%             
%         end        
%     end
%     
%     close all;
%     % store session average
%     session_pow.session_avg = session_avg;
%     
%     % save mat files
%     save(fullfile(results_folder_psd, ['LFP_Power_' session_pow.session '.mat']), 'session_pow');
%     % save settings file
%     save(fullfile(results_folder_psd, 'lfp_tfa_settings.mat'), 'lfp_tfa_cfg');
end
        
