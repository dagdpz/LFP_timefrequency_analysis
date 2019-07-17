function [ sitepair_syncspctrm ] = lfp_tfa_sitepair_avg_syncspctrum( site1_lfp, site2_lfp, site_conditions, lfp_tfa_cfg ) 

% lfp_tfa_plot_average_powspctrum  - calculate and plot the average lfp power spectrum for
% different hand-space tuning conditions for each site and across sites for
% a session
%
% USAGE:
%	[ session_pow ] = lfp_tfa_plot_average_powspctrum( states_lfp, lfp_tfa_cfg )
%
% INPUTS:
%		states_lfp  	- structure containing processed lfp data for all
%		sites of one session, output of lfp_tfa_process_lfp or
%		lfp_tfa_reject_noisy_lfp or lfp_tfa_compute_baseline_power
%       lfp_tfa_cfg     - struct containing configuration for LFP TFR analysis 
%           Required fields: see settings\lfp_tfa_settings
%               session_results_fldr            - folder to which the
%               results of the session should be saved
%               mintrials_percondition          - minimum number of trials
%               required per condition for considering the site for
%               averaging
%               analyse_epochs                  - epochs to analyse 
%               
% OUTPUTS:
%		session_pow     - output structure which saves the average LFP power spectrum for  
%       trials of a given condition for different handspace 
%       tunings and periods around the epochs analysed
%
% REQUIRES:	lfp_tfa_compare_conditions, lfp_tfa_plot_hs_tuned_psd
%
% See also lfp_tfa_process_lfp, lfp_tfa_settings,
% lfp_tfa_compare_conditions, lfp_tfa_plot_hs_tuned_psd
    
    % suppress warning for xticklabel
    warning ('off', 'MATLAB:hg:willberemoved');

    % make a folder to save figures
    results_folder_session = lfp_tfa_cfg.session_info...
        (strcmp([lfp_tfa_cfg.session_info.session], site1_lfp.session)).lfp_syncspctrm_results_fldr;
    if ~exist(results_folder_session, 'dir')
        mkdir(results_folder_session);
    end
       
    % condition based TFS
    sitepair_syncspctrm = struct();
           
        
    % folder to save sitewise results
    sitepair_results_folder = fullfile(results_folder_session, [site1_lfp.site_ID, ...
        '-', site2_lfp.site_ID]);
    if ~exist(sitepair_results_folder, 'dir')
        mkdir(sitepair_results_folder);
    end
    % struct to store condition-wise LFP power spectra average
    sitepair_syncspctrm.condition = struct();
    sitepair_syncspctrm.sites = {site1_lfp.site_ID, site2_lfp.site_ID};
    sitepair_syncspctrm.session = site1_lfp.session;
    sitepair_syncspctrm.targets = {site1_lfp.target, site2_lfp.target};
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
            cond_trials = lfp_tfa_get_condition_trials(site1_lfp, site_conditions(cn));
            % get trial indices for this condition and hand-space label
            cond_trials = cond_trials & ...
                strcmp({site1_lfp.trials.hndspc_lbl}, hs_labels(hs));
            sitepair_syncspctrm.condition(cn).ntrials(hs) = sum(cond_trials);

            fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
            fprintf('Total number of trials %g\n', sum(cond_trials));
            
            % reject noisy trials for both sites
            sitepair_syncspctrm.condition(cn).noisytrials(hs) = ...
                sum(cond_trials & ([site1_lfp.trials.noisy] | [site1_lfp.trials.noisy])); 

            % consider only non noisy trials
            fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                & ([site1_lfp.trials.noisy] | [site1_lfp.trials.noisy])));
            cond_trials = cond_trials & ~([site1_lfp.trials.noisy] | [site1_lfp.trials.noisy]);
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

                epoch_lfp = struct();
                epoch_lfp.time = {};
                epoch_lfp.trial = {};
                epoch_tfs.psd = {}; % power spectrum
                epoch_tfs.psd_f = {}; % power spectrum freq

                for t = find(cond_trials)

                    % get timing information of epoch, same for both sites
                    states          = site1_lfp.trials(t).states;
                    state_onset_t   = states([states(:).id] == ...
                        epoch_refstate).onset_t;
                    epoch_start_t   = states([states(:).id] == ...
                        epoch_refstate).onset_t + epoch_reftstart;
                    epoch_end_t     = states([states(:).id] == ...
                        epoch_refstate).onset_t + epoch_reftend;
                    % sampling frequency
                    fs = site1_lfp.trials(t).fsample; 

                    % LFP data
                    trial_lfp = [site1_lfp.trials(t).lfp_data; site2_lfp.trials(t).lfp_data];
                    trial_time = site1_lfp.trials(t).time; % same for both sites
                    trial_epoch_lfp = trial_lfp(:, trial_time >= epoch_start_t & ...
                        trial_time <= epoch_end_t);
                    trial_epoch_timestamps = trial_time(trial_time >= epoch_start_t & ...
                        trial_time <= epoch_end_t);
%                     if lfp_tfa_cfg.tfr.foi(1) < ceil(1/(epoch_reftend - epoch_reftstart))
%                         trial_epoch_lfp_zeropad = zeros(size(trial_epoch_lfp, 1), ...
%                             round(fs/lfp_tfa_cfg.tfr.foi(1)));
%                         trial_epoch_time_zeropad = zeros(1, ...
%                             round(fs/lfp_tfa_cfg.tfr.foi(1)));
%                         trial_epoch_lfp_zeropad(:,1:size(trial_epoch_lfp, 2)) = trial_epoch_lfp;
%                         trial_epoch_lfp = trial_epoch_lfp_zeropad;
%                         trial_epoch_time_zeropad(1:length(trial_epoch_timestamps)) = ...
%                             trial_epoch_timestamps;
%                         trial_epoch_time_zeropad(length(trial_epoch_timestamps) + 1:end) = ...
%                             (1/fs) : (1/fs) : (1/fs) * ...
%                             (length(trial_epoch_time_zeropad) - length(trial_epoch_timestamps));
%                         trial_epoch_time_zeropad(length(trial_epoch_timestamps) + 1:end) = ...
%                             trial_epoch_time_zeropad(length(trial_epoch_timestamps) + 1:end) + ...
%                             trial_epoch_timestamps(end);
%                         trial_epoch_timestamps = trial_epoch_time_zeropad;
%                     end                        
                    epoch_lfp.trial = [epoch_lfp.trial, ...
                        trial_epoch_lfp];  
                    epoch_lfp.time = [epoch_lfp.time, ...
                        trial_epoch_timestamps - state_onset_t];                    
                    
                end
                
                % same bins for all epochs
                ntimebins = min(cellfun('length', epoch_lfp.trial));
                for k = 1:length(epoch_lfp.trial)
                    epoch_lfp.trial{k} = epoch_lfp.trial{k}(:,1:ntimebins);
                    epoch_lfp.time{k} = epoch_lfp.time{k}(1:ntimebins);
                end
                
                epoch_lfp.fsample = fs;
                epoch_lfp.label = sitepair_syncspctrm.sites;
                
                % find the minimum frequency that can be analysed using FFT
                % based on epoch length
                min_freq = ceil(1/(epoch_reftend - epoch_reftstart));
                
                if ~isempty(epoch_lfp)
                    % calculate the cross spectrum between 2 sites for this
                    % epoch
                    cfg                     = [];
                    cfg.method              = 'mtmfft';
                    cfg.output              = 'powandcsd';
                    cfg.taper               = 'hanning';
                    cfg.foilim              = [2 120];   % frequencies of interest
                    cfg.pad                 = 0.5;
                    %cfg.t_ftimwin          = ones(length(cfg.foi),1).*0.5; %lfp_tfa_cfg.tfr.t_ftimwin;           % number of cycles per time window
                    cfg.keeptrials          = 'yes';
                    epoch_tfr               = ft_freqanalysis(cfg, epoch_lfp); 
                    
                    % calculate LFP-LFP phase sync average
                    cfg = [];
                    cfg.method              = lfp_tfa_cfg.sync.measure;
                    epoch_sync              = ft_connectivityanalysis(cfg, epoch_tfr);                    
                    
                    
                    % save LFP power spectrum
                    sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).ppc = epoch_sync;
                    sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).trials = find(cond_trials);
                    sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).hs_label = hs_labels(hs);
                    sitepair_syncspctrm.condition(cn).hs_tuned_syncsp(ep, hs).epoch_name = epoch_name;

                end

            end

        end

        % plot site average power spectrum
        if ~isempty(sitepair_syncspctrm.condition(cn).hs_tuned_syncsp)

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
        
