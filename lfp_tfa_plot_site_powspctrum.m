function [ session_pow ] = lfp_tfa_plot_site_powspctrum( states_lfp, lfp_tfa_cfg ) 

% lfp_tfa_plot_site_powspctrum  - calculate and plot the average lfp power spectrum for
% different trial conditions for each site and across sites for a session
% (Condition is a combination of
% perturbation/choice/type-effector/hand-space tuning)
%
% USAGE:
%	[ session_pow ] = lfp_tfa_plot_site_powspctrum( states_lfp, lfp_tfa_cfg )
%
% INPUTS:
%		states_lfp  	- 1xN structure containing processed lfp data for all
%		sites of one session (N=number of sites), see lfp_tfa_process_lfp
%       lfp_tfa_cfg     - struct containing configuration for LFP TFR analysis 
%           Required fields: see settings\lfp_tfa_settings_example
%               session_results_fldr            - folder to which the
%               results of the session should be saved
%               perturbation_groups             - 1x2 cell array containing
%               the blocks to be considered as pre- and post- injection
%               mintrials_percondition          - minimum number of trials
%               required per condition for considering the site for
%               averaging
%               analyse_epochs                  - epochs to analyse 
%               ref_hemisphere                  - reference hemispehere for
%               ipsi- and contra-labeling
%               
% OUTPUTS:
%		session_pow	- output structure which saves the condition-wise
%                   average LFP power spectra for the specified epochs
%           Fields:
%           sites       - 1xN struct containing condition-wise average LFP spectra 
%                       across trials from a single site
%           session_avg - 1xT struct containing condition-wise average LFP spectra 
%                       averaged across muliple sites in a target area in
%                       one session (T = number of target areas)
%
% REQUIRES:	lfp_tfa_compare_conditions, lfp_tfa_get_condition_trials,
% lfp_tfa_plot_hs_tuned_psd_2 
%
% See also lfp_tfa_process_lfp, settings/lfp_tfa_settings_example,
% lfp_tfa_compare_conditions, lfp_tfa_plot_hs_tuned_psd_2
    
    % suppress warning for xticklabel
    warning ('off', 'MATLAB:hg:willberemoved');

    % make a folder to save figures
    results_folder_psd = fullfile(lfp_tfa_cfg.session_results_fldr, 'Condition_based_LFP_Power');
    if ~exist(results_folder_psd, 'dir')
        mkdir(results_folder_psd);
    end
       
    % condition based TFS
    sites_pow = struct();
    session_pow = struct();
    session_pow.session = states_lfp(1).session;
    % perturbation groups for this session
    perturbation_groups = lfp_tfa_cfg.perturbation_groups;
    site_conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg, perturbation_groups);
    
    % loop through each site
    for i = 1:length(states_lfp)  
        
        
        % folder to save sitewise results
%         site_results_folder = fullfile(results_folder_psd, states_lfp(i).site_ID);
%         if ~exist(site_results_folder, 'dir')
%             mkdir(site_results_folder);
%         end
        % struct to store condition-wise LFP power spectra average
        sites_pow(i).condition = struct();
        sites_pow(i).site_ID = states_lfp(i).site_ID;
        sites_pow(i).session = states_lfp(i).session;
        sites_pow(i).target = states_lfp(i).target;
        % flag to indicate if this site should be used for
        % averaging based on minimum no:of trials per condition
        sites_pow(i).use_for_avg = 1;
        
        % loop through conditions
        for cn = 1:length(site_conditions)

            % hand-space tuning of LFP
            hs_labels = site_conditions(cn).hs_labels;

            % store details of condition
            sites_pow(i).condition(cn).label = site_conditions(cn).label;
            sites_pow(i).condition(cn).cfg_condition = site_conditions(cn);
            sites_pow(i).condition(cn).hs_tuned_tfs = struct(); 
            sites_pow(i).condition(cn).ntrials = zeros(1,length(hs_labels));           

            for hs = 1:length(hs_labels)
                % get trial indices for this condition
                cond_trials = lfp_tfa_get_condition_trials(states_lfp(i), site_conditions(cn));
                % get trial indices for this condition and hand-space label
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
                
                sites_pow(i).condition(cn).ntrials(hs) = sum(cond_trials);

                fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
                fprintf('Total number of trials %g\n', sum(cond_trials));

                sites_pow(i).condition(cn).noisytrials(hs) = ...
                    sum(cond_trials & [states_lfp(i).trials.noisy]); 

                % consider only non noisy trials
                fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                    & [states_lfp(i).trials.noisy]));
                cond_trials = cond_trials & ~[states_lfp(i).trials.noisy];
                % check if the site contains a specified minimum number
                % of trials for all conditions
                if sum(cond_trials) < lfp_tfa_cfg.mintrials_percondition
                    sites_pow(i).use_for_avg = 0;
                end


                
                % loop through epochs to analyse
                for ep = 1:size(lfp_tfa_cfg.analyse_epochs, 1)
                    epoch_refstate   = lfp_tfa_cfg.analyse_epochs{ep, 1};
                    epoch_name       = lfp_tfa_cfg.analyse_epochs{ep, 2};
                    epoch_reftstart  = lfp_tfa_cfg.analyse_epochs{ep, 3};
                    epoch_reftend    = lfp_tfa_cfg.analyse_epochs{ep, 4};
                    
                    epoch_tfs.psd = {}; % power spectrum
                    epoch_tfs.psd_f = {}; % power spectrum freq

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
                        epoch_powspctrm = states_lfp(i).trials(t).tfs.powspctrm(1, :, ...
                            (states_lfp(i).trials(t).tfs.time >= epoch_start_t & ...
                            states_lfp(i).trials(t).tfs.time <= epoch_end_t));  
                        epoch_tfs.psd = [epoch_tfs.psd sum(epoch_powspctrm, 3)];
                        epoch_tfs.psd_f = states_lfp(i).trials(t).tfs.freq;

                    end

                    if ~isempty(epoch_tfs.psd)
                        % power spectrum average
                        arr_state_psd = vertcat(epoch_tfs.psd{:});
                        epoch_psd_mean = nanmean(arr_state_psd, 1);

                        % save LFP power spectrum
                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).psd = arr_state_psd;
                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).dimord = 'ntrials_freq';
                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).mean = epoch_psd_mean;
                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).freq = epoch_tfs.psd_f;
                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).trials = find(cond_trials);
                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).hs_label = hs_labels(hs);
                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).epoch_name = epoch_name;

                    end

                end

            end
            
            % plot site average power spectrum
            if ~isempty(fieldnames(sites_pow(i).condition(cn).hs_tuned_power))

                plottitle = ['Site ID: ', sites_pow(i).site_ID ...
                    ', Target = ' sites_pow(i).target ' (ref_' lfp_tfa_cfg.ref_hemisphere '), ' ...
                'Perturb ' num2str(site_conditions(cn).perturbation_group{1}) ', '];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                else
                    plottitle = [plottitle 'Choice trials'];
                end
                site_results_folder = fullfile(results_folder_psd, 'sites', sites_pow(i).site_ID);
                if ~exist(site_results_folder, 'dir')
                    mkdir(site_results_folder);
                end
                result_file = fullfile(site_results_folder, ...
                    ['LFP_Power_' states_lfp(i).site_ID '_' site_conditions(cn).label '.png']);
                lfp_tfa_plot_hs_tuned_psd_2(sites_pow(i).condition(cn).hs_tuned_power, ...
                    lfp_tfa_cfg, plottitle, result_file);
            end

        end
        % save mat file for site
        site_pow = sites_pow(i);
        save(fullfile(site_results_folder, ...
            ['LFP_Power_' site_pow.site_ID '.mat']), 'site_pow');
        % save into a mother struct
        session_pow.sites(i) = site_pow;
    end
        
    
    
    % Calculate average power spectrum across all sites
    session_avg = struct();
    % targets for this session
    targets = unique({states_lfp.target});
    for t = 1:length(targets)
        session_avg(t).target = targets{t};
        for cn = 1:length(site_conditions)
            session_avg(t).condition(cn).hs_tuned_power = struct();
            session_avg(t).condition(cn).condition = site_conditions(cn);
            session_avg(t).condition(cn).session = states_lfp(i).session;
            session_avg(t).condition(cn).target = states_lfp(i).target;
            session_avg(t).condition(cn).label = site_conditions(cn).label;% variable to store no:of sites with trials satisfying this
            % condition
            % initialize number of site pairs for each handspace
            % label
            for st = 1:size(sites_pow(1).condition(cn).hs_tuned_power, 1)
                for hs = 1:size(sites_pow(1).condition(cn).hs_tuned_power, 2)
                    session_avg(t).condition(cn).hs_tuned_power(st, hs).nsites = 0;
                    session_avg(t).condition(cn).hs_tuned_power(st, hs).psd = [];
                end
            end
            isite = 0;
            for i = 1:length(states_lfp)
                % if site's target is the target beong considered
                if ~strcmp(states_lfp(i).target, targets{t})
                    continue;
                end
                if sites_pow(i).use_for_avg                    
                % calculate the average LFP power spectrum across sites for this condition 
                    if ~isempty(sites_pow(i).condition(cn).hs_tuned_power) && ...
                            isfield(sites_pow(i).condition(cn).hs_tuned_power, 'mean')
                        %isite = isite + 1;


                        % hand-space labels
                        for hs = 1:size(sites_pow(i).condition(cn).hs_tuned_power, 2)
                            % epochs
                            for ep = 1:size(sites_pow(i).condition(cn).hs_tuned_power, 1)
                                if ~isempty(sites_pow(i).condition(cn).hs_tuned_power(ep, hs).lfp)
                                    session_avg(t).condition(cn).hs_tuned_power(st, hs).nsites = ...
                                        session_avg(t).condition(cn).hs_tuned_power(st, hs).nsites + 1;

                                    if session_avg(t).condition(cn).hs_tuned_power(st, hs).nsites == 1
                                        session_avg(t).condition(cn).hs_tuned_power(ep, hs).psd = ...
                                            sites_pow(i).condition(cn).hs_tuned_power(ep, hs).mean ;
                                        session_avg(t).condition(cn).hs_tuned_power(ep, hs).freq = ...
                                            sites_pow(i).condition(cn).hs_tuned_power(ep, hs).freq;
                                    else
                                        if ~isempty(session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean)
                                            nfreqs = length(session_avg(t).condition(cn).hs_tuned_power(ep, hs).freq);
                                            % average same number fo frequency
                                            % bins
                                            if nfreqs > length(sites_pow(i).condition(cn).hs_tuned_power(ep, hs).freq)
                                                nfreqs = length(sites_pow(i).condition(cn).hs_tuned_power(ep, hs).freq);
                                            end                               
                                            session_avg(t).condition(cn).hs_tuned_power(ep, hs).psd = ...
                                                [session_avg(t).condition(cn).hs_tuned_power(ep, hs).psd(1:nfreqs); ...
                                                sites_pow(i).condition(cn).hs_tuned_power(ep, hs).mean(1:nfreqs)];
                                            session_avg(t).condition(cn).hs_tuned_power(ep, hs).psd = ...
                                                sites_pow(i).condition(cn).hs_tuned_power(ep, hs).freq(1:nfreqs);
                                        else
                                            session_avg(t).condition(cn).hs_tuned_power(ep, hs).psd = ...
                                                sites_pow(i).condition(cn).hs_tuned_power(ep, hs).mean ;
                                        end
                                    end
                                    % store lfp power spectra average for
                                    % session
                                    session_avg(t).condition(cn).hs_tuned_power(ep, hs).hs_label = ...
                                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).hs_label;
                                    session_avg(t).condition(cn).hs_tuned_power(ep, hs).epoch_name = ...
                                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).epoch_name;
                                    
                                end
                            end
                        end
                    end
                end
            end
            
            % average TFR across sites for a session
            if isfield(session_avg(t).condition(cn).hs_tuned_power, 'mean')
                for hs = 1:size(session_avg(t).condition(cn).hs_tuned_power, 2)
                    for ep = 1:size(session_avg(t).condition(cn).hs_tuned_power, 1)
                        session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean = ...
                            nanmean(session_avg(t).condition(cn).hs_tuned_power(ep, hs).psd);
                        session_avg(t).condition(cn).hs_tuned_power(ep, hs).dimord = 'nsites_freq';
                    end
                end
            end 
            
            % plot average power spectrum across sites for this session
            if ~isempty(fieldnames(session_avg(t).condition(cn).hs_tuned_power))
                plottitle = ['Session: ', session_avg(t).condition(cn).session ...
                    ', Target = ' session_avg(t).condition(cn).target '(ref_' lfp_tfa_cfg.ref_hemisphere '), '  ...
                    'Perturb ' num2str(site_conditions(cn).perturbation_group{1}) ', '];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                else
                    plottitle = [plottitle 'Choice trials'];
                end
                results_file = fullfile(results_folder_psd, ...
                    ['LFP_Power_' session_avg(t).condition(cn).session '_' site_conditions(cn).label '.png']);
                lfp_tfa_plot_hs_tuned_psd_2(session_avg(t).condition(cn).hs_tuned_power, ...
                            lfp_tfa_cfg, plottitle, results_file);
            end
            
        end        
    end
    
    close all;
    % store session average
    session_pow.session_avg = session_avg;
    
    % save mat files
    save(fullfile(results_folder_psd, ['LFP_Power_' session_pow.session '.mat']), 'session_pow');
    % save settings file
    save(fullfile(results_folder_psd, 'lfp_tfa_settings.mat'), 'lfp_tfa_cfg');
end
        
