function [session_band] = lfp_tfa_plot_site_band_average( states_lfp, site_conditions, lfp_tfa_cfg )
% lfp_tfa_plot_site_band_average  - compute and plot average lfp
% averaged accross frequency bands over time
% for different hand-space tuning and trial conditions for each site and
% across all sites of a session
%
% USAGE:
%	[ session_tfs ] = lfp_tfa_plot_site_band_average( states_lfp, ...
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
results_folder_tfr = fullfile(lfp_tfa_cfg.session_results_fldr, 'Condition_based_band_average');
if ~exist(results_folder_tfr, 'dir')
    mkdir(results_folder_tfr);
end

% condition based TFS
% struct to store TFR for each site
sites_band = struct();
% struct to accumulate TFR for each site and store session average
session_band = struct();
% session name
session_band.session = states_lfp(1).session;

% Define frequency indices for band average
band_freq_index.gamma = find(lfp_tfa_cfg.band.gamma(1) < lfp_tfa_cfg.tfr.foi ...
    &  lfp_tfa_cfg.tfr.foi < lfp_tfa_cfg.band.gamma(2));
band_freq_index.beta = find(lfp_tfa_cfg.band.beta(1) < lfp_tfa_cfg.tfr.foi ...
    &  lfp_tfa_cfg.tfr.foi < lfp_tfa_cfg.band.beta(2));
band_freq_index.alpha = find(lfp_tfa_cfg.band.alpha(1) < lfp_tfa_cfg.tfr.foi ...
    &  lfp_tfa_cfg.tfr.foi < lfp_tfa_cfg.band.alpha(2));
band_freq_index.theta = find(lfp_tfa_cfg.band.theta(1) < lfp_tfa_cfg.tfr.foi ...
    &  lfp_tfa_cfg.tfr.foi < lfp_tfa_cfg.band.theta(2));



% loop through each site
for i = 1:length(states_lfp)
    
    %rng(lfp_tfa_cfg.random_seed); % set random seed for reproducibility
    
    % folder to save sitewise results
    site_results_folder = fullfile(results_folder_tfr, 'sites');
    if ~exist(site_results_folder, 'dir')
        mkdir(site_results_folder);
    end
    
    % structure to store condition-wise tfs
    sites_band(i).condition = struct();
    % info about session and site
    sites_band(i).site_ID = states_lfp(i).site_ID;
    sites_band(i).session = states_lfp(i).session;
    sites_band(i).target = states_lfp(i).target;
    % flag to indicate if this site should be used for
    % averaging based on minimum no:of trials per condition
    sites_band(i).use_for_avg = 1;
    
    % loop through each condition
    for cn = 1:length(site_conditions) % perturbation, choice, effector and so on
        
        % hand-space tuning of LFP
        hs_labels = site_conditions(cn).hs_labels;
        
        % store details of condition analysed
        sites_band(i).condition(cn).label = site_conditions(cn).label;
        sites_band(i).condition(cn).cfg_condition = site_conditions(cn);
        sites_band(i).condition(cn).hs_tuned_band = struct();
        sites_band(i).condition(cn).ntrials = zeros(1,length(hs_labels));
        
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
            
            sites_band(i).condition(cn).ntrials(hs) = sum(cond_trials);
            
            fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
            fprintf('Total number of trials %g\n', sum(cond_trials));
            
            sites_band(i).condition(cn).noisytrials(hs) = ...
                sum(cond_trials & [states_lfp(i).trials.noisy]);
            
            % consider only non noisy trials
            fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                & [states_lfp(i).trials.noisy]));
            cond_trials = cond_trials & ~[states_lfp(i).trials.noisy];
            
            % check if the site contains a specified minimum number
            % of trials for all conditions
            if sum(cond_trials) < lfp_tfa_cfg.mintrials_percondition
                sites_band(i).use_for_avg = 0;
            end
            % loop through states to analyse
            
            for st = 1:size(lfp_tfa_cfg.analyse_states, 1)
                
                if strcmp(lfp_tfa_cfg.analyse_states{st, 1}, 'combined')
                    state_band = lfp_tfa_get_combined_tfs(states_lfp(i), ...
                        cond_trials, lfp_tfa_cfg.analyse_states(st, :), lfp_tfa_cfg, site_conditions(cn).perturbation);
                end
                
                if strcmp(lfp_tfa_cfg.analyse_states{st, 1}, 'single')
                    state_band = lfp_tfa_get_state_tfs(states_lfp(i), ...
                        cond_trials, lfp_tfa_cfg.analyse_states(st, :), lfp_tfa_cfg, site_conditions(cn).perturbation);
                end
                
                if ~isempty(state_band.powspctrm)
                    
                    % save average tfs for this condition, hand-space
                    % label, and state
                    sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.powspctrm = state_band.powspctrm;
                    %                         sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.norm_mean = state_band.powspctrm_normmean;
                    %                         sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.raw_mean = state_band.powspctrm_rawmean;
                    sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.time = state_band.time;
                    sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.freq = state_band.freq;
                    sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.cfg = state_band.cfg;
                    sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.baseline = state_band.baseline;
                    if isfield(state_band, 'state_id') && isfield(state_band, 'state_name')
                        sites_band(i).condition(cn).hs_tuned_band(st, hs).state = state_band.state_id;
                        sites_band(i).condition(cn).hs_tuned_band(st, hs).state_name = state_band.state_name;
                    end
                    sites_band(i).condition(cn).hs_tuned_band(st, hs).hs_label = hs_labels(hs);
                    sites_band(i).condition(cn).hs_tuned_band(st, hs).trials = find(cond_trials);
                    sites_band(i).condition(cn).hs_tuned_band(st, hs).ntrials = length(find(cond_trials));
                    
                    sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.gamma = nanmean(state_band.powspctrm(:,band_freq_index.gamma,:),2);
                    sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.beta = nanmean(state_band.powspctrm(:,band_freq_index.beta,:),2);
                    sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.alpha = nanmean(state_band.powspctrm(:,band_freq_index.alpha,:),2);
                    sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.theta = nanmean(state_band.powspctrm(:,band_freq_index.theta,:),2);
                    
                    
                end
                
            end
            
        end
        
        
        % plot TFR
        if ~isempty(fieldnames(sites_band(i).condition(cn).hs_tuned_band))
            if site_conditions(cn).perturbation == 0
                injection = 'Pre';
            elseif site_conditions(cn).perturbation == 1
                injection = 'Post';
            end
            plottitle = ['LFP band average (' injection '): Site ' sites_band(i).site_ID ...
                ', Target ' sites_band(i).target '(ref_' lfp_tfa_cfg.ref_hemisphere '), '  ...
                site_conditions(cn).label];
            if site_conditions(cn).choice == 0
                plottitle = [plottitle 'Instructed trials'];
            elseif site_conditions(cn).choice == 1
                plottitle = [plottitle 'Choice trials'];
            end
            
%             if lfp_tfa_cfg.plot_site_average
%                 result_file = fullfile(site_results_folder, ...
%                     ['LFP_band_average' sites_band(i).site_ID '_' site_conditions(cn).label ]);
%                                 lfp_tfa_plot_hs_tuned_tfr_band_average(sites_band(i).condition(cn).hs_tuned_band, ...
%                                     lfp_tfa_cfg, plottitle, result_file);
%             end
        end
        
    end
    
    %% plot evoked LFP for pre and post injection on same plot
    for trial_type = 1:2
        plottitle = ['Site ID: ', sites_band(i).site_ID ', Target = ' ...
            sites_band(i).target '(ref_' lfp_tfa_cfg.ref_hemisphere ') '];
        if trial_type == 1
            plottitle = [plottitle 'Instructed trials'];
            trial_title = 'Instructed trials';
        elseif trial_type == 2
            plottitle = [plottitle 'Choice trials'];
            trial_title = 'Choice trials';
        end
        result_file = fullfile(site_results_folder, ...
            ['LFP_band_average ' sites_band(i).site_ID '__combined_' trial_title ]);
                 if lfp_tfa_cfg.plot_site_average
               lfp_tfa_plot_hs_tuned_tfr_band_average_combined(sites_band(i).condition, lfp_tfa_cfg, ...
                     plottitle, result_file,trial_type);
                 end
    end
    
    %%
    
    
    % save mat file for each site
    site_band = sites_band(i);
    save(fullfile(site_results_folder, ...
        ['LFP_band_average' site_band.site_ID '.mat']), 'site_band');
    % save into a mother struct
    session_band.sites(i) = site_band;
    
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
        session_avg(t).condition(cn).hs_tuned_band = struct();
        session_avg(t).condition(cn).cfg_condition = site_conditions(cn);
        session_avg(t).condition(cn).label = site_conditions(cn).label;
        session_avg(t).condition(cn).session = states_lfp(i).session;
        session_avg(t).condition(cn).target = states_lfp(i).target;
        % initialize number of site pairs for each handspace
        % label
        for st = 1:size(sites_band(1).condition(cn).hs_tuned_band, 1)
            for hs = 1:size(sites_band(1).condition(cn).hs_tuned_band, 2)
                session_avg(t).condition(cn).hs_tuned_band(st, hs).nsites = 0;
                session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.powspctrm = [];
            end
        end
        isite = 0;
        
        for i = 1:length(states_lfp)
            if strcmp(states_lfp(i).target, targets{t})
                % check if this site should be used for averaging
                if sites_band(i).use_for_avg
                    % calculate the average across sites for this condition
                    if ~isempty(sites_band(i).condition(cn).hs_tuned_band) && ...
                            isfield(sites_band(i).condition(cn).hs_tuned_band, 'freq')% && ...
                        %isite = isite + 1;
                        
                        for hs = 1:size(sites_band(i).condition(cn).hs_tuned_band, 2)
                            for st = 1:size(sites_band(i).condition(cn).hs_tuned_band, 1)
                                if isfield(sites_band(i).condition(cn).hs_tuned_band(st, hs), 'freq') ...
                                        && ~isempty(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq) ...
                                        && isfield(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq, 'powspctrm') ...
                                        && ~isempty(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.powspctrm)
                                    
                                    % increment number of site pairs
                                    session_avg(t).condition(cn).hs_tuned_band(st, hs).nsites = ...
                                        session_avg(t).condition(cn).hs_tuned_band(st, hs).nsites + 1;
                                    
                                    if session_avg(t).condition(cn).hs_tuned_band(st, hs).nsites == 1
                                        session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.powspctrm = ...
                                            nanmean(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.powspctrm, 1);
                                        session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.time = ...
                                            sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.time;
                                        session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.freq = ...
                                            sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.freq;
                                        session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.cfg = ...
                                            sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.cfg;
                                        
                                        session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.gamma = ...
                                            nanmean(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.gamma,1);
                                        session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.beta = ...
                                            nanmean(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.beta,1);
                                        session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.alpha = ...
                                            nanmean(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.alpha,1);
                                        session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.theta = ...
                                            nanmean(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.theta,1);
                                        
                                        
                                        session_avg(t).condition(cn).hs_tuned_band(st, hs).hs_label = sites_band(i).condition(cn).hs_tuned_band(st, hs).hs_label;
                                        if isfield(sites_band(i).condition(cn).hs_tuned_band(st, hs), 'state') && ...
                                                isfield(sites_band(i).condition(cn).hs_tuned_band(st, hs), 'state_name')
                                            session_avg(t).condition(cn).hs_tuned_band(st, hs).state = sites_band(i).condition(cn).hs_tuned_band(st, hs).state;
                                            session_avg(t).condition(cn).hs_tuned_band(st, hs).state_name = sites_band(i).condition(cn).hs_tuned_band(st, hs).state_name;
                                        end
                                        
                                    else
                                        ntimebins = size(session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.powspctrm, 3);
                                        % average same number of time bins
                                        if ntimebins > length(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.time)
                                            ntimebins = length(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.time);
                                            session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.powspctrm = ...
                                                cat(1, session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.powspctrm(:,:,1:ntimebins), ...
                                                nanmean(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.powspctrm(:,:,1:ntimebins), 1)) ;
                                            
                                            session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.gamma = ...
                                                cat(1,sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.gamma(:,:,1:ntimebins)),...
                                                nanmean(session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.gamma(:,:,1:ntimebins),1);
                                            session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.beta = ...
                                                cat(1,ssession_avg(t).condition(cn).hs_tuned_band(st, hs).freq.beta(:,:,1:ntimebins)),...
                                                nanmean(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.beta(:,:,1:ntimebins),1);
                                            session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.alpha = ...
                                                cat(1,session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.alpha(:,:,1:ntimebins)),...
                                                nanmean(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.alpha(:,:,1:ntimebins),1);
                                            session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.theta = ...
                                                cat(1,session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.theta(:,:,1:ntimebins)),...
                                                nanmean(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.theta(:,:,1:ntimebins),1);
                                            
                                            
                                        else
                                            if ~isempty(session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.powspctrm)
                                                session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.powspctrm = ...
                                                    cat(1, session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.powspctrm, ...
                                                    nanmean(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.powspctrm(:,:,1:ntimebins), 1));
                                                session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.gamma = ...
                                                    cat(1,session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.gamma,...
                                                    nanmean(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.gamma(:,:,1:ntimebins),1));
                                                session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.beta = ...
                                                    cat(1,session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.beta,...
                                                    nanmean(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.beta(:,:,1:ntimebins),1));
                                                session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.alpha = ...
                                                    cat(1,session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.alpha,...
                                                    nanmean(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.alpha(:,:,1:ntimebins),1));
                                                session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.theta = ...
                                                    cat(1,session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.theta,...
                                                    nanmean(sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.theta(:,:,1:ntimebins),1));
                                                
                                            else
                                                session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.powspctrm = ...
                                                    sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.powspctrm(:,:,1:ntimebins) ;
                                                
                                                session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.gamma = ...
                                                    sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.gamma(:,:,1:ntimebins);
                                                session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.beta = ...
                                                    sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.beta(:,:,1:ntimebins);
                                                session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.alpha = ...
                                                    sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.alpha(:,:,1:ntimebins);
                                                session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.theta = ...
                                                    sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.theta(:,:,1:ntimebins);
                                                
                                                
                                            end
                                        end
                                        session_avg(t).condition(cn).hs_tuned_band(st, hs).freq.time = ...
                                            sites_band(i).condition(cn).hs_tuned_band(st, hs).freq.time(1:ntimebins);
                                    end
                                    
                                    
                                end
                            end
                        end
                    end
                end
            else
                continue;
            end
        end
        
        
        % plot average TFR for this condition and target
        
        if ~isempty(session_avg(t).condition(cn).hs_tuned_band)
            if isfield(session_avg(t).condition(cn).hs_tuned_band, 'freq')
                plottitle = ['LFP band average: Target = ' session_avg(t).target ...
                    ', (ref_', lfp_tfa_cfg.ref_hemisphere, ') ',  ...
                    'Session ', session_avg(t).condition(cn).session, ...
                    ' ', site_conditions(cn).label];
                
                result_file = fullfile(results_folder_tfr, ...
                    ['LFP_band_average' session_avg(t).target '_'...
                    session_avg(t).condition(cn).session '_' site_conditions(cn).label ]);
                 lfp_tfa_plot_hs_tuned_tfr_band_average(session_avg(t).condition(cn).hs_tuned_band, ...
                     lfp_tfa_cfg, plottitle, result_file);
            end
        end
        
    end
    
    for trial_type = 1:2
        plottitle = ['LFP band average: Target = ' session_avg(t).target ...
            ', (ref_', lfp_tfa_cfg.ref_hemisphere, ') ',  ...
            'Session ', session_avg(t).condition(cn).session, ...
            ' ', site_conditions(cn).label];
        if trial_type == 1
            plottitle = [plottitle 'Instructed'];
            trial_title = 'Instructed';
        elseif trial_type == 2
            plottitle = [plottitle 'Choice'];
            trial_title = 'Choice';
        end
        result_file = fullfile(results_folder_tfr, ...
            ['LFP_band_average, ' session_avg(t).target ...
                    ', (ref_', lfp_tfa_cfg.ref_hemisphere, ') ',  ...
                    'Session ', session_avg(t).condition(cn).session, ...
                    ' ', '__comb_' trial_title]);
        
            lfp_tfa_plot_hs_tuned_tfr_band_average_combined(session_avg(t).condition, lfp_tfa_cfg, ...
                plottitle, result_file,trial_type);
      
        
    end
    
    session_band.session_avg = session_avg;
    
    % close figures
    close all;
    
    % save session average tfs
    save(fullfile(results_folder_tfr, ['LFP_band_' session_band.session '.mat']), 'session_band');
    % save settings file
    save(fullfile(results_folder_tfr, 'lfp_tfa_settings.mat'), 'lfp_tfa_cfg');
    
end
