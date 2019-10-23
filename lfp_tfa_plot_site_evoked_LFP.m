function [ session_evoked ] = lfp_tfa_plot_site_evoked_LFP( sites_lfp, analyse_states, lfp_tfa_cfg ) 

% lfp_tfa_plot_site_evoked_LFP  - computes and plots average evoked LFP for
% different trial conditions for each site and across all sites
% of a session (Condition is a combination of
% perturbation/choice/type-effector/hand-space tuning)
%
% USAGE:
%	[ session_evoked ] = lfp_tfa_plot_site_evoked_LFP( sites_lfp, analyse_states, lfp_tfa_cfg ) 
%
% INPUTS:
%		sites_lfp     	- 1xN struct containing raw lfp data for all sites of a 
%       session, see lfp_tfa_process_lfp 
%       analyse_states  - cell array containing states to be
%       analysed and corresponding time windows
%       lfp_tfa_cfg     - struct containing configuration for TFR 
%           Required fields:
%               session_results_fldr            - folder to which the
%               results of the session should be saved
%               perturbation_groups             - 1x2 cell array containing
%               the blocks to be considered as pre- and post- injection
%               random_seed                     - for reproducibility of
%               random numbers, see rng
%               mintrials_percondition          - minimum number of trials
%               required per condition for considering the site for
%               averaging (Condition is a combination of perturbation,
%               choice, type-effector, and hand-space tuning)
%               ref_hemisphere                  - reference hemispehere for
%               ipsi- and contra-labeling
%               mintrials_percondition          - minimum number of trials
%               required per condition for considering the site for
%               averaging
%
% OUTPUTS:
%		session_evoked	- output structure which saves the condition-wise
%                       average evoked for the specified periods around the 
%                       states analysed
%           Fields:
%           sites       - 1xN struct containing condition-wise average evoked 
%                       across trials from a single site
%           session_avg - 1xT struct containing condition-wise average evoked 
%                       averaged across muliple sites in a target area in
%                       one session (T = number of target areas)
% 
% REQUIRES:	lfp_tfa_compare_conditions, lfp_tfa_get_condition_trials, 
% lfp_tfa_get_combined_evoked_lfp, lfp_tfa_get_state_evoked_lfp,
% lfp_tfa_plot_evoked_lfp
%
% See also lfp_tfa_process_lfp, lfp_tfa_compare_conditions, 
% lfp_tfa_plot_evoked_lfp, lfp_tfa_plot_site_average_tfr,
% lfp_tfa_plot_site_powspctrum
    
    % suppress warning for xticklabel
    warning ('off', 'MATLAB:hg:willberemoved');

    % make a folder to save figures
    results_folder_evoked = fullfile(lfp_tfa_cfg.session_results_fldr, 'Condition_based_Evoked_LFP');
    if ~exist(results_folder_evoked, 'dir')
        mkdir(results_folder_evoked);
    end
       
    % condition based Evoked
    sites_evoked = struct();
    session_evoked = struct();
    session_evoked.session = sites_lfp(1).session;
    % perturbation groups for this session
    perturbation_groups = lfp_tfa_cfg.perturbation_groups;
    % get trial conditions for this session
    site_conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg, perturbation_groups);
    
    % loop through each site
    for i = 1:length(sites_lfp) 

        rng(lfp_tfa_cfg.random_seed); % set random seed for reproducibility
        
        % folder to save sitewise results
        site_results_folder = fullfile(results_folder_evoked, 'sites');
        if ~exist(site_results_folder, 'dir')
            mkdir(site_results_folder);
        end
        % struct to store condition-wise evoked
        sites_evoked(i).condition = struct();
        sites_evoked(i).site_ID = sites_lfp(i).site_ID;
        sites_evoked(i).session = sites_lfp(i).session;
        sites_evoked(i).target = sites_lfp(i).target;
        % flag to indicate if this site should be used for
        % averaging based on minimum no:of trials per condition
        sites_evoked(i).use_for_avg = 1;
        
        % loop through conditions
        for cn = 1:length(site_conditions)

            % hand-space tuning of LFP
            hs_labels = site_conditions(cn).hs_labels;

            % num sites
            nsites = length(sites_lfp);                 
            
            % store details of analysed condition
            sites_evoked(i).condition(cn).label = site_conditions(cn).label;
            sites_evoked(i).condition(cn).cfg_condition = site_conditions(cn);
            sites_evoked(i).condition(cn).hs_tuned_evoked = struct(); 
            sites_evoked(i).condition(cn).ntrials = zeros(1,length(hs_labels));        

            % loop through hand space labels
            for hs = 1:length(hs_labels)
                % get trial indices for the given condition
                cond_trials = lfp_tfa_get_condition_trials(sites_lfp(i), site_conditions(cn));
                % filter trials by hand-space labels
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
                
                sites_evoked(i).condition(cn).ntrials(hs) = sum(cond_trials);

                fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
                fprintf('Total number of trials %g\n', sum(cond_trials));

                sites_evoked(i).condition(cn).noisytrials(hs) = ...
                    sum(cond_trials & [sites_lfp(i).trials.noisy]); 

                % consider only non noisy trials
                fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                    & [sites_lfp(i).trials.noisy]));
                cond_trials = cond_trials & ~[sites_lfp(i).trials.noisy];

                % check if the site contains a specified minimum number
                % of trials for all conditions
                if sum(cond_trials) < lfp_tfa_cfg.mintrials_percondition
                    sites_evoked(i).use_for_avg = 0;
                end


                % loop through time windows around the states to analyse
                for st = 1:size(analyse_states, 1)
                    
                    cond_trials_lfp = sites_lfp(i).trials(cond_trials);
                    
                    if strcmp(analyse_states{st, 1}, 'single')
                        state_tfs = lfp_tfa_get_state_evoked_lfp(cond_trials_lfp, ...
                            analyse_states(st, :));
                    elseif strcmp(analyse_states{st, 1}, 'combined')
                        state_tfs = lfp_tfa_get_combined_evoked_lfp(cond_trials_lfp, ...
                            analyse_states(st, :));
                    end                        


                    if ~isempty(state_tfs.lfp)

                        % save evoked LFP
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).lfp = state_tfs.lfp;
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean = state_tfs.mean;
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).std = state_tfs.std; 
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time = state_tfs.lfp_time;
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).dimord = 'ntrials_time';
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).trials = find(cond_trials);
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).hs_label = hs_labels(hs);
                        if isfield(state_tfs, 'state_id') && isfield(state_tfs, 'state_name')
                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state = state_tfs.state_id;
                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state_name = state_tfs.state_name;
                        end

                    end

                end

            end
            
            % plots
            % Evoked LFP
            if ~isempty(fieldnames(sites_evoked(i).condition(cn).hs_tuned_evoked))
                plottitle = ['Site ID: ', sites_evoked(i).site_ID ', Target = ' ...
                    sites_evoked(i).target '(ref_' lfp_tfa_cfg.ref_hemisphere '), '  ...
                    site_conditions(cn).label '), '];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                elseif site_conditions(cn).choice == 1
                    plottitle = [plottitle 'Choice trials'];
                end
                result_file = fullfile(site_results_folder, ...
                    ['LFP_Evoked_' sites_evoked(i).site_ID '_' site_conditions(cn).label '.png']);

                lfp_tfa_plot_evoked_lfp (sites_evoked(i).condition(cn).hs_tuned_evoked, lfp_tfa_cfg, ...
                    plottitle, result_file);
            end

        end
        site_evoked_lfp = sites_evoked(i);
        % save mat file for site
        save(fullfile(site_results_folder, ['LFP_evoked_' sites_evoked(i).site_ID '.mat']), 'site_evoked_lfp');
        % save to a mother struct
        session_evoked.sites(i) = site_evoked_lfp;
    end
        
    % Average across sites for a session
    session_avg = struct();
    % targets for this session
    targets = unique({sites_lfp.target});
    % average each target separately
    for t = 1:length(targets)
        session_avg(t).target = targets{t};
        % conditions
        for cn = 1:length(site_conditions)
            session_avg(t).condition(cn).hs_tuned_evoked = struct();
            session_avg(t).condition(cn).condition = site_conditions(cn);
            session_avg(t).condition(cn).label = site_conditions(cn).label;
            session_avg(t).condition(cn).session = sites_lfp(i).session;
            session_avg(t).condition(cn).target = sites_lfp(i).target;
            % initialize number of site pairs for each handspace
            % label
            for st = 1:size(sites_evoked(1).condition(cn).hs_tuned_evoked, 1)
                for hs = 1:size(sites_evoked(1).condition(cn).hs_tuned_evoked, 2)
                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = 0;
                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp = [];
                end
            end
            isite = 0;
            for i = 1:length(sites_evoked)
                % if the site's target is same as target being considered
                if ~strcmp(sites_lfp(i).target, targets{t})
                    continue;
                end
                if sites_evoked(i).use_for_avg
                    % calculate the average evoked LFP across sites for this condition 
                    if ~isempty(sites_evoked(i).condition(cn).hs_tuned_evoked) && ...
                        isfield(sites_evoked(i).condition(cn).hs_tuned_evoked, 'mean')
                        isite = isite + 1;
                        for hs = 1:size(sites_evoked(i).condition(cn).hs_tuned_evoked, 2)
                            for st = 1:size(sites_evoked(i).condition(cn).hs_tuned_evoked, 1)
                                if ~isempty(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).lfp)
                                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = ...
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites + 1;
                                    if session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites == 1
                                        
                                        % struct to store average evoked LFP across sites
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).hs_label = ...
                                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).hs_label;
                                        if isfield(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs), 'state') && ...
                                                isfield(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs), 'state_name')
                                            session_avg(t).condition(cn).hs_tuned_evoked(st, hs).state = ...
                                                sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state;
                                            session_avg(t).condition(cn).hs_tuned_evoked(st, hs).state_name = ...
                                                sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state_name;
                                        end
                                        
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp = ...
                                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean ;
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time = ...
                                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time;

                                    else
                                        
                                        nsamples = length(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time);
                                        if nsamples > length(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time)
                                            nsamples = length(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time);
                                        end
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp = ...
                                            cat(1, ...
                                            session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp(:, 1:nsamples), ...
                                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean(1:nsamples)) ;
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time = ...
                                            session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time(1:nsamples) ;

                                    end
                                    
                                end

                            end
                        end
                    end
                end
            end
            % average TFR across sites for a session
            if isfield(session_avg(t).condition(cn).hs_tuned_evoked, 'lfp') 
                for st = 1:size(session_avg(t).condition(cn).hs_tuned_evoked, 1)
                    for hs = 1:size(session_avg(t).condition(cn).hs_tuned_evoked, 2)                    
                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).mean = ...
                            nanmean(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp, 1);
                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = ...
                            nanstd(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp, 0, 1);
                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).dimord = 'nsites_time';
                    end
                end
            end           
            % plot average evoked LFP across sites for this session
            if ~isempty(fieldnames(session_avg(t).condition(cn).hs_tuned_evoked))
                plottitle = ['Session: ', session_avg(t).condition(cn).session ', Target = ' ...
                    session_avg(t).condition(cn).target ' (ref_' lfp_tfa_cfg.ref_hemisphere '), '  ...
                    site_conditions(cn).label ', '];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                elseif site_conditions(cn).choice == 1
                    plottitle = [plottitle 'Choice trials'];
                end
                result_file = fullfile(results_folder_evoked, ['LFP_Evoked_' ...
                    session_avg(t).condition(cn).session '_' site_conditions(cn).label '.png']);
                lfp_tfa_plot_evoked_lfp (session_avg(t).condition(cn).hs_tuned_evoked, lfp_tfa_cfg, ...
                    plottitle, result_file);
            end
        end 

    end
        
    close all;
    
    % store session average data
    session_evoked.session_avg = session_avg;
    
    % save mat files
    save(fullfile(results_folder_evoked, ['LFP_evoked_' session_evoked.session '.mat']), 'session_evoked');
    % save settings file
    save(fullfile(results_folder_evoked, 'lfp_tfa_settings.mat'), 'lfp_tfa_cfg');
end
        
