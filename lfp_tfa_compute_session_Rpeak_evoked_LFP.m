function [ session_evoked ] = lfp_tfa_compute_session_Rpeak_evoked_LFP( session_proc_lfp, analyse_states, lfp_tfa_cfg ) 

% lfp_tfa_plot_average_evoked_LFP  - plots average evoked LFP for
% different hand-space tuning conditions for each site and across all sites
% of a session
%
% USAGE:
%	[ session_evoked ] = lfp_tfa_plot_site_evoked_LFP( states_lfp, analyse_states, lfp_tfa_cfg ) 
%
% INPUTS:
%		states_lfp  	- struct containing raw lfp data for all sites of a 
%       session, output from lfp_tfa_process_lfp or
%       lfp_tfa_compute_baseline or lfp_tfa_reject_noisy_lfp
%       analyse_states  - cell array containing states to be
%       analysed and corresponding time windows
%       lfp_tfa_cfg     - struct containing configuration for TFR 
%           Required fields:
%               session_results_fldr            - folder to which the
%               results of the session should be saved
%               mintrials_percondition          - minimum number of trials
%               required per condition for considering the site for
%               averaging
%               analyse_states                  - states to analyse 
%
% OUTPUTS:
%		session_evoked	- output structure which saves the average evoked LFP for  
%                         trials of a given condition for different handspace 
%                         tunings and periods around the states analysed
% 
% REQUIRES:	lfp_tfa_compare_conditions, lfp_tfa_plot_evoked_lfp
%
% See also lfp_tfa_process_lfp, lfp_tfa_compute_baseline, lfp_tfa_reject_noisy_lfp, 
% lfp_tfa_compare_conditions, lfp_tfa_plot_evoked_lfp
    
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
    session_evoked.session = session_proc_lfp(1).session;
    % perturbation groups for this session
    %perturbation_groups = lfp_tfa_cfg.perturbation_groups;
    % get trial conditions for this session
    site_conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg, {0, 1});
    
    % loop through each site
    for i = 1:length(session_proc_lfp) 

        rng(lfp_tfa_cfg.random_seed); % set random seed for reproducibility
        
        % folder to save sitewise results
        site_results_folder = fullfile(results_folder_evoked, 'sites');
        if ~exist(site_results_folder, 'dir')
            mkdir(site_results_folder);
        end
        % struct to store condition-wise evoked
        sites_evoked(i).condition = struct();
        sites_evoked(i).site_ID = session_proc_lfp(i).site_ID;
        sites_evoked(i).session = session_proc_lfp(i).session;
        sites_evoked(i).target = session_proc_lfp(i).target;
        % flag to indicate if this site should be used for
        % averaging based on minimum no:of trials per condition
        sites_evoked(i).use_for_avg = 1;
        
        % loop through conditions
        for cn = 1:length(site_conditions)

            % hand-space tuning of LFP
            hs_labels = site_conditions(cn).hs_labels;

            % num sites
            nsites = length(session_proc_lfp);                 
            
            % store details of analysed condition
            sites_evoked(i).condition(cn).label = site_conditions(cn).label;
            sites_evoked(i).condition(cn).cfg_condition = site_conditions(cn);
            sites_evoked(i).condition(cn).hs_tuned_evoked = struct(); 
            sites_evoked(i).condition(cn).ntrials = zeros(1,length(hs_labels));        

            % loop through hand space labels
            for hs = 1:length(hs_labels)
                % get trial indices for the given condition
                cond_trials = lfp_tfa_get_condition_trials(session_proc_lfp(i), site_conditions(cn));
                % filter trials by hand-space labels
                if ~strcmp(site_conditions(cn).reach_hands{hs}, 'any')
                    cond_trials = cond_trials & ...
                        strcmp({session_proc_lfp(i).trials.reach_hand}, ...
                        site_conditions(cn).reach_hands{hs});
                end
                if ~strcmp(site_conditions(cn).reach_spaces{hs}, 'any')
                    cond_trials = cond_trials & ...
                        strcmp({session_proc_lfp(i).trials.reach_space}, ...
                        site_conditions(cn).reach_spaces{hs});
                end
                
                sites_evoked(i).condition(cn).ntrials(hs) = sum(cond_trials);

                fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
                fprintf('Total number of trials %g\n', sum(cond_trials));

                sites_evoked(i).condition(cn).noisytrials(hs) = ...
                    sum(cond_trials & [session_proc_lfp(i).trials.noisy]); 

                % consider only non noisy trials
                fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                    & [session_proc_lfp(i).trials.noisy]));
                cond_trials = cond_trials & ~[session_proc_lfp(i).trials.noisy];

                % check if the site contains a specified minimum number
                % of trials for all conditions
                if sum(cond_trials) < lfp_tfa_cfg.mintrials_percondition
                    sites_evoked(i).use_for_avg = 0;
                end
                
                if sum(cond_trials) == 0
                    continue;
                end


                % loop through time windows around the states to analyse
                for st = 1:size(analyse_states, 1)
                    
                    cond_trials_lfp = session_proc_lfp(i).trials(cond_trials);
                    
                    if strcmp(analyse_states{st, 1}, 'ecg')
                        state_evoked = lfp_tfa_get_ECG_triggered_evoked(cond_trials_lfp, ...
                            analyse_states(st, :));
                        if isfield(lfp_tfa_cfg, 'random_permute_triggers') && ...
                                lfp_tfa_cfg.random_permute_triggers
                            nshuffles = 100;
                            if isfield(lfp_tfa_cfg, 'n_shuffles')
                                nshuffles = lfp_tfa_cfg.n_shuffles;
                            end
                            shuffled_Rpeak_evoked = lfp_tfa_shuffled_Rpeak_evoked_LFP(cond_trials_lfp, ...
                                analyse_states(st, :), nshuffles);
                        end
                        
%                         state_evoked = lfp_tfa_get_ECG_based_STA(cond_trials_lfp, ...
%                             site_lfp(i).site_ID, analyse_states(st, :), 'lfp');
%                         state_evoked_ecg = lfp_tfa_get_ECG_based_STA(cond_trials_lfp, ...
%                             site_lfp(i).site_ID, analyse_states(st, :), 'ecg');
%                     else 
%                         state_evoked = lfp_tfa_get_state_evoked_lfp(cond_trials_lfp, ...
%                             analyse_states(st, :));                    
                    end                        


                    if ~isempty(state_evoked.lfp)

                        % save evoked LFP
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).lfp = state_evoked.lfp;
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean = state_evoked.mean;
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).std = state_evoked.std;                          
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time = state_evoked.lfp_time;
                        
                        if isfield(lfp_tfa_cfg, 'random_permute_triggers') && ...
                                lfp_tfa_cfg.random_permute_triggers
                            
                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).shuffled_lfp = ...
                                shuffled_Rpeak_evoked.lfp;
                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).shuffled_mean = ...
                                shuffled_Rpeak_evoked.mean;
                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).shuffled_std = ...
                                shuffled_Rpeak_evoked.std;
                            
                        end
                        
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).trials = find(cond_trials);
                        sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).hs_label = hs_labels(hs);
                        if isfield(state_evoked, 'state_id') && isfield(state_evoked, 'state_name')
                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state = state_evoked.state_id;
                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state_name = state_evoked.state_name;
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
        
        % difference between conditions
        sites_evoked(i).difference = [];
        for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
            diff_condition = lfp_tfa_cfg.diff_condition{diff};
            diff_color = []; diff_legend = [];
            if isfield(lfp_tfa_cfg, 'diff_color')
                diff_color = lfp_tfa_cfg.diff_color{diff};
            end
            if isfield(lfp_tfa_cfg, 'diff_legend')
                diff_legend = lfp_tfa_cfg.diff_legend{diff};
            end
            sites_evoked(i).difference = [sites_evoked(i).difference, ...
                lfp_tfa_compute_diff_condition_evoked(sites_evoked(i).condition, diff_condition, diff_color, diff_legend)];
        end
        % plot Difference TFR
        for dcn = 1:length(sites_evoked(i).difference)
            if ~isempty(sites_evoked(i).difference(dcn).hs_tuned_evoked)
                if isfield(sites_evoked(i).difference(dcn).hs_tuned_evoked,... 
                        'mean')
                    plottitle = ['Site ID: ' sites_evoked(i).site_ID ', Target = ' ...
                        sites_evoked(i).target '(ref_' lfp_tfa_cfg.ref_hemisphere '), ' ...
                        sites_evoked(i).difference(dcn).label];
                    result_file = fullfile(site_results_folder, ...
                        ['LFP_DiffEvoked_' sites_evoked(i).site_ID ...
                        '_' 'diff_condition' num2str(dcn) '.png']);
                        %sites_avg(t).difference(dcn).label '.png']);
                    lfp_tfa_plot_evoked_lfp(sites_evoked(i).difference(dcn).hs_tuned_evoked, ...
                        lfp_tfa_cfg, plottitle, result_file);
                end
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
    targets = unique({session_proc_lfp.target});
    % average each target separately
    for t = 1:length(targets)
        session_avg(t).target = targets{t};
        
        % conditions
        for cn = 1:length(site_conditions)
            session_avg(t).condition(cn).hs_tuned_evoked = struct();
            isite = 0;
            session_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
            session_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
            for st = 1:size(sites_evoked(1).condition(1).hs_tuned_evoked, 1)
                for hs = 1:size(sites_evoked(1).condition(1).hs_tuned_evoked, 2)
                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = 0;
                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp = [];
                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).mean = [];
                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = [];
                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time = [];
                end
            end
            for i = 1:length(sites_evoked)
                % if the site's target is same as target being considered
                if ~strcmp(session_proc_lfp(i).target, targets{t})
                    continue;
                end
                if sites_evoked(i).use_for_avg
                    % calculate the average evoked LFP across sites for this condition 
                    if ~isempty(sites_evoked(i).condition(cn).hs_tuned_evoked) && ...
                        isfield(sites_evoked(i).condition(cn).hs_tuned_evoked, 'mean')
                        isite = isite + 1;
                        for hs = 1:size(sites_evoked(i).condition(cn).hs_tuned_evoked, 2)
                            for st = 1:size(sites_evoked(i).condition(cn).hs_tuned_evoked, 1)
                                if ~isempty(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean)
                                    session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = ...
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites + 1;
                                    if session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites == 1
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp = ...
                                            [session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp; ...
                                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean] ;
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = ...
%                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).std ;
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time = ...
                                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time;

                                    else
                                        nsamples = length(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time);
                                        if nsamples > length(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time)
                                            nsamples = length(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time);
                                        end
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp = ...
                                            [session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp(:,1:nsamples); ...
                                            sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean(1:nsamples)] ;
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = ...
%                                             session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std(1:nsamples) + ...
%                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).std(1:nsamples) ;
                                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time = ...
                                            session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time(1:nsamples) ;

                                    end
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
                                    session_avg(t).condition(cn).condition = site_conditions(cn);
                                    session_avg(t).condition(cn).label = site_conditions(cn).label;
                                    session_avg(t).condition(cn).session = session_proc_lfp(i).session;
                                    session_avg(t).condition(cn).target = session_proc_lfp(i).target;
                                    %session_avg(t).condition(cn).nsites = nsites;
                                end

                            end
                        end
                    end
                end
            end
            % average TFR across sites for a session
            if isfield(session_avg(t).condition(cn).hs_tuned_evoked, 'lfp') && ...
                    ~isempty([session_avg(t).condition(cn).hs_tuned_evoked.lfp])
                
                for hs = 1:size(session_avg(t).condition(cn).hs_tuned_evoked, 2)
                    for st = 1:size(session_avg(t).condition(cn).hs_tuned_evoked, 1)
                        if ~isempty(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp)
                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).mean = ...
                            nanmean(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp, 1);
                        session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = ...
                            nanstd(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp, 0, 1);
                        %session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = isite;
                        end
                    end
                end
            end           
            % plot average evoked LFP across sites for this session
            if ~isempty([session_avg(t).condition(cn).hs_tuned_evoked.lfp])
                plottitle = ['Session: ', session_avg(t).condition(cn).session ', Target = ' ...
                    session_avg(t).condition(cn).target ' (ref_' lfp_tfa_cfg.ref_hemisphere '), '  ...
                    site_conditions(cn).label ', '];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                elseif site_conditions(cn).choice == 1
                    plottitle = [plottitle 'Choice trials'];
                end
                result_file = fullfile(results_folder_evoked, ['LFP_Evoked_' ...
                    session_avg(t).condition(cn).session '_' ...
                    session_avg(t).condition(cn).target '_' site_conditions(cn).label '.png']);
                lfp_tfa_plot_evoked_lfp (session_avg(t).condition(cn).hs_tuned_evoked, lfp_tfa_cfg, ...
                    plottitle, result_file);
            end
        end 
        
        
        % difference between conditions
        session_avg(t).difference = [];
        for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
            diff_condition = lfp_tfa_cfg.diff_condition{diff};
            diff_color = []; diff_legend = [];
            if isfield(lfp_tfa_cfg, 'diff_color')
                diff_color = lfp_tfa_cfg.diff_color{diff};
            end
            if isfield(lfp_tfa_cfg, 'diff_legend')
                diff_legend = lfp_tfa_cfg.diff_legend{diff};
            end
            session_avg(t).difference = [session_avg(t).difference, ...
                lfp_tfa_compute_diff_condition_evoked(session_avg(t).condition, diff_condition, diff_color, diff_legend)];
        end
        % plot Difference TFR
        for dcn = 1:length(session_avg(t).difference)
            if ~isempty(session_avg(t).difference(dcn).hs_tuned_evoked)
                if isfield(session_avg(t).difference(dcn).hs_tuned_evoked,... 
                        'mean')
                    plottitle = ['Target ', lfp_tfa_cfg.compare.targets{t}, ...
                    ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') ', ...
                    session_avg(t).difference(dcn).label];
                    result_file = fullfile(results_folder_evoked, ...
                        ['LFP_DiffEvoked_' lfp_tfa_cfg.compare.targets{t} ...
                        '_' 'diff_condition' num2str(dcn) '.png']);
                        %session_avg(t).difference(dcn).label '.png']);
                    lfp_tfa_plot_evoked_lfp(session_avg(t).difference(dcn).hs_tuned_evoked, ...
                        lfp_tfa_cfg, plottitle, result_file);
                end
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
        
