function [ session_ecg_tfs ] = lfp_tfa_plot_session_tfs_ECG( session_ecg, session_info, analyse_states, lfp_tfa_cfg ) 

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
    results_folder_tfs = fullfile(session_info.analyse_lfp_fldr, 'ECG_analysis');
    if ~exist(results_folder_tfs, 'dir')
        mkdir(results_folder_tfs);
    end
       
    % condition based Evoked
    session_ecg_tfs = struct();
    %session_ecg_tfs = struct();
    %session_ecg_tfs.session = session_ecg(1).session;
    
    % get trial conditions for this session
    conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg, {0, 1});
    
    % loop through each site
    for i = 1:length(session_ecg) 

        rng(lfp_tfa_cfg.random_seed); % set random seed for reproducibility
        
        % struct to store condition-wise evoked
        session_ecg_tfs(i).condition = struct();
        session_ecg_tfs(i).session = session_ecg(i).session;
        % flag to indicate if this site should be used for
        % averaging based on minimum no:of trials per condition
        session_ecg_tfs(i).use_for_avg = 1;
        
        % loop through conditions
        for cn = 1:length(conditions)

            % hand-space tuning of LFP
            hs_labels = conditions(cn).hs_labels;

            % store details of analysed condition
            session_ecg_tfs(i).condition(cn).label = conditions(cn).label;
            session_ecg_tfs(i).condition(cn).cfg_condition = conditions(cn);
            session_ecg_tfs(i).condition(cn).hs_tuned_tfs = struct(); 
            session_ecg_tfs(i).condition(cn).ntrials = zeros(1,length(hs_labels));        

            % loop through hand space labels
            for hs = 1:length(hs_labels)
                % get trial indices for the given condition
                cond_trials = lfp_tfa_get_condition_trials(session_ecg(i), conditions(cn));
                % filter trials by hand-space labels
                if ~strcmp(conditions(cn).reach_hands{hs}, 'any')
                    cond_trials = cond_trials & ...
                        strcmp({session_ecg(i).trials.reach_hand}, ...
                        conditions(cn).reach_hands{hs});
                end
                if ~strcmp(conditions(cn).reach_spaces{hs}, 'any')
                    cond_trials = cond_trials & ...
                        strcmp({session_ecg(i).trials.reach_space}, ...
                        conditions(cn).reach_spaces{hs});
                end
                
                session_ecg_tfs(i).condition(cn).ntrials(hs) = sum(cond_trials);

                fprintf('Condition %s - %s\n', conditions(cn).label, hs_labels{hs});
                fprintf('Total number of trials %g\n', sum(cond_trials));

                session_ecg_tfs(i).condition(cn).noisytrials(hs) = ...
                    sum(cond_trials & [session_ecg(i).trials.noisy]); 

                % consider only non noisy trials
                fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                    & [session_ecg(i).trials.noisy]));
                cond_trials = cond_trials & ~[session_ecg(i).trials.noisy];

                % check if the site contains a specified minimum number
                % of trials for all conditions
                if sum(cond_trials) < lfp_tfa_cfg.mintrials_percondition
                    session_ecg_tfs(i).use_for_avg = 0;
                end


                % loop through time windows around the states to analyse
                for st = 1:size(analyse_states, 1)
                             
                    state_tfs = lfp_tfa_get_state_tfs_ECG(session_ecg, ...
                        cond_trials, analyse_states(st, :), lfp_tfa_cfg);                    


                    if ~isempty(state_tfs.powspctrm)

                        % save tfs ECG
                        session_ecg_tfs(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm = state_tfs.powspctrm_normmean;
                        session_ecg_tfs(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm_raw = state_tfs.powspctrm;
                        session_ecg_tfs(i).condition(cn).hs_tuned_tfs(st, hs).time = state_tfs.time;
                        session_ecg_tfs(i).condition(cn).hs_tuned_tfs(st, hs).freq = state_tfs.freq; 
                        session_ecg_tfs(i).condition(cn).hs_tuned_tfs(st, hs).cfg = state_tfs.cfg;
                        session_ecg_tfs(i).condition(cn).hs_tuned_tfs(st, hs).hs_label = hs_labels(hs);
                        if isfield(state_tfs, 'state_id') && isfield(state_tfs, 'state_name')
                            session_ecg_tfs(i).condition(cn).hs_tuned_tfs(st, hs).state = state_tfs.state_id;
                            session_ecg_tfs(i).condition(cn).hs_tuned_tfs(st, hs).state_name = state_tfs.state_name;
                        end
                        session_ecg_tfs(i).condition(cn).hs_tuned_tfs(st, hs).trials = find(cond_trials);
                        session_ecg_tfs(i).condition(cn).hs_tuned_tfs(st, hs).ntrials = length(find(cond_trials));

                    end

                end

            end
            
            % plots
            % Evoked LFP
            if ~isempty(fieldnames(session_ecg_tfs(i).condition(cn).hs_tuned_tfs))
                plottitle = ['Session: ', session_ecg_tfs(i).session ...
                    conditions(cn).label '), '];
                if conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                elseif conditions(cn).choice == 1
                    plottitle = [plottitle 'Choice trials'];
                end
                result_file = fullfile(results_folder_tfs, ...
                    ['TFS_ECG_' session_ecg_tfs(i).session '_' conditions(cn).label '.png']);

                lfp_tfa_plot_hs_tuned_tfr_multiple_img (session_ecg_tfs(i).condition(cn).hs_tuned_tfs, lfp_tfa_cfg, ...
                    plottitle, result_file);
            end

        end
        
        % difference between conditions
        session_ecg_tfs(i).difference = [];
        for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
            diff_condition = lfp_tfa_cfg.diff_condition{diff};
            session_ecg_tfs(i).difference = [session_ecg_tfs(i).difference, ...
                lfp_tfa_compute_difference_condition_tfr(session_ecg_tfs(i).condition, diff_condition)];
        end
        % plot Difference TFR
        for dcn = 1:length(session_ecg_tfs(i).difference)
            if ~isempty(session_ecg_tfs(i).difference(dcn).hs_tuned_tfs)
                if isfield(session_ecg_tfs(i).difference(dcn).hs_tuned_tfs,... 
                        'powspctrm')
                    plottitle = ['Session: ', session_ecg_tfs(i).session ...
                        session_ecg_tfs(i).difference(dcn).label];
                    result_file = fullfile(results_folder_tfs, ...
                        ['ECG_DiffTFS_' session_ecg_tfs(i).session ...
                        '_' 'diff_condition' num2str(dcn) '.png']);
                        %sites_avg(t).difference(dcn).label '.png']);
                    lfp_tfa_plot_hs_tuned_tfr_multiple_img(session_ecg_tfs(i).difference(dcn).hs_tuned_tfs, ...
                        lfp_tfa_cfg, plottitle, result_file, 'bluewhitered');
                end
            end
        end
        
        site_evoked_ecg = session_ecg_tfs(i);
        % save mat file for site
        save(fullfile(results_folder_tfs, ['TFS_ECG_' session_ecg_tfs(i).session '.mat']), 'site_evoked_ecg');
        % save to a mother struct
        session_ecg_tfs.sites(i) = site_evoked_ecg;
    end
        
%     % Average across sites for a session
%     session_avg = struct();
%     % targets for this session
%     targets = unique({session_ecg.target});
%     % average each target separately
%     for t = 1:length(targets)
%         session_avg(t).target = targets{t};
%         
%         % conditions
%         for cn = 1:length(site_conditions)
%             session_avg(t).condition(cn).hs_tuned_evoked = struct();
%             isite = 0;
%             session_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
%             session_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
%             for st = 1:size(sites_evoked(1).condition(1).hs_tuned_evoked, 1)
%                 for hs = 1:size(sites_evoked(1).condition(1).hs_tuned_evoked, 2)
%                     session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = 0;
%                     session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg = [];
%                     session_avg(t).condition(cn).hs_tuned_evoked(st, hs).mean = [];
%                     session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = [];
%                     session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time = [];
%                 end
%             end
%             for i = 1:length(sites_evoked)
%                 % if the site's target is same as target being considered
%                 if ~strcmp(session_ecg(i).target, targets{t})
%                     continue;
%                 end
%                 if sites_evoked(i).use_for_avg
%                     % calculate the average evoked ECG across sites for this condition 
%                     if ~isempty(sites_evoked(i).condition(cn).hs_tuned_evoked) && ...
%                         isfield(sites_evoked(i).condition(cn).hs_tuned_evoked, 'mean')
%                         isite = isite + 1;
%                         for hs = 1:size(sites_evoked(i).condition(cn).hs_tuned_evoked, 2)
%                             for st = 1:size(sites_evoked(i).condition(cn).hs_tuned_evoked, 1)
%                                 if ~isempty(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean)
%                                     session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = ...
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites + 1;
%                                     if session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites == 1
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg = ...
%                                             [session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg; ...
%                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean] ;
% %                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = ...
% %                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).std ;
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time = ...
%                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time;
% 
%                                     else
%                                         nsamples = length(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time);
%                                         if nsamples > length(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time)
%                                             nsamples = length(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).time);
%                                         end
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg = ...
%                                             [session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg(:,1:nsamples); ...
%                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).mean(1:nsamples)] ;
% %                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = ...
% %                                             session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std(1:nsamples) + ...
% %                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).std(1:nsamples) ;
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time = ...
%                                             session_avg(t).condition(cn).hs_tuned_evoked(st, hs).time(1:nsamples) ;
% 
%                                     end
%                                     % struct to store average evoked LFP across sites
%                                     session_avg(t).condition(cn).hs_tuned_evoked(st, hs).hs_label = ...
%                                         sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).hs_label;
%                                     if isfield(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs), 'state') && ...
%                                             isfield(sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs), 'state_name')
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).state = ...
%                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state;
%                                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).state_name = ...
%                                             sites_evoked(i).condition(cn).hs_tuned_evoked(st, hs).state_name;
%                                     end
%                                     session_avg(t).condition(cn).condition = site_conditions(cn);
%                                     session_avg(t).condition(cn).label = site_conditions(cn).label;
%                                     session_avg(t).condition(cn).session = session_ecg(i).session;
%                                     session_avg(t).condition(cn).target = session_ecg(i).target;
%                                     %session_avg(t).condition(cn).nsites = nsites;
%                                 end
% 
%                             end
%                         end
%                     end
%                 end
%             end
%             % average TFR across sites for a session
%             if isfield(session_avg(t).condition(cn).hs_tuned_evoked, 'ecg') && ...
%                     ~isempty([session_avg(t).condition(cn).hs_tuned_evoked.ecg])
%                 
%                 for hs = 1:size(session_avg(t).condition(cn).hs_tuned_evoked, 2)
%                     for st = 1:size(session_avg(t).condition(cn).hs_tuned_evoked, 1)
%                         if ~isempty(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg)
%                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).mean = ...
%                             nanmean(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg, 1);
%                         session_avg(t).condition(cn).hs_tuned_evoked(st, hs).std = ...
%                             nanstd(session_avg(t).condition(cn).hs_tuned_evoked(st, hs).ecg, 0, 1);
%                         %session_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = isite;
%                         end
%                     end
%                 end
%             end           
%             % plot average evoked ECG across sites for this session
%             if ~isempty([session_avg(t).condition(cn).hs_tuned_evoked.ecg])
%                 plottitle = ['Session: ', session_avg(t).condition(cn).session ', Target = ' ...
%                     session_avg(t).condition(cn).target ' (ref_' lfp_tfa_cfg.ref_hemisphere '), '  ...
%                     site_conditions(cn).label ', '];
%                 if site_conditions(cn).choice == 0
%                     plottitle = [plottitle 'Instructed trials'];
%                 elseif site_conditions(cn).choice == 1
%                     plottitle = [plottitle 'Choice trials'];
%                 end
%                 result_file = fullfile(results_folder_evoked, ['ECG_Evoked_' ...
%                     session_avg(t).condition(cn).session '_' ...
%                     session_avg(t).condition(cn).target '_' site_conditions(cn).label '.png']);
%                 lfp_tfa_plot_evoked_lfp (session_avg(t).condition(cn).hs_tuned_evoked, lfp_tfa_cfg, ...
%                     plottitle, result_file);
%             end
%         end 
%         
%         
%         % difference between conditions
%         session_avg(t).difference = [];
%         for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
%             diff_condition = lfp_tfa_cfg.diff_condition{diff};
%             session_avg(t).difference = [session_avg(t).difference, ...
%                 lfp_tfa_compute_diff_condition_evoked(session_avg(t).condition, diff_condition)];
%         end
%         % plot Difference TFR
%         for dcn = 1:length(session_avg(t).difference)
%             if ~isempty(session_avg(t).difference(dcn).hs_tuned_evoked)
%                 if isfield(session_avg(t).difference(dcn).hs_tuned_evoked,... 
%                         'mean')
%                     plottitle = ['Target ', lfp_tfa_cfg.compare.targets{t}, ...
%                     ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') ', ...
%                     session_avg(t).difference(dcn).label];
%                     result_file = fullfile(results_fldr, ...
%                         ['ECG_DiffEvoked_' lfp_tfa_cfg.compare.targets{t} ...
%                         '_' 'diff_condition' num2str(dcn) '.png']);
%                         %session_avg(t).difference(dcn).label '.png']);
%                     lfp_tfa_plot_evoked_lfp(session_avg(t).difference(dcn).hs_tuned_evoked, ...
%                         lfp_tfa_cfg, plottitle, result_file);
%                 end
%             end
%         end
% 
%     end
%         
%     close all;
%     
%     % store session average data
%     session_evoked_ecg.session_avg = session_avg;
%     
%     % save mat files
%     save(fullfile(results_folder_evoked, ['ECG_evoked_' session_evoked_ecg.session '.mat']), 'session_evoked');
%     % save settings file
%     save(fullfile(results_folder_evoked, 'lfp_tfa_settings.mat'), 'lfp_tfa_cfg');
end
        
