function sessions_avg = lfp_tfa_avg_evoked_LFP_across_sessions(lfp_evoked, lfp_tfa_cfg)

    % results folder
    results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sessions', 'LFP_Evoked');
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % Average TFR across sessions
    sessions_avg.condition = struct();
    for cn = 1:length(lfp_evoked.session(1).condition)
        fprintf('Condition %s\n', lfp_evoked.session(1).condition(cn).label);
        sessions_avg.condition(cn).avg_across_sessions = struct();
        nsessions = 0;
        %sessions_avg.cond_based_tfs(cn).tfs_across_sessions = struct();
        for i = 1:length(lfp_evoked.session)  
            if ~isempty(lfp_evoked.session(i).condition(cn).hs_tuned_evoked) && ... 
                isfield(lfp_evoked.session(i).condition(cn).hs_tuned_evoked, 'mean')
                nsessions = nsessions + 1;   
                %sessions_avg.cond_based_tfs(cn).tfs_across_sessions = struct();
                for st = 1:size(lfp_evoked.session(i).condition(cn).hs_tuned_evoked, 1)
                    for hs = 1:size(lfp_evoked.session(i).condition(cn).hs_tuned_evoked, 2)
                        if isfield(lfp_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs), 'mean') ...
                                && ~isempty(lfp_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).mean)
                            if nsessions == 1%~isfield(sessions_avg.condition(cn).avg_across_sessions, 'mean')
                                sessions_avg.condition(cn).avg_across_sessions(st,hs).time ...
                                = lfp_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).time;
                                sessions_avg.condition(cn).avg_across_sessions(st,hs).hs_label ...
                                    = lfp_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).hs_label;
                                sessions_avg.condition(cn).avg_across_sessions(st,hs).state ...
                                    = lfp_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).state;
                                sessions_avg.condition(cn).avg_across_sessions(st,hs).mean ...
                                    = lfp_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).mean;
                                sessions_avg.condition(cn).avg_across_sessions(st,hs).std ...
                                    = lfp_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).std;                                
                                if isfield(lfp_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs), 'nsites')
                                    sessions_avg.condition(cn).avg_across_sessions(st,hs).nsites ...
                                        = lfp_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).nsites;
                                end
                            else
                                ntimebins = length(sessions_avg.condition(cn).avg_across_sessions(st, hs).time);
                                % average same number of time bins
                                if ntimebins > length(lfp_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).time)
                                    ntimebins = length(lfp_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).time);
                                end
                                sessions_avg.condition(cn).avg_across_sessions(st,hs).mean ...
                                    = (lfp_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).mean(1:ntimebins)) + ...
                                    sessions_avg.condition(cn).avg_across_sessions(st,hs).mean(1:ntimebins);
                                sessions_avg.condition(cn).avg_across_sessions(st,hs).std ...
                                    = (lfp_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).std(1:ntimebins)) + ...
                                    sessions_avg.condition(cn).avg_across_sessions(st,hs).std(1:ntimebins);
                                if isfield(sessions_avg.condition(cn).avg_across_sessions(st,hs), 'nsites')
                                    sessions_avg.condition(cn).avg_across_sessions(st,hs).nsites ...
                                        = lfp_evoked.session(i).condition(cn).hs_tuned_evoked(st, hs).nsites + ...
                                        sessions_avg.condition(cn).avg_across_sessions(st,hs).nsites;
                                end                               
                            end
                        end
                    end
                end
            end
        end
        
        % compute average
        if isfield(sessions_avg.condition(cn).avg_across_sessions, 'mean')
            for st = 1:size(sessions_avg.condition(cn).avg_across_sessions, 1)
                for hs = 1:size(sessions_avg.condition(cn).avg_across_sessions, 2)
                    if ~isempty(sessions_avg.condition(cn).avg_across_sessions(st,hs).mean)
                        sessions_avg.condition(cn).avg_across_sessions(st,hs).nsessions = nsessions;                                
                        sessions_avg.condition(cn).avg_across_sessions(st,hs).mean = ...
                            (1/nsessions) * sessions_avg.condition(cn).avg_across_sessions(st,hs).mean;
                        sessions_avg.condition(cn).avg_across_sessions(st,hs).std = ...
                            (1/nsessions) * sessions_avg.condition(cn).avg_across_sessions(st,hs).std;
                    end
                end
            end
        end
        
        
        if ~isempty(sessions_avg.condition(cn).avg_across_sessions)
            if isfield(sessions_avg.condition(cn).avg_across_sessions,... 
                    'mean')
                plottitle = lfp_evoked.session(1).condition(cn).label;
                result_file = fullfile(results_fldr, ...
                                ['LFP_Evoked_' lfp_evoked.session(1).condition(cn).label '.png']);
                lfp_tfa_plot_evoked_lfp(sessions_avg.condition(cn).avg_across_sessions, ...
                            lfp_tfa_cfg, plottitle, result_file);
            end
        end
        % save session average tfs
        save(fullfile(results_fldr, 'sessions_average_evoked.mat'), 'sessions_avg');
    end
    close all;
end