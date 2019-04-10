function sessions_avg = lfp_tfa_avg_tfr_per_session(lfp_tfr, lfp_tfa_cfg)

    % results folder
    results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Average_TFR');
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % Average TFR across sessions
    sessions_avg.condition = struct();
    for cn = 1:length(lfp_tfr(1).session(1).condition)
        fprintf('Condition %s\n', lfp_tfr(1).session(1).condition(cn).label);
        sessions_avg.condition(cn).tfs_across_sessions = struct();
        nsessions = 0;
        %sessions_avg.condition(cn).tfs_across_sessions = struct();
        for i = 1:length(lfp_tfr.session)  
            if ~isempty(lfp_tfr.session(i).condition(cn).hs_tuned_tfs) && ... 
                isfield(lfp_tfr.session(i).condition(cn).hs_tuned_tfs, 'powspctrm')
                nsessions = nsessions + 1;   
                %sessions_avg.condition(cn).tfs_across_sessions = struct();
                for st = 1:size(lfp_tfr.session(i).condition(cn).hs_tuned_tfs, 1)
                    for hs = 1:size(lfp_tfr.session(i).condition(cn).hs_tuned_tfs, 2)
                        if isfield(lfp_tfr.session(i).condition(cn).hs_tuned_tfs(st, hs), 'powspctrm') ...
                                && ~isempty(lfp_tfr.session(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm)
                            if nsessions == 1%~isfield(sessions_avg.cond_based_tfs(cn).tfs_across_sessions, 'powspctrm')
                                sessions_avg.condition(cn).tfs_across_sessions(st,hs).time ...
                                = lfp_tfr.session(i).condition(cn).hs_tuned_tfs(st, hs).time;
                                sessions_avg.condition(cn).tfs_across_sessions(st,hs).hs_label ...
                                    = lfp_tfr.session(i).condition(cn).hs_tuned_tfs(st, hs).hs_label;
                                sessions_avg.condition(cn).tfs_across_sessions(st,hs).state ...
                                    = lfp_tfr.session(i).condition(cn).hs_tuned_tfs(st, hs).state;
                                sessions_avg.condition(cn).tfs_across_sessions(st,hs).cfg ...
                                    = lfp_tfr.session(i).condition(cn).hs_tuned_tfs(st, hs).cfg;
                                sessions_avg.condition(cn).tfs_across_sessions(st,hs).freq ...
                                    = lfp_tfr.session(i).condition(cn).hs_tuned_tfs(st, hs).freq;
                                sessions_avg.condition(cn).tfs_across_sessions(st,hs).powspctrm ...
                                    = lfp_tfr.session(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm;
                                if isfield(sessions_avg.condition(cn).tfs_across_sessions(st,hs), 'nsites')
                                    sessions_avg.condition(cn).tfs_across_sessions(st,hs).nsites ...
                                        = sessions_avg.condition(cn).tfs_across_sessions(st,hs).nsites;
                                end
                            else
                                ntimebins = size(sessions_avg.condition(cn).tfs_across_sessions(st, hs).powspctrm, 3);
                                % average same number of time bins
                                if ntimebins > length(lfp_tfr.session(i).condition(cn).hs_tuned_tfs(st, hs).time)
                                    ntimebins = length(lfp_tfr.session(i).condition(cn).hs_tuned_tfs(st, hs).time);
                                end
                                sessions_avg.condition(cn).tfs_across_sessions(st,hs).powspctrm ...
                                    = (lfp_tfr.session(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins)) + ...
                                    sessions_avg.condition(cn).tfs_across_sessions(st,hs).powspctrm(1,:,1:ntimebins);
                                if isfield(sessions_avg.condition(cn).tfs_across_sessions(st,hs), 'nsites')
                                    sessions_avg.condition(cn).tfs_across_sessions(st,hs).nsites ...
                                        = lfp_tfr.session(i).condition(cn).hs_tuned_tfs(st, hs).nsites + ...
                                        sessions_avg.condition(cn).tfs_across_sessions(st,hs).nsites;
                                end                               
                            end
                        end
                    end
                end
            end
        end
        
        % compute average
        for st = 1:size(sessions_avg.condition(cn).tfs_across_sessions, 1)
            for hs = 1:size(sessions_avg.condition(cn).tfs_across_sessions, 2)
                sessions_avg.condition(cn).tfs_across_sessions(st,hs).nsessions = nsessions;                                
                sessions_avg.condition(cn).tfs_across_sessions(st,hs).powspctrm = ...
                    (1/nsessions) * sessions_avg.condition(cn).tfs_across_sessions(st,hs).powspctrm;
            end
        end
        
        
        if ~isempty(sessions_avg.condition(cn).tfs_across_sessions)
            if isfield(sessions_avg.condition(cn).tfs_across_sessions,... 
                    'powspctrm')
                plottitle = lfp_tfr.session(1).condition(cn).label;
                result_file = fullfile(results_fldr, ...
                                ['Avg_Spectrogram_' lfp_tfr.session(1).condition(cn).label '.png']);
                lfp_tfa_plot_hs_tuned_tfr(sessions_avg.condition(cn).tfs_across_sessions, ...
                            lfp_tfa_cfg, plottitle, result_file);
            end
        end
        % save session average tfs
        save(fullfile(results_fldr, 'sessions_average_tfs.mat'), 'sessions_avg');
    end
end