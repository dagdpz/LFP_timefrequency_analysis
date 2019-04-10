function sessions_avg = lfp_tfa_avg_pow_per_session(lfp_pow, lfp_tfa_cfg)

    % results folder
    results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Average_Power');
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % Average TFR across sessions
    sessions_avg.cond_based_pow = struct();
    for cn = 1:length(lfp_pow.session(1).condition)
        fprintf('Condition %s\n', lfp_pow.session(1).condition(cn).label);
        sessions_avg.condition(cn).avg_across_sessions = struct();
        nsessions = 0;
        %sessions_avg.cond_based_tfs(cn).tfs_across_sessions = struct();
        for i = 1:length(lfp_pow.session)  
            if ~isempty(lfp_pow.session(i).condition(cn).session_lfp_psd) && ... 
                isfield(lfp_pow.session(i).condition(cn).session_lfp_psd, 'mean')
                nsessions = nsessions + 1;   
                %sessions_avg.cond_based_tfs(cn).tfs_across_sessions = struct();
                for st = 1:size(lfp_pow.session(i).condition(cn).session_lfp_psd, 1)
                    for hs = 1:size(lfp_pow.session(i).condition(cn).session_lfp_psd, 2)
                        if isfield(lfp_pow.session(i).condition(cn).session_lfp_psd(st, hs), 'mean') ...
                                && ~isempty(lfp_pow.session(i).condition(cn).session_lfp_psd(st, hs).mean)
                            if nsessions == 1%~isfield(sessions_avg.cond_based_tfs(cn).tfs_across_sessions, 'powspctrm')
                                sessions_avg.condition(cn).avg_across_sessions(st,hs).freq ...
                                = lfp_pow.session(i).condition(cn).session_lfp_psd(st, hs).freq;
                                sessions_avg.condition(cn).avg_across_sessions(st,hs).hs_label ...
                                    = lfp_pow.session(i).condition(cn).session_lfp_psd(st, hs).hs_label;
                                sessions_avg.condition(cn).avg_across_sessions(st,hs).mean ...
                                    = lfp_pow.session(i).condition(cn).session_lfp_psd(st, hs).mean;
                                if isfield(lfp_pow.session(i).condition(cn).session_lfp_psd(st, hs), 'nsites')
                                    sessions_avg.condition(cn).avg_across_sessions(st,hs).nsites ...
                                        = lfp_pow.session(i).condition(cn).session_lfp_psd(st, hs).nsites;
                                end
                            else
                                nfreqbins = length(sessions_avg.condition(cn).avg_across_sessions(st, hs).freq);
                                % average same number of time bins
                                if nfreqbins > length(lfp_pow.session(i).condition(cn).session_lfp_psd(st, hs).freq)
                                    nfreqbins = length(lfp_pow.session(i).condition(cn).session_lfp_psd(st, hs).freq);
                                end
                                sessions_avg.condition(cn).avg_across_sessions(st,hs).mean ...
                                    = (lfp_pow.session(i).condition(cn).session_lfp_psd(st, hs).mean(1:nfreqbins)) + ...
                                    sessions_avg.condition(cn).avg_across_sessions(st,hs).mean(1:nfreqbins);
                                if isfield(sessions_avg.condition(cn).avg_across_sessions(st,hs), 'nsites')
                                    sessions_avg.condition(cn).avg_across_sessions(st,hs).nsites ...
                                        = lfp_pow.session(i).condition(cn).session_lfp_psd(st, hs).nsites + ...
                                        sessions_avg.condition(cn).avg_across_sessions(st,hs).nsites;
                                end                               
                            end
                        end
                    end
                end
            end
        end
        
        % compute average
        for st = 1:size(sessions_avg.condition(cn).avg_across_sessions, 1)
            for hs = 1:size(sessions_avg.condition(cn).avg_across_sessions, 2)
                if isfield(sessions_avg.condition(cn).avg_across_sessions(st,hs), 'mean')
                    sessions_avg.condition(cn).avg_across_sessions(st,hs).nsessions = nsessions;                                
                    sessions_avg.condition(cn).avg_across_sessions(st,hs).mean = ...
                        (1/nsessions) * sessions_avg.condition(cn).avg_across_sessions(st,hs).mean;
                end
            end
        end
        
        
        if ~isempty(sessions_avg.condition(cn).avg_across_sessions)
            if isfield(sessions_avg.condition(cn).avg_across_sessions,... 
                    'mean')
                plottitle = lfp_pow.session(1).condition(cn).label;
                result_file = fullfile(results_fldr, ...
                                ['Avg_LFP_Power_' lfp_pow.session(1).condition(cn).label '.png']);
                lfp_tfa_plot_hs_tuned_psd(sessions_avg.condition(cn).avg_across_sessions, ...
                            lfp_tfa_cfg, plottitle, result_file);
            end
        end
        % save session average tfs
        save(fullfile(results_fldr, 'sessions_average_evoked.mat'), 'sessions_avg');
    end
end