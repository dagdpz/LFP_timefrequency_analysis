function sites_avg = lfp_tfa_avg_pow_across_sites(lfp_pow, lfp_tfa_cfg)

    % results folder
    results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sites', 'LFP_Power');
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % Average TFR across sessions
    sites_avg.condition = struct();
    for cn = 1:length(lfp_pow.session(1).condition)
        fprintf('Condition %s\n', lfp_pow.session(1).condition(cn).label);
        sites_avg.condition(cn).label = lfp_pow.session(1).condition(cn).label;
        sites_avg.condition(cn).avg_across_sessions = struct();
        nsites = 0;
        %sessions_avg.cond_based_tfs(cn).tfs_across_sessions = struct();
        for i = 1:length(lfp_pow.session)  
            for j = 1:length(lfp_pow.session(i).sites)
                if ~isempty(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power) && ... 
                    isfield(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power, 'mean')
                    nsites = nsites + 1;   
                    %sessions_avg.cond_based_tfs(cn).tfs_across_sessions = struct();
                    for ep = 1:size(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power, 1)
                        for hs = 1:size(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power, 2)
                            if isfield(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs), 'mean') ...
                                    && ~isempty(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).mean)
                                if nsites == 1%~isfield(sessions_avg.cond_based_tfs(cn).tfs_across_sessions, 'powspctrm')
                                    sites_avg.condition(cn).avg_across_sessions(ep,hs).freq ...
                                    = lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).freq;
                                    sites_avg.condition(cn).avg_across_sessions(ep,hs).hs_label ...
                                        = lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).hs_label;
                                    sites_avg.condition(cn).avg_across_sessions(ep,hs).mean ...
                                        = lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).mean;
                                    if isfield(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs), 'nsites')
                                        sites_avg.condition(cn).avg_across_sessions(ep,hs).nsites ...
                                            = lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).nsites;
                                    end
                                else
                                    nfreqbins = length(sites_avg.condition(cn).avg_across_sessions(ep, hs).freq);
                                    % average same number of time bins
                                    if nfreqbins > length(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).freq)
                                        nfreqbins = length(lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).freq);
                                    end
                                    sites_avg.condition(cn).avg_across_sessions(ep,hs).mean ...
                                        = (lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).mean(1:nfreqbins)) + ...
                                        sites_avg.condition(cn).avg_across_sessions(ep,hs).mean(1:nfreqbins);
                                    if isfield(sites_avg.condition(cn).avg_across_sessions(ep,hs), 'nsites')
                                        sites_avg.condition(cn).avg_across_sessions(ep,hs).nsites ...
                                            = lfp_pow.session(i).sites(j).condition(cn).hs_tuned_power(ep, hs).nsites + ...
                                            sites_avg.condition(cn).avg_across_sessions(ep,hs).nsites;
                                    end                               
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % compute average
        for ep = 1:size(sites_avg.condition(cn).avg_across_sessions, 1)
            for hs = 1:size(sites_avg.condition(cn).avg_across_sessions, 2)
                if isfield(sites_avg.condition(cn).avg_across_sessions(ep,hs), 'mean')
                    sites_avg.condition(cn).avg_across_sessions(ep,hs).nsites = nsites;                                
                    sites_avg.condition(cn).avg_across_sessions(ep,hs).mean = ...
                        (1/nsites) * sites_avg.condition(cn).avg_across_sessions(ep,hs).mean;
                end
            end
        end
        
        
        if ~isempty(sites_avg.condition(cn).avg_across_sessions)
            if isfield(sites_avg.condition(cn).avg_across_sessions,... 
                    'mean')
                plottitle = lfp_pow.session(1).condition(cn).label;
                result_file = fullfile(results_fldr, ...
                                ['LFP_Power_' lfp_pow.session(1).condition(cn).label '.png']);
                lfp_tfa_plot_hs_tuned_psd(sites_avg.condition(cn).avg_across_sessions, ...
                            lfp_tfa_cfg, plottitle, result_file);
            end
        end
        % save session average tfs
        save(fullfile(results_fldr, 'LFP_Power_sites_avg.mat'), 'sites_avg');
    end
end