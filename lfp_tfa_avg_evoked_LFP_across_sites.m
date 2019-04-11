function sites_avg = lfp_tfa_avg_evoked_LFP_across_sites(lfp_evoked, lfp_tfa_cfg)

    % results folder
    results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sites', 'LFP_Evoked');
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % Average TFR across sessions
    sites_avg.condition = struct();
    for cn = 1:length(lfp_evoked.session(1).condition)
        fprintf('Condition %s\n', lfp_evoked.session(1).condition(cn).label);
        sites_avg.condition(cn).hs_tuned_evoked = struct();
        sites_avg.condition(cn).label = lfp_evoked.session(1).condition(cn).label;
        nsites = 0;
        %sessions_avg.cond_based_tfs(cn).tfs_across_sessions = struct();
        for i = 1:length(lfp_evoked.session) 
            for j = 1:length(lfp_evoked.session(i).sites)
                if ~isempty(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked) && ... 
                    isfield(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked, 'mean')
                    nsites = nsites + 1;   
                    %sessions_avg.cond_based_tfs(cn).tfs_across_sessions = struct();
                    for st = 1:size(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked, 1)
                        for hs = 1:size(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked, 2)
                            if isfield(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs), 'mean') ...
                                    && ~isempty(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).mean)
                                if nsites == 1%~isfield(sessions_avg.condition(cn).hs_tuned_evoked, 'mean')
                                    sites_avg.condition(cn).hs_tuned_evoked(st,hs).time ...
                                    = lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).time;
                                    sites_avg.condition(cn).hs_tuned_evoked(st,hs).hs_label ...
                                        = lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).hs_label;
                                    sites_avg.condition(cn).hs_tuned_evoked(st,hs).state ...
                                        = lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).state;
                                    sites_avg.condition(cn).hs_tuned_evoked(st,hs).mean ...
                                        = lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).mean;
                                    sites_avg.condition(cn).hs_tuned_evoked(st,hs).std ...
                                        = lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).std;                                
                                    if isfield(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs), 'nsites')
                                        sites_avg.condition(cn).hs_tuned_evoked(st,hs).nsites ...
                                            = lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).nsites;
                                    end
                                else
                                    ntimebins = length(sites_avg.condition(cn).hs_tuned_evoked(st, hs).time);
                                    % average same number of time bins
                                    if ntimebins > length(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).time)
                                        ntimebins = length(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).time);
                                    end
                                    sites_avg.condition(cn).hs_tuned_evoked(st,hs).mean ...
                                        = (lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).mean(1:ntimebins)) + ...
                                        sites_avg.condition(cn).hs_tuned_evoked(st,hs).mean(1:ntimebins);
                                    sites_avg.condition(cn).hs_tuned_evoked(st,hs).std ...
                                        = (lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).std(1:ntimebins)) + ...
                                        sites_avg.condition(cn).hs_tuned_evoked(st,hs).std(1:ntimebins);
                                    if isfield(sites_avg.condition(cn).hs_tuned_evoked(st,hs), 'nsites')
                                        sites_avg.condition(cn).hs_tuned_evoked(st,hs).nsites ...
                                            = lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).nsites + ...
                                            sites_avg.condition(cn).hs_tuned_evoked(st,hs).nsites;
                                    end                               
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % compute average
        for st = 1:size(sites_avg.condition(cn).hs_tuned_evoked, 1)
            for hs = 1:size(sites_avg.condition(cn).hs_tuned_evoked, 2)
                if isfield(sites_avg.condition(cn).hs_tuned_evoked(st,hs), 'mean')
                    sites_avg.condition(cn).hs_tuned_evoked(st,hs).nsites = nsites;                                
                    sites_avg.condition(cn).hs_tuned_evoked(st,hs).mean = ...
                        (1/nsites) * sites_avg.condition(cn).hs_tuned_evoked(st,hs).mean;
                    sites_avg.condition(cn).hs_tuned_evoked(st,hs).std = ...
                        (1/nsites) * sites_avg.condition(cn).hs_tuned_evoked(st,hs).std;
                end
            end
        end
        
        
        if ~isempty(sites_avg.condition(cn).hs_tuned_evoked)
            if isfield(sites_avg.condition(cn).hs_tuned_evoked,... 
                    'mean')
                plottitle = lfp_evoked.session(1).condition(cn).label;
                result_file = fullfile(results_fldr, ...
                                ['LFP_Evoked_' lfp_evoked.session(1).condition(cn).label '.png']);
                lfp_tfa_plot_evoked_lfp(sites_avg.condition(cn).hs_tuned_evoked, ...
                            lfp_tfa_cfg, plottitle, result_file);
            end
        end
        % save session average tfs
        save(fullfile(results_fldr, 'LFP_Evoked_sites_average.mat'), 'sites_avg');
    end
    close all;
end