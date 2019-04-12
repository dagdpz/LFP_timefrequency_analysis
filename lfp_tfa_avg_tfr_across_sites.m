function sites_avg = lfp_tfa_avg_tfr_across_sites(lfp_tfr, lfp_tfa_cfg)

    % results folder
    results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sites', 'LFP_TFR');
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % Average TFR across sessions
    sites_avg.condition = struct();
    for cn = 1:length(lfp_tfr(1).session(1).condition)
        fprintf('Condition %s\n', lfp_tfr.session(1).condition(cn).label);
        sites_avg.condition(cn).hs_tuned_tfs = struct();
        nsites = 0;
        %sessions_avg.condition(cn).hs_tuned_tfs = struct();
        for i = 1:length(lfp_tfr.session)  
            for j = 1:length(lfp_tfr.session(i).sites)
                if ~isempty(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs) && ... 
                    isfield(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs, 'powspctrm')
                    nsites = nsites + 1;   
                    %sessions_avg.condition(cn).tfs_across_sessions = struct();
                    for st = 1:size(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs, 1)
                        for hs = 1:size(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs, 2)
                            if isfield(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs), 'powspctrm') ...
                                    && ~isempty(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs).powspctrm)
                                if nsites == 1%~isfield(sessions_avg.cond_based_tfs(cn).tfs_across_sessions, 'powspctrm')
                                    sites_avg.condition(cn).hs_tuned_tfs(st,hs).time ...
                                    = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs).time;
                                    sites_avg.condition(cn).hs_tuned_tfs(st,hs).hs_label ...
                                        = lfp_tfr.session(i).condition(cn).hs_tuned_tfs(st, hs).hs_label;
                                    sites_avg.condition(cn).hs_tuned_tfs(st,hs).state ...
                                        = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs).state;
                                    sites_avg.condition(cn).hs_tuned_tfs(st,hs).cfg ...
                                        = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs).cfg;
                                    sites_avg.condition(cn).hs_tuned_tfs(st,hs).freq ...
                                        = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs).freq;
                                    sites_avg.condition(cn).hs_tuned_tfs(st,hs).powspctrm ...
                                        = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs).powspctrm;
                                    sites_avg.condition(cn).label = ...
                                        lfp_tfr.session(i).condition(cn).label;
                                    sites_avg.condition(cn).cfg_condition = ...
                                        lfp_tfr.session(i).condition(cn).cfg_condition;
                                    if isfield(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs), 'nsites')
                                        sites_avg.condition(cn).hs_tuned_tfs(st,hs).nsites ...
                                            = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs).nsites;
                                    end
                                else
                                    ntimebins = size(sites_avg.condition(cn).hs_tuned_tfs(st, hs).powspctrm, 3);
                                    % average same number of time bins
                                    if ntimebins > length(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs).time)
                                        ntimebins = length(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs).time);
                                    end
                                    sites_avg.condition(cn).hs_tuned_tfs(st,hs).powspctrm ...
                                        = (lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins)) + ...
                                        sites_avg.condition(cn).hs_tuned_tfs(st,hs).powspctrm(1,:,1:ntimebins);
                                    if isfield(sites_avg.condition(cn).hs_tuned_tfs(st,hs), 'nsites')
                                        sites_avg.condition(cn).hs_tuned_tfs(st,hs).nsites ...
                                            = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_tfs(st, hs).nsites + ...
                                            sites_avg.condition(cn).hs_tuned_tfs(st,hs).nsites;
                                    end                               
                                end
                            end
                        end
                    end
                end
            end
        end
        
        % compute average
        if ~isempty(sites_avg.condition(cn).hs_tuned_tfs)
            if isfield(sites_avg.condition(cn).hs_tuned_tfs, 'powspctrm')
                for st = 1:size(sites_avg.condition(cn).hs_tuned_tfs, 1)
                    for hs = 1:size(sites_avg.condition(cn).hs_tuned_tfs, 2)
                        sites_avg.condition(cn).hs_tuned_tfs(st,hs).nsites = nsites;                                
                        sites_avg.condition(cn).hs_tuned_tfs(st,hs).powspctrm = ...
                            (1/nsites) * sites_avg.condition(cn).hs_tuned_tfs(st,hs).powspctrm;
                    end
                end
            end
        end
        
        
        if ~isempty(sites_avg.condition(cn).hs_tuned_tfs)
            if isfield(sites_avg.condition(cn).hs_tuned_tfs,... 
                    'powspctrm')
                plottitle = sites_avg.condition(cn).label;
                result_file = fullfile(results_fldr, ...
                                ['LFP_TFR_' sites_avg.condition(cn).label '.png']);
                lfp_tfa_plot_hs_tuned_tfr(sites_avg.condition(cn).hs_tuned_tfs, ...
                            lfp_tfa_cfg, plottitle, result_file);
            end
        end        
    end
    % difference between pre- and post-injection
    sites_avg.difference = lfp_tfa_compute_diff_tfr(sites_avg, lfp_tfa_cfg);
    
    % plot Difference TFR
    for dcn = 1:length(sites_avg.difference)
        if ~isempty(sites_avg.difference(dcn).hs_tuned_tfs)
            if isfield(sites_avg.difference(dcn).hs_tuned_tfs,... 
                    'powspctrm')
                plottitle = sites_avg.difference(dcn).label;
                result_file = fullfile(results_fldr, ...
                                ['LFP_DiffTFR_' sites_avg.difference(dcn).label '.png']);
                lfp_tfa_plot_hs_tuned_tfr(sites_avg.difference(dcn).hs_tuned_tfs, ...
                            lfp_tfa_cfg, plottitle, result_file);
            end
        end
    end
    
    % save session average tfs
    save(fullfile(results_fldr, 'LFP_TFR_sites_avg.mat'), 'sites_avg');
    
end