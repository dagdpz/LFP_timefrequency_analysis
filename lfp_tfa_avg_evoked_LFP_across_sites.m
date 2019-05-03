function sites_avg = lfp_tfa_avg_evoked_LFP_across_sites(lfp_evoked, lfp_tfa_cfg)
%lfp_tfa_avg_evoked_LFP_across_sites  - Condition-based evoked LFP response
% average across many site averages
%
% USAGE:
%	sites_avg = lfp_tfa_avg_evoked_LFP_across_sessions(lfp_evoked, lfp_tfa_cfg)
%
% INPUTS:
%		lfp_evoked		- struct containing the condition-based evoked LFP response for
%		indiviual sites, output of lfp_tfa_plot_site_evoked_LFP.m
%           Required Fields:
%               1. session.sites       - 1xN struct containing condition-based
%               average evoked LFP response for N sites
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields:
%               1. conditions          - trial conditions to compare, see
%               lfp_tfa_settings.m, lfp_tfa_define_settings.m and lfp_tfa_compare_conditions.m
%               2. root_results_fldr   - root folder where results are saved
%               3. compare.targets     - targets to compare, see lfp_tfa_settings.m
%               4. 
% OUTPUTS:
%		sites_avg    - structure containing condition-based evoked LFP
%		response averaged across multiple sites
%
% REQUIRES:	lfp_tfa_plot_evoked_lfp
%
% See also lfp_tfa_settings, lfp_tfa_define_settings, lfp_tfa_compare_conditions, 
% lfp_tfa_plot_site_evoked_LFP
%
% Author(s):	S.Nair, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2019-02-15:	Created function (Sarath Nair)
% 2019-03-05:	First Revision
% ...
% $Revision: 1.0 $  $Date: 2019-03-05 17:18:00 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%

    % results folder
    results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sites', 'LFP_Evoked');
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % Average Evoked LFP across sites
    sites_avg = struct();
    
    for t = 1:length(lfp_tfa_cfg.compare.targets)
        sites_avg(t).target = lfp_tfa_cfg.compare.targets{t};
        for cn = 1:length(lfp_tfa_cfg.conditions)
            fprintf('Condition %s\n', lfp_tfa_cfg.conditions(cn).label);
            sites_avg(t).condition(cn).hs_tuned_evoked = struct();
            sites_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
            nsites = 0;
            %sessions_avg.cond_based_tfs(cn).tfs_across_sessions = struct();
            for i = 1:length(lfp_evoked.session) 
                for j = 1:length(lfp_evoked.session(i).sites)
                    if ~strcmp(lfp_evoked.session(i).sites(j).target, lfp_tfa_cfg.compare.targets{t})
                        continue;
                    end
                    if ~lfp_evoked.session(i).sites(j).use_for_avg
                        continue;
                    end
                    if ~isempty(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked) && ... 
                        isfield(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked, 'mean')
                        nsites = nsites + 1;   
                        %sessions_avg.cond_based_tfs(cn).tfs_across_sessions = struct();
                        for st = 1:size(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked, 1)
                            for hs = 1:size(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked, 2)
                                if isfield(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs), 'mean') ...
                                        && ~isempty(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).mean)
                                    if nsites == 1%~isfield(sessions_avg.condition(cn).hs_tuned_evoked, 'mean')
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).time ...
                                        = lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).time;
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).hs_label ...
                                            = lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).hs_label;
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).state ...
                                            = lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).state;
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).mean ...
                                            = lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).mean;
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).std ...
                                            = lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).std;                                
                                        if isfield(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs), 'nsites')
                                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).nsites ...
                                                = lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).nsites;
                                        end
                                    else
                                        ntimebins = length(sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).time);
                                        % average same number of time bins
                                        if ntimebins > length(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).time)
                                            ntimebins = length(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).time);
                                        end
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).time = ...
                                            lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).time(1:ntimebins);
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).mean ...
                                            = (lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).mean(1:ntimebins)) + ...
                                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).mean(1:ntimebins);
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).std ...
                                            = (lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).std(1:ntimebins)) + ...
                                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).std(1:ntimebins);
                                        if isfield(sites_avg(t).condition(cn).hs_tuned_evoked(st,hs), 'nsites')
                                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).nsites ...
                                                = lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).nsites + ...
                                                sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).nsites;
                                        end                               
                                    end
                                end
                            end
                        end
                    end
                end
            end

            % compute average
            for st = 1:size(sites_avg(t).condition(cn).hs_tuned_evoked, 1)
                for hs = 1:size(sites_avg(t).condition(cn).hs_tuned_evoked, 2)
                    if isfield(sites_avg(t).condition(cn).hs_tuned_evoked(st,hs), 'mean')
                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).nsites = nsites;                                
                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).mean = ...
                            (1/nsites) * sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).mean;
                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).std = ...
                            (1/nsites) * sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).std;
                    end
                end
            end

            % plot
            if ~isempty(sites_avg(t).condition(cn).hs_tuned_evoked)
                if isfield(sites_avg(t).condition(cn).hs_tuned_evoked,... 
                        'mean')
                    plottitle = [lfp_tfa_cfg.compare.targets{t} '_' lfp_tfa_cfg.conditions(cn).label];
                    result_file = fullfile(results_fldr, ...
                                    ['LFP_Evoked_' lfp_tfa_cfg.compare.targets{t} '_' lfp_tfa_cfg.conditions(cn).label '.png']);
                    lfp_tfa_plot_evoked_lfp(sites_avg(t).condition(cn).hs_tuned_evoked, ...
                                lfp_tfa_cfg, plottitle, result_file);
                end
            end            
        end
    end
    % save session average tfs
    save(fullfile(results_fldr, 'LFP_Evoked_sites_average.mat'), 'sites_avg');
    close all;
end