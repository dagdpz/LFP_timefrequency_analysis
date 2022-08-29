function sites_avg = lfp_tfa_avg_evoked_LFP_across_sites(Sessions, lfp_tfa_cfg, varargin)
%lfp_tfa_avg_evoked_LFP_across_sites  - Condition-based evoked LFP response
% grand average across many site averages
%
% USAGE:
%	sites_avg = lfp_tfa_avg_evoked_LFP_across_sites(lfp_evoked, lfp_tfa_cfg)
%
% INPUTS:
%		lfp_evoked		- struct containing the condition-based evoked LFP 
%       response average for indiviual sites, i.e., output of 
%       lfp_tfa_plot_site_evoked_LFP.m
%           Required Fields:
%               1. session.sites       - 1xN struct containing condition-based
%               average evoked LFP response for N sites
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields:
%               1. conditions          - trial conditions to compare, see
%               settings/lfp_tfa_settings_example, 
%               lfp_tfa_define_settings.m and lfp_tfa_compare_conditions.m
%               2. root_results_fldr   - root folder where results are 
%               saved. Results will be saved under 
%               [lfp_tfa_cfg.root_results_fldr ...
%               '/Avg_across_sites/LFP_Evoked']
%               3. compare.targets     - targets to compare, see lfp_tfa_settings.m
%               4. ref_hemisphere      - reference hemisphere for contra-
%               and ipsi- labelling, see settings/lfp_tfa_settings_example.m
% OUTPUTS:
%		sites_avg    - structure containing condition-based evoked LFP
%		response averaged across multiple sites
%
% REQUIRES:	lfp_tfa_plot_evoked_lfp
%
% See also see settings/lfp_tfa_settings_example.m, lfp_tfa_define_settings, 
% lfp_tfa_compare_conditions, lfp_tfa_plot_site_evoked_LFP
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
    if nargin < 3
        results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sites', 'LFP_Evoked');
    else
        results_fldr = varargin{1};
    end
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
            sites_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
            sites_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
            % initialize number of site pairs for each handspace
            % label
            for st = 1:size(Sessions(1).sites(1).condition(cn).hs_tuned_evoked, 1)
                for hs = 1:size(Sessions(1).sites(1).condition(cn).hs_tuned_evoked, 2)
                    sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = 0;
                    sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).lfp = [];
                end
            end
            for i = 1:length(Sessions) 
                for j = 1:length(Sessions(i).sites)
                    %% LS 2021 
                    if ~ismember(lfp_tfa_cfg.compare.targets{t}, Sessions(i).sites(j).target) 
                        continue;
                    end
                    if ~Sessions(i).sites(j).use_for_avg
                        continue;
                    end
                    if ~isempty(Sessions(i).sites(j).condition(cn).hs_tuned_evoked) && ... 
                        isfield(Sessions(i).sites(j).condition(cn).hs_tuned_evoked, 'mean')
                        for st = 1:size(Sessions(i).sites(j).condition(cn).hs_tuned_evoked, 1)
                            for hs = 1:size(Sessions(i).sites(j).condition(cn).hs_tuned_evoked, 2)
                                if isfield(Sessions(i).sites(j).condition(cn).hs_tuned_evoked(st, hs), 'mean') ...
                                        && ~isempty(Sessions(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).mean)
                                    sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites = ...
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites + 1;                                    
                                    
                                    if sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).nsites == 1%~isfield(sessions_avg.condition(cn).hs_tuned_evoked, 'mean')
                                        
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).hs_label ...
                                            = Sessions(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).hs_label;
                                        if isfield(Sessions(i).sites(j).condition(cn).hs_tuned_evoked(st, hs), 'state') && ...
                                                isfield(Sessions(i).sites(j).condition(cn).hs_tuned_evoked(st, hs), 'state_name')
                                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).state ...
                                                = Sessions(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).state;
                                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).state_name ...
                                                = Sessions(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).state_name;
                                        end
                                        
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).time ...
                                            = Sessions(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).time;
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).lfp ...
                                            = Sessions(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).mean;
                                        
                                    else
                                        
                                        ntimebins = length(sites_avg(t).condition(cn).hs_tuned_evoked(st, hs).time);
                                        % average same number of time bins
                                        if ntimebins > length(Sessions(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).time)
                                            ntimebins = length(Sessions(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).time);
                                        end
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).time = ...
                                            Sessions(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).time(1:ntimebins);
                                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).lfp ...
                                            = cat(1, ...
                                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).lfp(:, 1:ntimebins), ...
                                            (Sessions(i).sites(j).condition(cn).hs_tuned_evoked(st, hs).mean(1:ntimebins)));
                                                                   
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
                    if isfield(sites_avg(t).condition(cn).hs_tuned_evoked(st,hs), 'lfp') && ...
                            ~isempty(sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).lfp)
                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).dimord = 'nsites_time';                                
                        [sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).mean, ...
                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).error] = ...
                            lfp_tfa_compute_statistics(...
                            sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).lfp, lfp_tfa_cfg.error_measure);
                        sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).std = ...
                            nanstd(sites_avg(t).condition(cn).hs_tuned_evoked(st,hs).lfp, 0, 1);
                    end
                end
            end

            % plot
            if ~isempty(sites_avg(t).condition(cn).hs_tuned_evoked)
                if isfield(sites_avg(t).condition(cn).hs_tuned_evoked,... 
                        'mean')
                    plottitle = [lfp_tfa_cfg.monkey lfp_tfa_cfg.compare.targets{t},...
                        ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') ', ...
                        '_' lfp_tfa_cfg.conditions(cn).label];
                    result_file = fullfile(results_fldr, ...
                                    [lfp_tfa_cfg.monkey 'LFP_Evoked_' lfp_tfa_cfg.compare.targets{t} '_' lfp_tfa_cfg.conditions(cn).label ]);
                    lfp_tfa_plot_evoked_lfp(sites_avg(t).condition(cn).hs_tuned_evoked, ...
                                lfp_tfa_cfg, plottitle, result_file);
                end
            end            
        end
                % difference between conditions
        sites_avg(t).difference = [];
        for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
            diff_condition = lfp_tfa_cfg.diff_condition{diff};
            sites_avg(t).difference = [sites_avg(t).difference, ...
                lfp_tfa_compute_difference_condition_ev(sites_avg(t).condition, ...
                diff_condition, 1, lfp_tfa_cfg)];
        end
    end
    % save session average tfs
    save(fullfile(results_fldr, 'LFP_Evoked_sites_average.mat'), 'sites_avg');
    close all;
end