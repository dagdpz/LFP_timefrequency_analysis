function sites_avg = lfp_tfa_avg_evoked_LFP_across_sites(Input, lfp_tfa_cfg, varargin)
%lfp_tfa_avg_evoked_LFP_across_sites  - Condition-based evoked LFP response
% grand average across many site averages
%
% USAGE:
%	sites_avg = lfp_tfa_avg_evoked_LFP_across_sites(Input, lfp_tfa_cfg)
%
% INPUTS:
%		Input		- struct containing the condition-based evoked LFP 
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
%               '/Avg_across_sites/Input']
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
        results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sites', 'Input');
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
            sites_avg(t).condition(cn).tuned = struct();
            sites_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
            sites_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
            % initialize number of site pairs for each handspace
            % label
            for st = 1:size(Input.session(1).sites(1).condition(cn).tuned, 1)
                for hs = 1:size(Input.session(1).sites(1).condition(cn).tuned, 2)
                    sites_avg(t).condition(cn).tuned(st, hs).nsites = 0;
                    sites_avg(t).condition(cn).tuned(st, hs).lfp = [];
                end
            end
            for i = 1:length(Input.session) 
                for j = 1:length(Input.session(i).sites)
                    %% LS 2021 
                    if ~ismember(lfp_tfa_cfg.compare.targets{t}, Input.session(i).sites(j).target) 
                        continue;
                    end
                    if ~Input.session(i).sites(j).use_for_avg
                        continue;
                    end
                    if ~isempty(Input.session(i).sites(j).condition(cn).tuned) && ... 
                        isfield(Input.session(i).sites(j).condition(cn).tuned, 'mean')
                        for st = 1:size(Input.session(i).sites(j).condition(cn).tuned, 1)
                            for hs = 1:size(Input.session(i).sites(j).condition(cn).tuned, 2)
                                if isfield(Input.session(i).sites(j).condition(cn).tuned(st, hs), 'mean') ...
                                        && ~isempty(Input.session(i).sites(j).condition(cn).tuned(st, hs).mean)
                                    sites_avg(t).condition(cn).tuned(st, hs).nsites = ...
                                        sites_avg(t).condition(cn).tuned(st, hs).nsites + 1;                                    
                                    
                                    if sites_avg(t).condition(cn).tuned(st, hs).nsites == 1%~isfield(sessions_avg.condition(cn).tuned, 'mean')
                                        
                                        sites_avg(t).condition(cn).tuned(st,hs).hs_label ...
                                            = Input.session(i).sites(j).condition(cn).tuned(st, hs).hs_label;
                                        if isfield(Input.session(i).sites(j).condition(cn).tuned(st, hs), 'state') && ...
                                                isfield(Input.session(i).sites(j).condition(cn).tuned(st, hs), 'state_name')
                                            sites_avg(t).condition(cn).tuned(st,hs).state ...
                                                = Input.session(i).sites(j).condition(cn).tuned(st, hs).state;
                                            sites_avg(t).condition(cn).tuned(st,hs).state_name ...
                                                = Input.session(i).sites(j).condition(cn).tuned(st, hs).state_name;
                                        end
                                        
                                        sites_avg(t).condition(cn).tuned(st,hs).time ...
                                            = Input.session(i).sites(j).condition(cn).tuned(st, hs).time;
                                        sites_avg(t).condition(cn).tuned(st,hs).lfp ...
                                            = Input.session(i).sites(j).condition(cn).tuned(st, hs).mean;
                                        
                                    else
                                        
                                        ntimebins = length(sites_avg(t).condition(cn).tuned(st, hs).time);
                                        % average same number of time bins
                                        if ntimebins > length(Input.session(i).sites(j).condition(cn).tuned(st, hs).time)
                                            ntimebins = length(Input.session(i).sites(j).condition(cn).tuned(st, hs).time);
                                        end
                                        sites_avg(t).condition(cn).tuned(st,hs).time = ...
                                            Input.session(i).sites(j).condition(cn).tuned(st, hs).time(1:ntimebins);
                                        sites_avg(t).condition(cn).tuned(st,hs).lfp ...
                                            = cat(1, ...
                                            sites_avg(t).condition(cn).tuned(st,hs).lfp(:, 1:ntimebins), ...
                                            (Input.session(i).sites(j).condition(cn).tuned(st, hs).mean(1:ntimebins)));
                                                                   
                                    end
                                end
                            end
                        end
                    end
                end
            end

            % compute average
            for st = 1:size(sites_avg(t).condition(cn).tuned, 1)
                for hs = 1:size(sites_avg(t).condition(cn).tuned, 2)
                    if isfield(sites_avg(t).condition(cn).tuned(st,hs), 'lfp') && ...
                            ~isempty(sites_avg(t).condition(cn).tuned(st,hs).lfp)
                        sites_avg(t).condition(cn).tuned(st,hs).dimord = 'nsites_time';                                
                        [sites_avg(t).condition(cn).tuned(st,hs).mean, ...
                            sites_avg(t).condition(cn).tuned(st,hs).error] = ...
                            lfp_tfa_compute_statistics(...
                            sites_avg(t).condition(cn).tuned(st,hs).lfp, lfp_tfa_cfg.error_measure);
                        sites_avg(t).condition(cn).tuned(st,hs).std = ...
                            nanstd(sites_avg(t).condition(cn).tuned(st,hs).lfp, 0, 1);
                    end
                end
            end

            % plot
            if ~isempty(sites_avg(t).condition(cn).tuned)
                if isfield(sites_avg(t).condition(cn).tuned,... 
                        'mean')
                    plottitle = [lfp_tfa_cfg.compare.targets{t},...
                        ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') ', ...
                        '_' lfp_tfa_cfg.conditions(cn).label];
                    result_file = fullfile(results_fldr, ...
                                    ['Input_' lfp_tfa_cfg.compare.targets{t} '_' lfp_tfa_cfg.conditions(cn).label ]);
                    lfp_tfa_plot_evoked_lfp(sites_avg(t).condition(cn).tuned, ...
                                lfp_tfa_cfg, plottitle, result_file);
                end
            end            
        end
    end
    % save session average tfs
    save(fullfile(results_fldr, 'Input_sites_average.mat'), 'sites_avg');
    close all;
end