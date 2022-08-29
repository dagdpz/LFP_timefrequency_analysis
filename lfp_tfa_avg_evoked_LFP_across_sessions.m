function sessions_avg = lfp_tfa_avg_evoked_LFP_across_sessions(Sessions, lfp_tfa_cfg)
%lfp_tfa_avg_evoked_LFP_across_sessions  - Condition-based evoked LFP response
% grand average of averages across all sites within a session
%
% USAGE:
%	sessions_avg = lfp_tfa_avg_evoked_LFP_across_sessions(lfp_evoked, lfp_tfa_cfg)
%
% INPUTS:
%		lfp_evoked		- struct containing the condition-based evoked LFP
%                         response average across all sites of each session
%                         analysed, i.e., the output of
%                         lfp_tfa_plot_site_evoked_LFP.m
%           Required Fields:
%               session.session_avg - 1xN struct containing condition-based
%               average evoked LFP response for N sessions (session_avg =
%               Average of site averages for one session)
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields:
%               1. conditions          - trial conditions to compare, see
%               lfp_tfa_settings.m and lfp_tfa_compare_conditions.m
%               2. root_results_fldr   - root folder where results are
%               saved. Results will be saved under
%               [lfp_tfa_cfg.root_results_fldr ...
%               '/Avg_across_sessions/LFP_Evoked']
%               3. compare.targets     - targets to compare, see lfp_tfa_settings.m
%               4. ref_hemisphere      - reference hemisphere for contra-
%               and ipsi- labelling, see settings/lfp_tfa_settings_example.m
% OUTPUTS:
%		sessions_avg    - structure containing condition-based evoked LFP
%		response averaged across multiple sessions
%
% REQUIRES:	lfp_tfa_plot_evoked_lfp
%
% See also lfp_tfa_define_settings, lfp_tfa_compare_conditions,
% lfp_tfa_plot_site_evoked_LFP
%
% Author(s):	S.Nair, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2019-02-15:	Created function (Sarath Nair)
% 2019-03-05:	First Revision
% 2019-10-22:   Edited documentation
% ...
% $Revision: 1.1 $  $Date: 2019-10-22 11:03:00 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%


% results folder
results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sessions', 'LFP_Evoked');
if ~exist(results_fldr, 'dir')
    mkdir(results_fldr);
end

% Average Evoked LFP response across sessions
sessions_avg = struct();

for t = 1:length(lfp_tfa_cfg.compare.targets)
    sessions_avg(t).target = lfp_tfa_cfg.compare.targets{t};
    for cn = 1:length(lfp_tfa_cfg.conditions)
        fprintf('Condition %s\n', lfp_tfa_cfg.conditions(cn).label);
        sessions_avg(t).condition(cn).avg_across_sessions = struct();
        sessions_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
        sessions_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
        % initialize number of site pairs for each handspace
        % label
        for st = 1:size(Sessions(1).session_avg(1).condition(cn).hs_tuned_evoked, 1)
            for hs = 1:size(Sessions(1).session_avg(1).condition(cn).hs_tuned_evoked, 2)
                sessions_avg(t).condition(cn).avg_across_sessions(st, hs).nsessions = 0;
                sessions_avg(t).condition(cn).avg_across_sessions(st, hs).lfp = [];
            end
        end
        
        for i = 1:length(Sessions)
            for k = 1:length(Sessions(i).session_avg)
                %% LS 2021: question is !? combine hemispheres BEFORE per session averaging??
%                 if ismember(lfp_tfa_cfg.compare.targets{t}, Sessions(i).session_avg(k).target)
                    %if strcmp(Sessions(i).session_avg(k).target, lfp_tfa_cfg.compare.targets{t})
                    
                    if ~isempty(Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked) && ...
                            isfield(Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked, 'mean')
                        for st = 1:size(Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked, 1)
                            for hs = 1:size(Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked, 2)
                                if isfield(Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs), 'mean') ...
                                        && ~isempty(Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).mean)
                                    
                                    sessions_avg(t).condition(cn).avg_across_sessions(st, hs).nsessions = ...
                                        sessions_avg(t).condition(cn).avg_across_sessions(st, hs).nsessions + 1;
                                    
                                    if sessions_avg(t).condition(cn).avg_across_sessions(st, hs).nsessions == 1%~isfield(sessions_avg(t).condition(cn).avg_across_sessions, 'mean')
                                        sessions_avg(t).condition(cn).avg_across_sessions(st,hs).time ...
                                            = Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).time;
                                        sessions_avg(t).condition(cn).avg_across_sessions(st,hs).hs_label ...
                                            = Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).hs_label;
                                        if isfield(Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs), 'state') && ...
                                                isfield(Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs), 'state_name')
                                            sessions_avg(t).condition(cn).avg_across_sessions(st,hs).state ...
                                                = Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).state;
                                            sessions_avg(t).condition(cn).avg_across_sessions(st,hs).state_name ...
                                                = Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).state_name;
                                        end
                                        sessions_avg(t).condition(cn).avg_across_sessions(st,hs).lfp ...
                                            = Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).mean;
                                        if isfield(Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs), 'nsites')
                                            sessions_avg(t).condition(cn).avg_across_sessions(st,hs).nsites ...
                                                = Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).nsites;
                                        end
                                    else
                                        ntimebins = length(sessions_avg(t).condition(cn).avg_across_sessions(st, hs).time);
                                        % average same number of time bins
                                        if ntimebins > length(Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).time)
                                            ntimebins = length(Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).time);
                                        end
                                        sessions_avg(t).condition(cn).avg_across_sessions(st,hs).time = ...
                                            sessions_avg(t).condition(cn).avg_across_sessions(st,hs).time(1:ntimebins);
                                        sessions_avg(t).condition(cn).avg_across_sessions(st,hs).lfp ...
                                            = cat(1, ...
                                            sessions_avg(t).condition(cn).avg_across_sessions(st,hs).lfp(:, 1:ntimebins), ...
                                            Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).mean(1:ntimebins));
                                        if isfield(sessions_avg(t).condition(cn).avg_across_sessions(st,hs), 'nsites')
                                            sessions_avg(t).condition(cn).avg_across_sessions(st,hs).nsites ...
                                                = Sessions(i).session_avg(k).condition(cn).hs_tuned_evoked(st, hs).nsites + ...
                                                sessions_avg(t).condition(cn).avg_across_sessions(st,hs).nsites;
                                        end
                                    end
                                end
                            end
                        end
                    end
%                 else
%                     continue;
%                 end
            end
        end
        
        % compute average
        if isfield(sessions_avg(t).condition(cn).avg_across_sessions, 'lfp')
            for st = 1:size(sessions_avg(t).condition(cn).avg_across_sessions, 1)
                for hs = 1:size(sessions_avg(t).condition(cn).avg_across_sessions, 2)
                    if ~isempty(sessions_avg(t).condition(cn).avg_across_sessions(st,hs).lfp)
                        sessions_avg(t).condition(cn).avg_across_sessions(st,hs).dimord = 'nsessions_time';
                        [sessions_avg(t).condition(cn).avg_across_sessions(st,hs).mean, ...
                            sessions_avg(t).condition(cn).avg_across_sessions(st,hs).error] = ...
                            lfp_tfa_compute_statistics(...
                            sessions_avg(t).condition(cn).avg_across_sessions(st,hs).lfp);
                        %                             sessions_avg(t).condition(cn).avg_across_sessions(st,hs).std = ...
                        %                                 nanstd(sessions_avg(t).condition(cn).avg_across_sessions(st,hs).lfp, 0, 1);
                    end
                end
            end
        end
        sessions_avg(t).difference = [];
        for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
            diff_condition = lfp_tfa_cfg.diff_condition{diff};
            
            sessions_avg(t).difference = [sessions_avg(t).difference, ...
                lfp_tfa_compute_difference_condition_ev(sessions_avg(t).condition, diff_condition,1,lfp_tfa_cfg)];
            
%             if ~isempty(sessions_avg(t).condition(cn).avg_across_sessions)
%                 if isfield(sessions_avg(t).condition(cn).avg_across_sessions,...
%                         'mean')
%                     plottitle = [lfp_tfa_cfg.monkey lfp_tfa_cfg.compare.targets{t} ...
%                         lfp_tfa_cfg.conditions(cn).label];
%                     result_file = fullfile(results_fldr, ...
%                         [lfp_tfa_cfg.monkey 'LFP_Evoked_' sessions_avg(t).target '_' ...
%                         lfp_tfa_cfg.conditions(cn).label]);
%                     lfp_tfa_plot_evoked_lfp(sessions_avg(t).condition(cn).avg_across_sessions, ...
%                         lfp_tfa_cfg, plottitle, result_file);
%                 end
%             end
            % save session average tfs
            save(fullfile(results_fldr, [lfp_tfa_cfg.monkey 'sessions_average_evoked.mat']), 'sessions_avg');
        end
    end
    close all;
    
end