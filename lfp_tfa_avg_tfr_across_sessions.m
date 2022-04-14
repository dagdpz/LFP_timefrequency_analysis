function sessions_avg = lfp_tfa_avg_tfr_across_sessions(Sessions, lfp_tfa_cfg,varargin)
%lfp_tfa_avg_tfr_across_sessions  - Condition-based LFP time frequency
%response average across many session averages (A session average is
%the LFP TFR average across site averages recorded in a session. A site
%the LFP TFR average across multiple trials recorded at a site in a session)
%
% USAGE:
%	sessions_avg = lfp_tfa_avg_tfr_across_sessions(lfp_tfr, lfp_tfa_cfg)
%
% INPUTS:
%		lfp_tfr     	- struct containing the condition-based LFP time freq spectrogram for
%		indiviual sites, i.e., the output of lfp_tfa_plot_site_average_tfr.m
%           Required Fields:
%               session.session_avg - session is a 1xM struct (M is the
%               number of sessions) and session_avg is a 1xK struct (K is
%               the number of target areas) containing average LFP TFR
%               results for a session
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields:
%               1. conditions          - trial conditions to compare, see
%               lfp_tfa_settings.m and lfp_tfa_compare_conditions.m
%               2. root_results_fldr   - root folder where results are
%               saved. Results will be saved under
%               [lfp_tfa_cfg.root_results_fldr ...
%               '/Avg_across_sessions/LFP_TFR']
%               3. compare.targets     - targets to compare, see lfp_tfa_settings.m
%               4. ref_hemisphere      - reference hemisphere for contra-
%               and ipsi- labelling, see settings/lfp_tfa_settings_example.m
% OUTPUTS:
%		sessions_avg    - structure containing condition-based LFP
%		spectrogram response averaged across multiple session averages
%
% REQUIRES:	lfp_tfa_plot_hs_tuned_tfr_multiple_img,
% lfp_tfa_compute_difference_condition_tfr
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_define_settings,
% lfp_tfa_compare_conditions, lfp_tfa_plot_site_average_tfr,
% lfp_tfa_compute_difference_condition_tfr,
% lfp_tfa_plot_hs_tuned_tfr_multiple_img, lfp_tfa_avg_tfr_across_sites
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

% results folder % results folder
    if nargin < 3
results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sessions');
    else
        results_fldr = varargin{1};
    end
if ~exist(results_fldr, 'dir')
    mkdir(results_fldr);
end

% conditions to compare
sessions_avg = struct(); % why was this inside the loop?
lfp_tfa_cfg.conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg);
for t = 1:length(lfp_tfa_cfg.compare.targets)
% Average TFR across sessions
    sessions_avg(t).target = lfp_tfa_cfg.compare.targets{t};
    for cn = 1:length(lfp_tfa_cfg.conditions)
        fprintf('Condition %s\n', lfp_tfa_cfg.conditions(cn).label);
        sessions_avg(t).condition(cn).hs_tuned_tfs = struct();
        sessions_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
        sessions_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
        % initialize number of site pairs for each handspace
        % label
        for st = 1:size(lfp_tfa_cfg.analyse_states, 1)
            for hs = 1:size(lfp_tfa_cfg.conditions(1).hs_labels, 2)
                sessions_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsessions = 0;
                sessions_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm = [];
            end
        end
        %nsessions = 0;
        for i = 1:length(Sessions)
            for k = 1:length(Sessions(i).session_avg)
                %% LS 2021: question is !? combine hemispheres BEFORE per session averaging??
                if ismember(lfp_tfa_cfg.compare.targets{t}, Sessions(i).session_avg(k).target)
                %if strcmp(Sessions(i).session_avg(k).target, lfp_tfa_cfg.compare.targets{t})
                    if ~isempty(Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs) && ...
                            isfield(Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs, 'freq')
                        %nsessions = nsessions + 1;
                        for st = 1:size(Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs, 1)
                            for hs = 1:size(Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs, 2)
                                if isfield(Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).freq, 'powspctrm') ...
                                        && ~isempty(Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm)
                                    sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsessions = ...
                                        sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsessions + 1;
                                    if sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsessions == 1
                                        
                                        sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).freq.time ...
                                            = Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).freq.time;
                                        sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).hs_label ...
                                            = Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).hs_label;
                                        if isfield(Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs), 'state') ...
                                                && isfield(Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs), 'state_name')
                                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).state ...
                                                = Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).state;
                                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).state_name ...
                                                = Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).state_name;
                                        end
                                        sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).freq.cfg ...
                                            = Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).freq.cfg;
                                        sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).freq.freq ...
                                            = Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).freq.freq;
                                        sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).freq.powspctrm ...
                                            = nanmean(Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm, 1);
                                        if isfield(sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs), 'nsites')
                                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsites ...
                                                = sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsites;
                                        end
                                    else
                                        ntimebins = size(sessions_avg(t).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm, 3);
                                        % average same number of time bins
                                        if ntimebins > length(Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).freq.time)
                                            ntimebins = length(Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).freq.time);
                                        end
                                        sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).freq.powspctrm ...
                                            = cat(1, sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).freq.powspctrm(:,:,1:ntimebins), ...
                                            nanmean(Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).freq.powspctrm(:,:,1:ntimebins), 1));
                                        sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).freq.time = ...
                                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).freq.time(1:ntimebins);
                                        if isfield(sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs), 'nsites')
                                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsites ...
                                                = Sessions(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).nsites + ...
                                                sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsites;
                                        end
                                    end
                                end
                            end
                        end
                    end
                else
                    continue;
                end
            end
        end
        
        % compute average
        %             if isfield(sessions_avg(t).condition(cn).hs_tuned_tfs, 'powspctrm')
        %                 for st = 1:size(sessions_avg(t).condition(cn).hs_tuned_tfs, 1)
        %                     for hs = 1:size(sessions_avg(t).condition(cn).hs_tuned_tfs, 2)
        %                         if ~isempty(sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm)
        %                             sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsessions = nsessions;
        %                             sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm = ...
        %                                 (1/nsessions) * sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm;
        %                         end
        %                     end
        %                 end
        %             end
        
        
        if ~isempty(sessions_avg(t).condition(cn).hs_tuned_tfs)
            if isfield(sessions_avg(t).condition(cn).hs_tuned_tfs,...
                    'freq')
                plottitle = [lfp_tfa_cfg.monkey 'Target ' lfp_tfa_cfg.compare.targets{t}, ...
                    ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') ', ...
                    sessions_avg(t).condition(cn).label];
                result_file = fullfile(results_fldr, ...
                    [lfp_tfa_cfg.monkey '_LFP_TFR_' lfp_tfa_cfg.compare.targets{t}, '_', ...
                    sessions_avg(t).condition(cn).label]);
                lfp_tfa_plot_hs_tuned_tfr_multiple_img(sessions_avg(t).condition(cn).hs_tuned_tfs, ...
                    lfp_tfa_cfg, plottitle, result_file);
            end
        end
        
    end
    % Difference TFR
    % difference between pre- and post-injection
    %sessions_avg(t).difference = lfp_tfa_compute_diff_condition_tfr(sessions_avg(t), lfp_tfa_cfg.difference);
    sessions_avg(t).difference = [];
    for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
        diff_condition = lfp_tfa_cfg.diff_condition{diff};
        
        sessions_avg(t).difference = [sessions_avg(t).difference, ...
            lfp_tfa_compute_difference_condition_tfr(sessions_avg(t).condition, diff_condition,1,lfp_tfa_cfg)];
        
    end
    % plot Difference TFR
    for dcn = 1:length(sessions_avg(t).difference)
        if ~isempty(sessions_avg(t).difference(dcn).hs_tuned_tfs)
            if isfield(sessions_avg(t).difference(dcn).hs_tuned_tfs,...
                    'freq')
                plottitle = [lfp_tfa_cfg.monkey 'Target' lfp_tfa_cfg.compare.targets{t}, ...
                    ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') ', ...
                    sessions_avg(t).difference(dcn).label];
                result_file = fullfile(results_fldr, ...
                    [lfp_tfa_cfg.monkey 'LFP_DiffTFR_' lfp_tfa_cfg.compare.targets{t} ...
                    '_' sessions_avg(t).difference(dcn).label ...
                    '_' sessions_avg(t).difference(dcn).cfg_condition.diff...
                    '_c_' num2str(sessions_avg(t).difference(dcn).cfg_condition.choice) '_p_' ...
                    num2str(sessions_avg(t).difference(dcn).cfg_condition.perturbation)]);
                if lfp_tfa_cfg.plot_significant
                    plottitle = [lfp_tfa_cfg.monkey 'Target ', lfp_tfa_cfg.compare.targets{t}, ...
                        ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') ', ...
                        sessions_avg(t).difference(dcn).label];
                    result_file = fullfile(results_fldr, ...
                        [lfp_tfa_cfg.monkey 'LFP_DiffTFR_' lfp_tfa_cfg.compare.targets{t} ...
                        '_' sessions_avg(t).difference(dcn).label ...
                        '_' sessions_avg(t).difference(dcn).cfg_condition.diff ...
                         '_c_' num2str(sessions_avg(t).difference(dcn).cfg_condition.choice) '_p_' ...
                         num2str(sessions_avg(t).difference(dcn).cfg_condition.perturbation) ...
                        '_' lfp_tfa_cfg.correction_method]);
                end
                %sessions_avg(t).difference(dcn).label '.png']);
                if lfp_tfa_cfg.plot_significant
                    lfp_tfa_plot_hs_tuned_tfr_multiple_img(sessions_avg(t).difference(dcn).hs_tuned_tfs, ...
                        lfp_tfa_cfg, plottitle, result_file, 'bluewhitered',1);
                else
                    lfp_tfa_plot_hs_tuned_tfr_multiple_img(sessions_avg(t).difference(dcn).hs_tuned_tfs, ...
                        lfp_tfa_cfg, plottitle, result_file, 'bluewhitered');
                end
            end
        end
    end
    
end

% save session average tfs
save(fullfile(results_fldr, [lfp_tfa_cfg.monkey  '_LFP_TFR_sessions_avg.mat']), 'sessions_avg');
end