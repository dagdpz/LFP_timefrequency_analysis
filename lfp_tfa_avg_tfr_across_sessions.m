function sessions_avg = lfp_tfa_avg_tfr_across_sessions(lfp_tfr, lfp_tfa_cfg)
%lfp_tfa_avg_tfr_across_sessions  - Condition-based LFP time frequency
%response average across many session averages
%
% USAGE:
%	sites_avg = lfp_tfa_avg_tfr_across_sessions(lfp_tfr, lfp_tfa_cfg)
%
% INPUTS:
%		lfp_tfr     	- struct containing the condition-based LFP time freq spectrogram for
%		indiviual sites, output of lfp_tfa_plot_site_average_tfr.m
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields:
%               1. conditions          - trial conditions to compare, see
%               lfp_tfa_settings.m and lfp_tfa_compare_conditions.m
%               2. root_results_fldr   - root folder where results are saved
%               3. compare.targets     - targets to compare, see lfp_tfa_settings.m
%               4. 
% OUTPUTS:
%		sessions_avg    - structure containing condition-based LFP
%		spectrogram response averaged across multiple sessions
%
% REQUIRES:	lfp_tfa_plot_hs_tuned_tfr, lfp_tfa_compute_diff_tfr
%
% See also lfp_tfa_settings, lfp_tfa_define_settings, lfp_tfa_compare_conditions, 
% lfp_tfa_plot_site_average_tfr, lfp_tfa_compute_diff_tfr
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
    results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sessions', 'LFP_TFR');
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % Average TFR across sessions
    sessions_avg = struct();
    % conditions to compare
    lfp_tfa_cfg.conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg);
    for t = 1:length(lfp_tfa_cfg.compare.targets)
        sessions_avg(t).target = lfp_tfa_cfg.compare.targets{t};
        for cn = 1:length(lfp_tfa_cfg.conditions)
            fprintf('Condition %s\n', lfp_tfa_cfg.conditions(cn).label);
            sessions_avg(t).condition(cn).hs_tuned_tfs = struct();
            sessions_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
            sessions_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
            nsessions = 0;
            for i = 1:length(lfp_tfr.session)  
                for k = 1:length(lfp_tfr.session(i).session_avg)
                    if strcmp(lfp_tfr.session(i).session_avg(k).target, lfp_tfa_cfg.compare.targets{t})
                        if ~isempty(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs) && ... 
                            isfield(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs, 'powspctrm')
                            nsessions = nsessions + 1;   
                            for st = 1:size(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs, 1)
                                for hs = 1:size(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs, 2)
                                    if isfield(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs), 'powspctrm') ...
                                            && ~isempty(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).powspctrm)
                                        if nsessions == 1
                                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).time ...
                                            = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).time;
                                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).hs_label ...
                                                = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).hs_label;
                                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).state ...
                                                = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).state;
                                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).state_name ...
                                                = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).state_name;
                                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).cfg ...
                                                = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).cfg;
                                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).freq ...
                                                = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).freq;
                                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm ...
                                                = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).powspctrm;
                                            if isfield(sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs), 'nsites')
                                                sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsites ...
                                                    = sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsites;
                                            end
                                        else
                                            ntimebins = size(sessions_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm, 3);
                                            % average same number of time bins
                                            if ntimebins > length(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).time)
                                                ntimebins = length(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).time);
                                            end
                                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm ...
                                                = (lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins)) + ...
                                                sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm(1,:,1:ntimebins);
                                            if isfield(sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs), 'nsites')
                                                sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsites ...
                                                    = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs).nsites + ...
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
            if isfield(sessions_avg(t).condition(cn).hs_tuned_tfs, 'powspctrm')
                for st = 1:size(sessions_avg(t).condition(cn).hs_tuned_tfs, 1)
                    for hs = 1:size(sessions_avg(t).condition(cn).hs_tuned_tfs, 2)
                        if ~isempty(sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm)
                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).nsessions = nsessions;                                
                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm = ...
                                (1/nsessions) * sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm;
                        end
                    end
                end
            end


            if ~isempty(sessions_avg(t).condition(cn).hs_tuned_tfs)
                if isfield(sessions_avg(t).condition(cn).hs_tuned_tfs,... 
                        'powspctrm')
                    plottitle = ['Target ' lfp_tfa_cfg.compare.targets{t}, ...
                        ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') ', ...
                        sessions_avg(t).condition(cn).label];
                    result_file = fullfile(results_fldr, ...
                        ['LFP_TFR_' lfp_tfa_cfg.compare.targets{t}, '_', ...
                        sessions_avg(t).condition(cn).label '.png']);
                    lfp_tfa_plot_hs_tuned_tfr_multiple_img(sessions_avg(t).condition(cn).hs_tuned_tfs, ...
                        lfp_tfa_cfg, plottitle, result_file);
                end
            end

        end
        % Difference TFR
        % difference between pre- and post-injection
        %sessions_avg(t).difference = lfp_tfa_compute_diff_condition_tfr(sessions_avg(t), lfp_tfa_cfg.difference);
        sessions_avg(t).difference = [];
        for diff = 1:size(lfp_tfa_cfg.difference, 1)
            diff_field = lfp_tfa_cfg.difference{diff, 1};
            diff_values = lfp_tfa_cfg.difference{diff, 2};
            % check if both pre- and post- injection blocks exist
            if strcmp(diff_field, 'perturbation')
                if sum(lfp_tfa_cfg.compare.perturbations == [diff_values{:}]) <= 1
                    continue;                        
                end
            elseif strcmp(diff_field, 'choice')
                if sum(lfp_tfa_cfg.compare.choice_trials == [diff_values{:}]) <= 1
                    continue;
                end
            end
            sessions_avg(t).difference = [sessions_avg(t).difference, ...
                lfp_tfa_compute_diff_condition_tfr(sessions_avg(t), diff_field, diff_values)];
        end
        % plot Difference TFR
        for dcn = 1:length(sessions_avg(t).difference)
            if ~isempty(sessions_avg(t).difference(dcn).hs_tuned_tfs)
                if isfield(sessions_avg(t).difference(dcn).hs_tuned_tfs,... 
                        'powspctrm')
                    plottitle = ['Target' lfp_tfa_cfg.compare.targets{t}, ...
                        ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') ', ...
                        sessions_avg(t).difference(dcn).label];
                    result_file = fullfile(results_fldr, ...
                                ['LFP_DiffTFR_' lfp_tfa_cfg.compare.targets{t} ...
                                '_' sessions_avg(t).difference(dcn).label '.png']);
                    lfp_tfa_plot_hs_tuned_tfr_multiple_img(sessions_avg(t).difference(dcn).hs_tuned_tfs, ...
                                lfp_tfa_cfg, plottitle, result_file, 'bluewhitered');
                end
            end
        end
        
    end
    
    % save session average tfs
    save(fullfile(results_fldr, 'LFP_TFR_sessions_avg.mat'), 'sessions_avg');
end