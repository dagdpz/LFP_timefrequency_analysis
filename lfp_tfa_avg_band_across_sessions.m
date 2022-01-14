function sessions_avg = lfp_tfa_avg_band_across_sessions(lfp_tfr, lfp_tfa_cfg)
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

% results folder
results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sessions', 'LFP_band');
if ~exist(results_fldr, 'dir')
    mkdir(results_fldr);
end

% Define frequency indices for band average
band_freq_index.gamma = find(lfp_tfa_cfg.band.gamma(1) < lfp_tfa_cfg.tfr.foi ...
    &  lfp_tfa_cfg.tfr.foi < lfp_tfa_cfg.band.gamma(2));
band_freq_index.beta = find(lfp_tfa_cfg.band.beta(1) < lfp_tfa_cfg.tfr.foi ...
    &  lfp_tfa_cfg.tfr.foi < lfp_tfa_cfg.band.beta(2));
band_freq_index.alpha = find(lfp_tfa_cfg.band.alpha(1) < lfp_tfa_cfg.tfr.foi ...
    &  lfp_tfa_cfg.tfr.foi < lfp_tfa_cfg.band.alpha(2));
band_freq_index.theta = find(lfp_tfa_cfg.band.theta(1) < lfp_tfa_cfg.tfr.foi ...
    &  lfp_tfa_cfg.tfr.foi < lfp_tfa_cfg.band.theta(2));

% Average TFR across sessions
sessions_avg = struct();
% conditions to compare
lfp_tfa_cfg.conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg);
for t = 1:length(lfp_tfa_cfg.compare.targets)
    sessions_avg(t).target = lfp_tfa_cfg.compare.targets{t};
    for cn = 1:length(lfp_tfa_cfg.conditions)
        fprintf('Condition %s\n', lfp_tfa_cfg.conditions(cn).label);
        sessions_avg(t).condition(cn).hs_tuned_band = struct();
        sessions_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
        sessions_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
        % initialize number of site pairs for each handspace
        % label
        for st = 1:size(lfp_tfa_cfg.analyse_states, 1)
            for hs = 1:size(lfp_tfa_cfg.conditions(1).hs_labels, 2)
                sessions_avg(t).condition(cn).hs_tuned_band(st, hs).nsessions = 0;
                sessions_avg(t).condition(cn).hs_tuned_band(st, hs).freq.powspctrm = [];
            end
        end
        %nsessions = 0;
        for i = 1:length(lfp_tfr.session)
            for k = 1:length(lfp_tfr.session(i).session_avg)
                if strcmp(lfp_tfr.session(i).session_avg(k).target, lfp_tfa_cfg.compare.targets{t})
                    if ~isempty(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band) && ...
                            isfield(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band, 'freq')
                        %nsessions = nsessions + 1;
                        for st = 1:size(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band, 1)
                            for hs = 1:size(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band, 2)
                                if isfield(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq, 'powspctrm') ...
                                        && ~isempty(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.powspctrm)
                                    sessions_avg(t).condition(cn).hs_tuned_band(st,hs).nsessions = ...
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).nsessions + 1;
                                    if sessions_avg(t).condition(cn).hs_tuned_band(st,hs).nsessions == 1
                                        
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.time ...
                                            = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.time;
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).hs_label ...
                                            = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).hs_label;
                                        if isfield(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs), 'state') ...
                                                && isfield(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs), 'state_name')
                                            sessions_avg(t).condition(cn).hs_tuned_band(st,hs).state ...
                                                = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).state;
                                            sessions_avg(t).condition(cn).hs_tuned_band(st,hs).state_name ...
                                                = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).state_name;
                                        end
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.cfg ...
                                            = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.cfg;
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.freq ...
                                            = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.freq;
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.powspctrm ...
                                            = nanmean(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.powspctrm, 1);
                                        
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.gamma ...
                                            = nanmean(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.gamma, 1);
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.beta ...
                                            = nanmean(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.beta, 1);
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.alpha ...
                                            = nanmean(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.alpha, 1);
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.theta ...
                                            = nanmean(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.theta, 1);
                                        if isfield(sessions_avg(t).condition(cn).hs_tuned_band(st,hs), 'nsites')
                                            sessions_avg(t).condition(cn).hs_tuned_band(st,hs).nsites ...
                                                = sessions_avg(t).condition(cn).hs_tuned_band(st,hs).nsites;
                                        end
                                    else
                                        ntimebins = size(sessions_avg(t).condition(cn).hs_tuned_band(st, hs).freq.powspctrm, 3);
                                        % average same number of time bins
                                        if ntimebins > length(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.time)
                                            ntimebins = length(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.time);
                                        end
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.powspctrm ...
                                            = cat(1, sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.powspctrm(:,:,1:ntimebins), ...
                                            nanmean(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.powspctrm(:,:,1:ntimebins), 1));
                                        
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.gamma ...
                                            = cat(1, sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.gamma(:,:,1:ntimebins), ...
                                            nanmean(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.gamma(:,:,1:ntimebins), 1));
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.beta ...
                                            = cat(1, sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.beta(:,:,1:ntimebins), ...
                                            nanmean(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.beta(:,:,1:ntimebins), 1));
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.alpha ...
                                            = cat(1, sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.alpha(:,:,1:ntimebins), ...
                                            nanmean(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.alpha(:,:,1:ntimebins), 1));
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.theta ...
                                            = cat(1, sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.theta(:,:,1:ntimebins), ...
                                            nanmean(lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).freq.theta(:,:,1:ntimebins), 1));
                                        
                                        sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.time = ...
                                            sessions_avg(t).condition(cn).hs_tuned_band(st,hs).freq.time(1:ntimebins);
                                        if isfield(sessions_avg(t).condition(cn).hs_tuned_band(st,hs), 'nsites')
                                            sessions_avg(t).condition(cn).hs_tuned_band(st,hs).nsites ...
                                                = lfp_tfr.session(i).session_avg(k).condition(cn).hs_tuned_band(st, hs).nsites + ...
                                                sessions_avg(t).condition(cn).hs_tuned_band(st,hs).nsites;
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
        
        
        
        if ~isempty(sessions_avg(t).condition(cn).hs_tuned_band)
            if isfield(sessions_avg(t).condition(cn).hs_tuned_band,...
                    'freq')
                plottitle = ['Target ' lfp_tfa_cfg.compare.targets{t}, ...
                    ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') ', ...
                    sessions_avg(t).condition(cn).label];
                result_file = fullfile(results_fldr, ...
                    ['LFP_band_' lfp_tfa_cfg.compare.targets{t}, '_', ...
                    sessions_avg(t).condition(cn).label]);
                lfp_tfa_plot_hs_tuned_tfr_band_average(sessions_avg(t).condition(cn).hs_tuned_band, ...
                    lfp_tfa_cfg, plottitle, result_file);
            end
        end
    end
    for trial_type = 1:2
        
        plottitle = ['Target ' lfp_tfa_cfg.compare.targets{t}, ...
            ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') '];
        if trial_type == 1
            plottitle = [plottitle 'Instructed'];
            trial_title = 'Instructed';
        elseif trial_type == 2
            plottitle = [plottitle 'Choice'];
            trial_title = 'Choice';
        end
        result_file = fullfile(results_fldr, ...
            ['LFP_band_' lfp_tfa_cfg.compare.targets{t}, '_combined_', ...
            trial_title]);
        lfp_tfa_plot_hs_tuned_tfr_band_average_combined(sessions_avg(t).condition, ...
            lfp_tfa_cfg, plottitle, result_file, trial_type);
    end
end







% save session average tfs
save(fullfile(results_fldr, 'LFP_band_sessions_avg.mat'), 'sessions_avg');
end