function sites_avg = lfp_tfa_avg_band_across_sites(lfp_tfr, lfp_tfa_cfg, varargin)
%lfp_tfa_avg_tfr_across_sites  - Condition-based LFP time frequency
%response average across many site averages (A site average is the LFP TFR
%average across multiple trials recorded at a site in a session)
%
% USAGE:
%	sites_avg = lfp_tfa_avg_tfr_across_sites(lfp_tfr, lfp_tfa_cfg)
%
% INPUTS:
%		lfp_tfr     	- struct containing the condition-based LFP time freq spectrogram average for
%		indiviual sites and sessions, i.e., the output of lfp_tfa_plot_site_average_tfr.m
%           Required Fields:
%               session.sites - session is a 1xM struct (M is the
%               number of sessions) and sites is a 1xN struct (N is
%               the number of sites) containing LFP TFR average across
%               sites for a session
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields:
%               1. conditions          - trial conditions to compare, see
%               settings/lfp_tfa_settings_example and lfp_tfa_compare_conditions.m
%               2. root_results_fldr   - root folder where results are
%               saved. Results will be saved under
%               [lfp_tfa_cfg.root_results_fldr ...
%               '/Avg_across_sites/LFP_TFR']
%               3. compare.targets     - targets to compare, see lfp_tfa_settings.m
%               4. ref_hemisphere      - reference hemisphere for contra-
%               and ipsi- labelling, see settings/lfp_tfa_settings_example.m
%               5. analyse_states      - time windows around the states to
%               be analysed, see settings/lfp_tfa_settings_example
% OUTPUTS:
%		sites_avg    - structure containing condition-based LFP
%		spectrogram response averaged across multiple sites of multiple
%		sessions
%
% REQUIRES:	lfp_tfa_plot_hs_tuned_band_multiple_img,
% lfp_tfa_compute_difference_condition_tfr
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_define_settings,
% lfp_tfa_compare_conditions, lfp_tfa_plot_site_average_tfr,
% lfp_tfa_compute_difference_condition_tfr,
% lfp_tfa_plot_hs_tuned_band_multiple_img, lfp_tfa_avg_tfr_across_sessions
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
    results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sites', 'LFP_band');
else
    results_fldr = varargin{1};
end
if ~exist(results_fldr, 'dir')
    mkdir(results_fldr);
end

% defaults
plot_significant = 0;
if isfield(lfp_tfa_cfg, 'plot_significant')
    plot_significant = lfp_tfa_cfg.plot_significant;
end

% Average TFR across sites
sites_avg = struct();
for t = 1:length(lfp_tfa_cfg.compare.targets)
    sites_avg(t).target = lfp_tfa_cfg.compare.targets{t};
    for cn = 1:length(lfp_tfa_cfg.conditions)
        fprintf('Condition %s\n', lfp_tfa_cfg.conditions(cn).label);
        sites_avg(t).condition(cn).hs_tuned_band = struct();
        sites_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
        sites_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
        % initialize number of site pairs for each handspace
        % label
        for st = 1:size(lfp_tfa_cfg.analyse_states, 1)
            for hs = 1:size(lfp_tfa_cfg.conditions(1).hs_labels, 2)
                sites_avg(t).condition(cn).hs_tuned_band(st, hs).nsites = 0;
                sites_avg(t).condition(cn).hs_tuned_band(st, hs).freq.powspctrm = [];
            end
        end
        nsites = 0;
        for i = 1:length(lfp_tfr.session)
            for j = 1:length(lfp_tfr.session(i).sites)
                if ~strcmp(lfp_tfr.session(i).sites(j).target, lfp_tfa_cfg.compare.targets{t})
                    continue;
                end
                if ~lfp_tfr.session(i).sites(j).use_for_avg
                    continue;
                end
                if ~isempty(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band) && ...
                        isfield(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band, 'freq')
                    %nsites = nsites + 1;
                    for st = 1:size(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band, 1)
                        for hs = 1:size(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band, 2)
                            if isfield(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq, 'powspctrm') ...
                                    && ~isempty(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.powspctrm)
                                sites_avg(t).condition(cn).hs_tuned_band(st,hs).nsites = ...
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).nsites + 1;
                                if sites_avg(t).condition(cn).hs_tuned_band(st,hs).nsites == 1%~isfield(sites_avg.cond_based_tfs(cn).tfs_across_sessions, 'powspctrm')
                                    
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.time ...
                                        = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.time;
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.cfg ...
                                        = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.cfg;
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.freq ...
                                        = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.freq;
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.powspctrm ...
                                        = nanmean(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.powspctrm, 1);
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.gamma ...
                                        = nanmean(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.gamma, 1);
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.beta ...
                                        = nanmean(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.beta, 1);
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.alpha ...
                                        = nanmean(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.alpha, 1);
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.theta ...
                                        = nanmean(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.theta, 1);
                                    
                                    sites_avg(t).condition(cn).label = ...
                                        lfp_tfr.session(i).sites(j).condition(cn).label;
                                    sites_avg(t).condition(cn).cfg_condition = ...
                                        lfp_tfr.session(i).sites(j).condition(cn).cfg_condition;
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).hs_label ...
                                        = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).hs_label;
                                    if isfield(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs), 'state') && ...
                                            isfield(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs), 'state_name')
                                        sites_avg(t).condition(cn).hs_tuned_band(st,hs).state ...
                                            = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).state;
                                        sites_avg(t).condition(cn).hs_tuned_band(st,hs).state_name ...
                                            = lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).state_name;
                                    end
                                    
                                else
                                    ntimebins = size(sites_avg(t).condition(cn).hs_tuned_band(st, hs).freq.powspctrm, 3);
                                    % average same number of time bins
                                    if ntimebins > length(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.time)
                                        ntimebins = length(lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.time);
                                    end
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.powspctrm ...
                                        = cat(1, sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.powspctrm(:,:,1:ntimebins), ...
                                        nanmean((lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.powspctrm(:,:,1:ntimebins)), 1));
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.time = ...
                                        sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.time(1:ntimebins);
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.gamma ...
                                        = cat(1, sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.gamma(:,:,1:ntimebins), ...
                                        nanmean((lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.gamma(:,:,1:ntimebins)), 1));
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.beta ...
                                        = cat(1, sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.beta(:,:,1:ntimebins), ...
                                        nanmean((lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.beta(:,:,1:ntimebins)), 1));
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.alpha ...
                                        = cat(1, sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.alpha(:,:,1:ntimebins), ...
                                        nanmean((lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.alpha(:,:,1:ntimebins)), 1));
                                    sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.theta ...
                                        = cat(1, sites_avg(t).condition(cn).hs_tuned_band(st,hs).freq.theta(:,:,1:ntimebins), ...
                                        nanmean((lfp_tfr.session(i).sites(j).condition(cn).hs_tuned_band(st, hs).freq.theta(:,:,1:ntimebins)), 1));
                                end
                            end
                        end
                    end
                end
            end
        end
        
        
        
        if ~isempty(sites_avg(t).condition(cn).hs_tuned_band)
            if isfield(sites_avg(t).condition(cn).hs_tuned_band,...
                    'freq')
                plottitle = ['Target ' lfp_tfa_cfg.compare.targets{t}, ...
                    ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') ', ...
                    sites_avg(t).condition(cn).label];
                result_file = fullfile(results_fldr, ...
                    ['LFP_band_' lfp_tfa_cfg.compare.targets{t}, '_', ...
                    sites_avg(t).condition(cn).label]);
                lfp_tfa_plot_hs_tuned_tfr_band_average(sites_avg(t).condition(cn).hs_tuned_band, ...
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
            ['LFP_band_' lfp_tfa_cfg.compare.targets{t}, '_combined', ...
            trial_title]);
        lfp_tfa_plot_hs_tuned_tfr_band_average_combined(sites_avg(t).condition, ...
            lfp_tfa_cfg, plottitle, result_file, trial_type);
    end
end




% save session average tfs
save(fullfile(results_fldr, 'LFP_TFR_sites_avg.mat'), 'sites_avg');
close all

end