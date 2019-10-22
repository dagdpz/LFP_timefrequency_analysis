function sessions_avg = lfp_tfa_avg_pow_across_sessions(lfp_pow, lfp_tfa_cfg)
%lfp_tfa_avg_pow_across_sessions  - Condition-based LFP power spectrum
%(Power vs. Frequency) grand average across many session averages (A
%session average is the average across the sites recorded in one session)
%
% USAGE:
%	sessions_avg = lfp_tfa_avg_pow_across_sessions(lfp_evoked, lfp_tfa_cfg)
%
% INPUTS:
%		lfp_pow		    - struct containing the condition-based average LFP
%       power spectrum for multiple sessions, i.e., the output of 
%       lfp_tfa_plot_site_powspctrum
%           Required Fields:
%               1. session.session_avg - 1xN struct containing condition-based
%               average LFP power spectrum for N sessions (session_avg =
%               Average of site averages for one session)
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields:
%               1. conditions          - trial conditions to compare, see
%               lfp_tfa_settings.m and lfp_tfa_compare_conditions.m
%               2. root_results_fldr   - root folder where results are 
%               saved. Results will be saved under 
%               [lfp_tfa_cfg.root_results_fldr ...
%               '/Avg_across_sessions/LFP_Power']
%               3. compare.targets     - target areas to compare, , 
%               see settings/lfp_tfa_settings_example
%               4. ref_hemisphere      - reference hemisphere for contra-
%               and ipsi- labelling, see settings/lfp_tfa_settings_example
% OUTPUTS:
%		sessions_avg    - structure containing condition-based LFP Power
%		spectrum response averaged across multiple sessions
%
% REQUIRES:	lfp_tfa_plot_hs_tuned_psd2
%
% See also lfp_tfa_settings, lfp_tfa_define_settings, lfp_tfa_compare_conditions, 
% lfp_tfa_plot_hs_tuned_psd
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
    results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'Avg_across_sessions', 'LFP_Power');
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % Average LFP Power spectrum across sessions
    sessions_avg = struct();

    for t = 1:length(lfp_tfa_cfg.compare.targets)
        sessions_avg(t).target = lfp_tfa_cfg.compare.targets{t};
        for cn = 1:length(lfp_tfa_cfg.conditions)
            fprintf('Condition %s\n', lfp_tfa_cfg.conditions(cn).label);
            sessions_avg(t).condition(cn).avg_across_sessions = struct();
            nsessions = 0;
            for i = 1:length(lfp_pow.session)
                for k = 1:length(lfp_pow.session(i).session_avg)
                    if strcmp(lfp_pow.session(i).session_avg(k).target, lfp_tfa_cfg.compare.targets{t})
                        if ~isempty(lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power) && ... 
                            isfield(lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power, 'mean')
                            nsessions = nsessions + 1;   
                            for ep = 1:size(lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power, 1)
                                for hs = 1:size(lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power, 2)
                                    if isfield(lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power(ep, hs), 'mean') ...
                                            && ~isempty(lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power(ep, hs).mean)
                                        if nsessions == 1%~isfield(sessions_avg(t).cond_based_tfs(cn).tfs_across_sessions, 'powspctrm')
                                            sessions_avg(t).condition(cn).avg_across_sessions(ep,hs).freq ...
                                            = lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power(ep, hs).freq;
                                            sessions_avg(t).condition(cn).avg_across_sessions(ep,hs).hs_label ...
                                                = lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power(ep, hs).hs_label;
                                            sessions_avg(t).condition(cn).avg_across_sessions(ep,hs).epoch_name ...
                                                = lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power(ep, hs).epoch_name;
                                            sessions_avg(t).condition(cn).avg_across_sessions(ep,hs).mean ...
                                                = lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power(ep, hs).mean;
                                            if isfield(lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power(ep, hs), 'nsites')
                                                sessions_avg(t).condition(cn).avg_across_sessions(ep,hs).nsites ...
                                                    = lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power(ep, hs).nsites;
                                            end
                                        else
                                            nfreqbins = length(sessions_avg(t).condition(cn).avg_across_sessions(ep, hs).freq);
                                            % average same number of time bins
                                            if nfreqbins > length(lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power(ep, hs).freq)
                                                nfreqbins = length(lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power(ep, hs).freq);
                                            end
                                            sessions_avg(t).condition(cn).avg_across_sessions(ep,hs).mean ...
                                                = (lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power(ep, hs).mean(1:nfreqbins)) + ...
                                                sessions_avg(t).condition(cn).avg_across_sessions(ep,hs).mean(1:nfreqbins);
                                            if isfield(sessions_avg(t).condition(cn).avg_across_sessions(ep,hs), 'nsites')
                                                sessions_avg(t).condition(cn).avg_across_sessions(ep,hs).nsites ...
                                                    = lfp_pow.session(i).session_avg(k).condition(cn).hs_tuned_power(ep, hs).nsites + ...
                                                    sessions_avg(t).condition(cn).avg_across_sessions(ep,hs).nsites;
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
            for ep = 1:size(sessions_avg(t).condition(cn).avg_across_sessions, 1)
                for hs = 1:size(sessions_avg(t).condition(cn).avg_across_sessions, 2)
                    if isfield(sessions_avg(t).condition(cn).avg_across_sessions(ep,hs), 'mean')
                        sessions_avg(t).condition(cn).avg_across_sessions(ep,hs).nsessions = nsessions;                                
                        sessions_avg(t).condition(cn).avg_across_sessions(ep,hs).mean = ...
                            (1/nsessions) * sessions_avg(t).condition(cn).avg_across_sessions(ep,hs).mean;
                    end
                end
            end

            % plot
            if ~isempty(sessions_avg(t).condition(cn).avg_across_sessions)
                if isfield(sessions_avg(t).condition(cn).avg_across_sessions,... 
                        'mean')
                    plottitle = [lfp_tfa_cfg.compare.targets{t}, ...
                        ' (ref_', lfp_tfa_cfg.ref_hemisphere, ') ', ...
                        lfp_tfa_cfg.conditions(cn).label];
                    result_file = fullfile(results_fldr, ...
                                    ['LFP_Power_' lfp_tfa_cfg.compare.targets{t} '_' lfp_tfa_cfg.conditions(cn).label '.png']);
                    lfp_tfa_plot_hs_tuned_psd_2(sessions_avg(t).condition(cn).avg_across_sessions, ...
                                lfp_tfa_cfg, plottitle, result_file);
                end
            end
            
        end
    end
    % save session average tfs
    save(fullfile(results_fldr, 'LFP_Power_sessions_avg.mat'), 'sessions_avg');
end