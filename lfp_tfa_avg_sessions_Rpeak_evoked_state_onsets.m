function sessions_avg = lfp_tfa_avg_sessions_Rpeak_evoked_state_onsets(Rpeak_state_onset, lfp_tfa_cfg)
%lfp_tfa_avg_evoked_LFP_across_sessions  - Condition-based evoked LFP response
% average across many session averages
%
% USAGE:
%	sessions_avg = lfp_tfa_avg_sessions_ECG_evoked(evoked_ecg, lfp_tfa_cfg)
%
% INPUTS:
%		lfp_evoked		- struct containing the condition-based evoked LFP response for
%		indiviual sites, output of lfp_tfa_plot_site_evoked_LFP.m
%           Required Fields:
%               1. session.session_avg - 1xN struct containing condition-based
%               average evoked LFP response for N sessions (session_avg =
%               Average of site averages for one session)
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields:
%               1. conditions          - trial conditions to compare, see
%               lfp_tfa_settings.m and lfp_tfa_compare_conditions.m
%               2. root_results_fldr   - root folder where results are saved
%               3. compare.targets     - targets to compare, see lfp_tfa_settings.m
%               4. ref_hemisphere      - reference hemisphere for ipsi and
%               contra labeling
% OUTPUTS:
%		sessions_avg    - structure containing condition-based evoked LFP
%		response averaged across multiple sessions
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
    results_fldr = fullfile(lfp_tfa_cfg.root_results_fldr, 'ECG analysis');
    if ~exist(results_fldr, 'dir')
        mkdir(results_fldr);
    end
    
    % Average Evoked LFP response across sessions
    sessions_avg = struct();
    t = 1;
    
    for cn = 1:length(lfp_tfa_cfg.conditions)
        fprintf('Condition %s\n', lfp_tfa_cfg.conditions(cn).label);
        sessions_avg(t).condition(cn).Rpeak_evoked = struct();
        sessions_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
        sessions_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
        
        % initialize number of site pairs for each handspace
        % label
        for st = 1:size(Rpeak_state_onset.session(end).condition(cn).Rpeak_evoked, 1)
            for hs = 1:size(Rpeak_state_onset.session(end).condition(cn).Rpeak_evoked, 2)
                sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).nsessions = 0;
                sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).abs_timefromRpeak = {}; 
                sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).rel_timefromRpeak = {}; 
                sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).abs_timeprob.timebins = {}; 
                sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).abs_timeprob.prob = {}; 
                sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).rel_timeprob.prob = []; 
            end
        end  

        for i = 1:length(Rpeak_state_onset.session)
            if isempty(Rpeak_state_onset.session(i).condition)
                continue;
            end
            if ~isempty(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked) && ... 
                isfield(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked, 'abs_timefromRpeak') && ...
                isfield(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked, 'rel_timefromRpeak')
                for st = 1:size(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked, 1)
                    for hs = 1:size(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked, 2)
                        if isfield(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs), 'abs_timefromRpeak') ...
                                && ~isempty(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).abs_timefromRpeak) ...
                                && isfield(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs), 'rel_timefromRpeak') ...
                                && ~isempty(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).rel_timefromRpeak)
                            sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).nsessions = ...
                                sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).nsessions + 1;
                            if sessions_avg(t).condition(cn).Rpeak_evoked(st, hs).nsessions == 1
                                sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).hs_label ...
                                    = Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).hs_label;
                                if isfield(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs), 'state') && ...
                                        isfield(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs), 'state_name')
                                    sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).state ...
                                        = Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).state;
                                    sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).state_name ...
                                        = Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).state_name;
                                end
                            end
                            sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timefromRpeak ...
                                = [sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timefromRpeak, ...
                                Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).abs_timefromRpeak];
                            sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).rel_timefromRpeak ...
                                = [sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).rel_timefromRpeak, ...
                                Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).rel_timefromRpeak];
                            % accumulate histogram counts
                            sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.timebins ...
                                = [sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.timebins, ...
                                Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).abs_timeprob.timebins];
                            sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.prob ...
                                = [sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.prob, ...
                                Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).abs_timeprob.prob];
                            sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).rel_timeprob.timebins ...
                                = Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).rel_timeprob.timebins;
                            sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).rel_timeprob.prob ...
                                = [sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).rel_timeprob.prob; ...
                                Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs).rel_timeprob.prob];
                        end
                    end
                end
            end                               
        end
        
        % concatenate probability histograms
        for st = 1:size(sessions_avg(t).condition(cn).Rpeak_evoked, 1)
            for hs = 1:size(sessions_avg(t).condition(cn).Rpeak_evoked, 2)
                if isfield(Rpeak_state_onset.session(i).condition(cn).Rpeak_evoked(st, hs), 'abs_timeprob')
                    nsessions = length(sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.timebins);
                    [max_ntimebins, idx] = max(cellfun('size', ...
                        sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.timebins, 2));
                    sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.timebins = ...
                        sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.timebins{idx};
                    state_probability = nan*ones(nsessions, max_ntimebins-1);
                    for i = 1:nsessions
                        nprob = length(sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.prob{i});
                        state_probability(i, 1:nprob) = ...
                            sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.prob{i};
                    end
                    sessions_avg(t).condition(cn).Rpeak_evoked(st,hs).abs_timeprob.prob = ...
                        state_probability;
                end
            end
        end                    

        if ~isempty(sessions_avg(t).condition(cn).Rpeak_evoked)
            if isfield(sessions_avg(t).condition(cn).Rpeak_evoked,... 
                    'abs_timeprob') && isfield(...
                    sessions_avg(t).condition(cn).Rpeak_evoked, 'rel_timeprob'); 
                plottitle = [lfp_tfa_cfg.compare.targets{t},...
                     lfp_tfa_cfg.conditions(cn).label];
                result_file = fullfile(results_fldr, ...
                                ['Rpeak_Evoked_state_onsets_' lfp_tfa_cfg.conditions(cn).label '.png']);
                lfp_tfa_plot_Rpeak_ref_state_onsets (sessions_avg(t).condition(cn).Rpeak_evoked, ...
                            lfp_tfa_cfg, plottitle, result_file);
            end
        end
        % save session average tfs
        save(fullfile(results_fldr, 'sessions_evoked_ECG.mat'), 'sessions_avg');
    end
 
        
    %close all;
end