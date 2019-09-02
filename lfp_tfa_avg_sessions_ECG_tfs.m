function sessions_avg = lfp_tfa_avg_sessions_ECG_tfs(evoked_ecg, lfp_tfa_cfg)
%lfp_tfa_avg_evoked_LFP_across_sessions  - Condition-based evoked LFP response
% average across many session averages
%
% USAGE:
%	sessions_avg = lfp_tfa_avg_sessions_ECG_tfs(evoked_ecg, lfp_tfa_cfg)
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
        sessions_avg(t).condition(cn).hs_tuned_evoked = struct();
        sessions_avg(t).condition(cn).cfg_condition = lfp_tfa_cfg.conditions(cn);
        sessions_avg(t).condition(cn).label = lfp_tfa_cfg.conditions(cn).label;
        
        % initialize number of site pairs for each handspace
        % label
        for st = 1:size(lfp_tfa_cfg.analyse_states, 1)
            for hs = 1:size(lfp_tfa_cfg.conditions(1).hs_labels, 2)
                sessions_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsessions = 0;
                sessions_avg(t).condition(cn).hs_tuned_tfs(st, hs).powspctrm = {};                
            end
        end  

        for i = 1:length(evoked_ecg.session)
            if ~isempty(evoked_ecg.session(i).condition(cn).hs_tuned_tfs) && ... 
                isfield(evoked_ecg.session(i).condition(cn).hs_tuned_tfs, 'powspctrm')
                for st = 1:size(evoked_ecg.session(i).condition(cn).hs_tuned_tfs, 1)
                    for hs = 1:size(evoked_ecg.session(i).condition(cn).hs_tuned_tfs, 2)
                        if isfield(evoked_ecg.session(i).condition(cn).hs_tuned_tfs(st, hs), 'powspctrm') ...
                                && ~isempty(evoked_ecg.session(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm)
                            sessions_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsessions = ...
                                sessions_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsessions + 1;
                            if sessions_avg(t).condition(cn).hs_tuned_tfs(st, hs).nsessions == 1
                                sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).time ...
                                    = evoked_ecg.session(i).condition(cn).hs_tuned_tfs(st, hs).time;
                                sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).hs_label ...
                                    = evoked_ecg.session(i).condition(cn).hs_tuned_tfs(st, hs).hs_label;
                                if isfield(evoked_ecg.session(i).session_avg(k).condition(cn).hs_tuned_tfs(st, hs), 'state') && ...
                                        isfield(evoked_ecg.session(i).condition(cn).hs_tuned_tfs(st, hs), 'state_name')
                                    sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).state ...
                                        = evoked_ecg.session(i).condition(cn).hs_tuned_tfs(st, hs).state;
                                    sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).state_name ...
                                        = evoked_ecg.session(i).condition(cn).hs_tuned_tfs(st, hs).state_name;
                                end
                            end
                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm ...
                                = [sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm; ...
                                evoked_ecg.session(i).condition(cn).hs_tuned_tfs(st, hs).powspctrm];                                                             
                        end
                    end
                end
            end                               
        end

        % compute average
        if isfield(sessions_avg(t).condition(cn).hs_tuned_tfs, 'powspctrm')
            for st = 1:size(sessions_avg(t).condition(cn).hs_tuned_tfs, 1)
                for hs = 1:size(sessions_avg(t).condition(cn).hs_tuned_tfs, 2)
                    if ~isempty(sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm)
                        sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm = ...
                            cat(1, sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm{:});
                        sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm_all = ...
                            sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm;        
                        sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).std = ...
                            std(sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm, 0, 1);
                        sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm = ...
                            mean(sessions_avg(t).condition(cn).hs_tuned_tfs(st,hs).powspctrm, 1);                           
                    end
                end
            end
        end


        if ~isempty(sessions_avg(t).condition(cn).hs_tuned_tfs)
            if isfield(sessions_avg(t).condition(cn).hs_tuned_tfs,... 
                    'mean')
                plottitle = [lfp_tfa_cfg.compare.targets{t},...
                     lfp_tfa_cfg.conditions(cn).label];
                result_file = fullfile(results_fldr, ...
                                ['ECG_Evoked_' lfp_tfa_cfg.conditions(cn).label '.png']);
                lfp_tfa_plot_hs_tuned_tfr_multiple_img(sessions_avg(t).condition(cn).hs_tuned_tfs, ...
                            lfp_tfa_cfg, plottitle, result_file);
            end
        end
        % save session average tfs
        save(fullfile(results_fldr, 'sessions_ECG_tfs.mat'), 'sessions_avg');
    end

    % difference between conditions
    sessions_avg(t).difference = [];
    for diff = 1:size(lfp_tfa_cfg.diff_condition, 2)
        diff_condition = lfp_tfa_cfg.diff_condition{diff};
        sessions_avg(t).difference = [sessions_avg(t).difference, ...
            lfp_tfa_compute_diff_condition_evoked(sessions_avg(t).condition, diff_condition)];
    end
    % plot Difference TFR
    for dcn = 1:length(sessions_avg(t).difference)
        if ~isempty(sessions_avg(t).difference(dcn).hs_tuned_tfs)
            if isfield(sessions_avg(t).difference(dcn).hs_tuned_tfs,... 
                    'mean')
                plottitle = ['Target ', lfp_tfa_cfg.compare.targets{t}, ...
                    sessions_avg(t).difference(dcn).label];
                result_file = fullfile(results_fldr, ...
                    ['ECG_DiffEvoked_' 'diff_condition' num2str(dcn) '.png']);
                    %sessions_avg(t).difference(dcn).label '.png']);
                lfp_tfa_plot_hs_tuned_tfr_multiple_img(sessions_avg(t).difference(dcn).hs_tuned_tfs, ...
                    lfp_tfa_cfg, plottitle, result_file);
            end
        end
    end
        
    close all;
end