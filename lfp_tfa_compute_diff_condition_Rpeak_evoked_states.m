function [ diff_Rpeak_evoked_states ] = lfp_tfa_compute_diff_condition_Rpeak_evoked_states( Rpeak_evoked_states, diff_condition, diff_color )
%lfp_tfa_compute_diff_tfr - function to compute the difference in time freq
%response between control and inactivation trials
%
% USAGE:
%	diff_tfr = lfp_tfa_compute_diff_tfr( lfp_tfr, lfp_tfa_cfg )
%
% INPUTS:
%       lfp_tfr         - struct containing the condition-based average LFP time freq spectrogram
%       for individual sites or average across sites or average across
%       sessions, see lfp_tfa_site_average_tfr,
%       lfp_tfa_avg_tfr_across_sessions, lfp_tfa_avg_across_sites
%		lfp_tfa_cfg     - struct containing the required settings
%
% OUTPUTS:
%		diff_tfr        - struct containing the condition-based LFP time freq spectrogram
%       average difference between post and pre injection
%
% REQUIRES:	
%
% See also lfp_tfa_site_average_tfr, lfp_tfa_avg_tfr_across_sessions, lfp_tfa_avg_across_sites 
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

    conditions = [Rpeak_evoked_states.cfg_condition];
    
    for i = 1:length(diff_condition)/2
        diff_Rpeak_evoked_states = struct();
        compare.field = diff_condition{1, 2*(i-1) + 1};
        compare.values = diff_condition{1, 2*(i-1) + 2};
        % check if conditions to compare exist
        if strcmp(compare.field, 'perturbation')
            if sum([compare.values{:}] == unique([conditions.perturbation])) < 2
                continue;
            end
        elseif strcmp(compare.field, 'choice')
            if sum([compare.values{:}] == unique([conditions.choice])) < 2
                continue;
            end
        elseif strcmp(compare.field, 'type_eff')
            if sum(ismember(vertcat(compare.values{:}), ...
                unique([conditions.type; conditions.effector]', 'rows'), 'rows')) < 2
                continue;
            end
        end
    
        dcn = 0;
        traversed_idx = [];
        for cn = 1:length(Rpeak_evoked_states)
            condition_found = false;
            if strcmp(compare.field, 'choice')
                condition_found = Rpeak_evoked_states(cn).cfg_condition.choice == compare.values{1};

            elseif strcmp(compare.field, 'perturbation')
                condition_found = Rpeak_evoked_states(cn).cfg_condition.perturbation == compare.values{1};
                
            elseif strcmp(compare.field, 'type_eff')
                condition_found = Rpeak_evoked_states(cn).cfg_condition.type == compare.values{1}(1) ...
                    & Rpeak_evoked_states(cn).cfg_condition.effector == compare.values{1}(2);

            end
            % initially load the pre-injection data structure
            if condition_found
                traversed_idx = [traversed_idx cn];            
            else
                continue;
            end
            for d = 1:length(Rpeak_evoked_states)
                if any(traversed_idx == d), continue; end
                comparison_pair_found = false;

                if strcmp(compare.field, 'choice')
                    comparison_pair_found = Rpeak_evoked_states(d).cfg_condition.type == Rpeak_evoked_states(cn).cfg_condition.type ...
                        & Rpeak_evoked_states(d).cfg_condition.effector == Rpeak_evoked_states(cn).cfg_condition.effector ...
                        & Rpeak_evoked_states(d).cfg_condition.choice == compare.values{2} ...
                        & Rpeak_evoked_states(d).cfg_condition.perturbation == Rpeak_evoked_states(cn).cfg_condition.perturbation;

                elseif strcmp(compare.field, 'perturbation')
                    comparison_pair_found = Rpeak_evoked_states(d).cfg_condition.type == Rpeak_evoked_states(cn).cfg_condition.type ...
                        & Rpeak_evoked_states(d).cfg_condition.effector == Rpeak_evoked_states(cn).cfg_condition.effector ...
                        & Rpeak_evoked_states(d).cfg_condition.choice == Rpeak_evoked_states(cn).cfg_condition.choice ...
                        & Rpeak_evoked_states(d).cfg_condition.perturbation == compare.values{2};
                    
                elseif strcmp(compare.field, 'type_eff')
                    comparison_pair_found = Rpeak_evoked_states(d).cfg_condition.type == compare.values{2}(1) ...
                        & Rpeak_evoked_states(d).cfg_condition.effector == compare.values{2}(2) ...
                        & Rpeak_evoked_states(d).cfg_condition.choice == Rpeak_evoked_states(cn).cfg_condition.choice ...
                        & Rpeak_evoked_states(d).cfg_condition.perturbation == Rpeak_evoked_states(cn).cfg_condition.perturbation;
                end
                if comparison_pair_found

                    dcn = dcn + 1;
                    %diff_tfr.difference(dcn) = struct();
                    % pre injection
                    cond1_evoked = Rpeak_evoked_states(cn);
                    if isempty(cond1_evoked.Rpeak_evoked) 
                        continue;
                    end
                    % post injection
                    cond2_evoked = Rpeak_evoked_states(d);
                    if isempty(cond2_evoked.Rpeak_evoked)
                        continue;
                    end

%                     if ~isfield(cond2_evoked.Rpeak_evoked, 'mean')
%                         continue;
%                     end
                    
                    diff_Rpeak_evoked_states.difference(dcn) = cond2_evoked;
                    
                    % change the condition label
                    diff_Rpeak_evoked_states.difference(dcn).label = ['( ' cond2_evoked.label, ' vs. ', ...
                                cond1_evoked.label ' )'];
                            
                    diff_Rpeak_evoked_states.difference(dcn).cfg_condition = cond2_evoked.cfg_condition;
                    legend = cell(1, 2);
                    if strcmp(compare.field, 'choice')                        
                        diff_Rpeak_evoked_states.difference(dcn).cfg_condition.choice = ['diff' num2str(i)];
                        if cond2_evoked.cfg_condition.choice == 0
                            legend(1) = {'Instructed'};
                        else
                            legend(1) = {'Choice'};
                        end
                        if cond1_evoked.cfg_condition.choice == 0
                            legend(2) = {'Instructed'};
                        else
                            legend(2) = {'Choice'};
                        end
                    elseif strcmp(compare.field, 'perturbation')
                        diff_Rpeak_evoked_states.difference(dcn).cfg_condition.perturbation = ['diff' num2str(i)];
                        if cond2_evoked.cfg_condition.perturbation == 0
                            legend(1) = {'Pre-injection'};
                        else
                            legend(1) = {'Post-injection'};
                        end
                        if cond1_evoked.cfg_condition.perturbation == 0
                            legend(2) = {'Pre-injection'};
                        else
                            legend(2) = {'Post-injection'};
                        end
                    elseif strcmp(compare.field, 'type_eff')
                        diff_Rpeak_evoked_states.difference(dcn).cfg_condition.type_eff = ['diff' num2str(i)];
                    end

                    % loop through handspace tunings
                    diff_Rpeak_evoked_states.difference(dcn).Rpeak_evoked = cond2_evoked.Rpeak_evoked;
                    for hs = 1:size(cond2_evoked.Rpeak_evoked, 2)
                        for st = 1:size(cond2_evoked.Rpeak_evoked, 1)
                            
                            if isfield(cond1_evoked.Rpeak_evoked(st, hs), 'abs_timeprob') && ...
                                    isfield(cond2_evoked.Rpeak_evoked(st, hs), 'abs_timeprob') && ...
                                    ~isempty(fieldnames(cond1_evoked.Rpeak_evoked(st, hs).abs_timeprob)) && ...
                                    ~isempty(fieldnames(cond2_evoked.Rpeak_evoked(st, hs).abs_timeprob)) && ...
                                    isfield(cond1_evoked.Rpeak_evoked(st, hs), 'rel_timeprob') && ...
                                    isfield(cond2_evoked.Rpeak_evoked(st, hs), 'rel_timeprob') && ...
                                    ~isempty(fieldnames(cond1_evoked.Rpeak_evoked(st, hs).rel_timeprob)) && ...
                                    ~isempty(fieldnames(cond2_evoked.Rpeak_evoked(st, hs).rel_timeprob))
                                diff_Rpeak_evoked_states.difference(dcn).Rpeak_evoked(st, hs).abs_timefromRpeak = ...
                                    {cond2_evoked.Rpeak_evoked(st, hs).abs_timefromRpeak, ...
                                    cond1_evoked.Rpeak_evoked(st, hs).abs_timefromRpeak};
                                diff_Rpeak_evoked_states.difference(dcn).Rpeak_evoked(st, hs).abs_timeprob = ...
                                    [cond2_evoked.Rpeak_evoked(st, hs).abs_timeprob, ...
                                    cond1_evoked.Rpeak_evoked(st, hs).abs_timeprob];
                                diff_Rpeak_evoked_states.difference(dcn).Rpeak_evoked(st, hs).rel_timefromRpeak = ...
                                    {cond2_evoked.Rpeak_evoked(st, hs).rel_timefromRpeak, ...
                                    cond1_evoked.Rpeak_evoked(st, hs).rel_timefromRpeak};
                                diff_Rpeak_evoked_states.difference(dcn).Rpeak_evoked(st, hs).rel_timeprob = ...
                                    [cond2_evoked.Rpeak_evoked(st, hs).rel_timeprob, ...
                                    cond1_evoked.Rpeak_evoked(st, hs).rel_timeprob];
                                if isfield(cond2_evoked.Rpeak_evoked(st, hs), 'ntrials')
                                    diff_Rpeak_evoked_states.difference(dcn).Rpeak_evoked(st, hs).ntrials = ...
                                        [];
                                end
                                diff_Rpeak_evoked_states.difference(dcn).Rpeak_evoked(st, hs).legend = ...
                                    legend;
                                if nargin > 2
                                    diff_Rpeak_evoked_states.difference(dcn).Rpeak_evoked(st, hs).colors = ...
                                        diff_color;
                                end
                                
                            else
                                %diff_evoked.difference(dcn).Rpeak_evoked(st, hs).lfp = [];
                                diff_Rpeak_evoked_states.difference(dcn).Rpeak_evoked(st, hs).abs_timefromRpeak = [];
                                diff_Rpeak_evoked_states.difference(dcn).Rpeak_evoked(st, hs).abs_timeprob = [];
                                diff_Rpeak_evoked_states.difference(dcn).Rpeak_evoked(st, hs).rel_timefromRpeak = [];
                                diff_Rpeak_evoked_states.difference(dcn).Rpeak_evoked(st, hs).rel_timeprob = [];
                                diff_Rpeak_evoked_states.difference(dcn).Rpeak_evoked(st, hs).ntrials = [];
                            end
                        end
                    end
                else
                    continue;
                end
            end
        end
        if isfield(diff_Rpeak_evoked_states, 'difference')
            Rpeak_evoked_states = diff_Rpeak_evoked_states.difference;
        end
    end
    % generate return variable
    if isfield(diff_Rpeak_evoked_states, 'difference')
        diff_Rpeak_evoked_states = diff_Rpeak_evoked_states.difference;
    else
        diff_Rpeak_evoked_states = [];
    end
end

