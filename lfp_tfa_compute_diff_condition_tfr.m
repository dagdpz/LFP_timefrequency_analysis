function [ diff_tfr ] = lfp_tfa_compute_diff_condition_tfr( lfp_tfr, varargin )
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

    diff_tfr = struct();
    
    % get the fields to compare and their values
    if nargin == 3
        compare.field = varargin{1};
        compare.values = varargin{2};
    else
        fprintf('Please enter a field and values to compare')
        return
    end
    
    dcn = 0;
    traversed_idx = [];
    for cn = 1:length(lfp_tfr.condition)-1
        condition_found = false;
        if strcmp(compare.field, 'choice')
            condition_found = lfp_tfr.condition(cn).cfg_condition.choice == compare.values{1};
                
        elseif strcmp(compare.field, 'perturbation')
            condition_found = lfp_tfr.condition(cn).cfg_condition.perturbation == compare.values{1};
            
        end
        % initially load the pre-injection data structure
        if condition_found
            traversed_idx = [traversed_idx cn];            
        else
            continue;
        end
        for d = 1:length(lfp_tfr.condition)
            if any(traversed_idx == d), continue; end
            comparison_pair_found = false;
            
            if strcmp(compare.field, 'choice')
                comparison_pair_found = lfp_tfr.condition(d).cfg_condition.type == lfp_tfr.condition(cn).cfg_condition.type ...
                    & lfp_tfr.condition(d).cfg_condition.effector == lfp_tfr.condition(cn).cfg_condition.effector ...
                    & lfp_tfr.condition(d).cfg_condition.choice == compare.values{2} ...
                    & lfp_tfr.condition(d).cfg_condition.perturbation == lfp_tfr.condition(cn).cfg_condition.perturbation;
            
            elseif strcmp(compare.field, 'perturbation')
                comparison_pair_found = lfp_tfr.condition(d).cfg_condition.type == lfp_tfr.condition(cn).cfg_condition.type ...
                    & lfp_tfr.condition(d).cfg_condition.effector == lfp_tfr.condition(cn).cfg_condition.effector ...
                    & lfp_tfr.condition(d).cfg_condition.choice == lfp_tfr.condition(cn).cfg_condition.choice ...
                    & lfp_tfr.condition(d).cfg_condition.perturbation == compare.values{2};
            end
            if comparison_pair_found
                
                dcn = dcn + 1;
                % pre injection
                preinj_sync = lfp_tfr.condition(cn);
                diff_tfr.difference(dcn) = preinj_sync;
                if isempty(preinj_sync.hs_tuned_tfs) 
                    continue;
                end
                % post injection
                postinj_tfr = lfp_tfr.condition(d);
                if isempty(postinj_tfr.hs_tuned_tfs)
                    diff_tfr.difference(dcn) = postinj_tfr;
                    continue;
                end
        
                if ~isfield(postinj_tfr.hs_tuned_tfs, 'powspctrm') || ~isfield(preinj_sync.hs_tuned_tfs, 'powspctrm')
                    continue;
                end
                
                % change the condition label
                if strcmp(compare.field, 'choice')
                    diff_tfr.difference(dcn).label = strrep(diff_tfr.difference(dcn).label, ...
                        'Instr', 'Diff_Choice');  
                    diff_tfr.difference(dcn).cfg_condition.choice = 'diff';
                elseif strcmp(compare.field, 'perturbation')
                    diff_tfr.difference(dcn).label = strrep(diff_tfr.difference(dcn).label, ...
                        'Pre', 'Diff_Perturb');  
                    diff_tfr.difference(dcn).cfg_condition.choice = 'diff';
                elseif strcmp(compare.field, 'type_eff')
                    diff_tfr.difference(dcn).label = strrep(diff_tfr.difference(dcn).label, ...
                        regexp(diff_tfr.difference(dcn).label, 'Type_\d_Eff_\d', 'match'), ...
                        'Diff_Task');  
                    diff_tfr.difference(dcn).cfg_condition.type_eff = 'diff';
   
                end

                % loop through handspace tunings
                for hs = 1:size(postinj_tfr.hs_tuned_tfs, 2)
                    for st = 1:size(postinj_tfr.hs_tuned_tfs, 1)
                        if isfield(preinj_sync.hs_tuned_tfs(st, hs), 'powspctrm') && ...
                                isfield(postinj_tfr.hs_tuned_tfs(st, hs), 'powspctrm') && ...
                                ~isempty(preinj_sync.hs_tuned_tfs(st, hs).powspctrm) && ...
                            ~isempty(postinj_tfr.hs_tuned_tfs(st, hs).powspctrm)
                            ntimebins = min([size(postinj_tfr.hs_tuned_tfs(st, hs).powspctrm, 3), ...
                                size(preinj_sync.hs_tuned_tfs(st, hs).powspctrm, 3), ...
                                size(diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).powspctrm, 3)]);
                            % calculate the difference
                            diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).powspctrm = ...
                                postinj_tfr.hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins) - ...
                                preinj_sync.hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins);
                            diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).time = ...
                                diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).time(1:ntimebins);
                            if isfield(diff_tfr.difference(dcn).hs_tuned_tfs(st, hs), 'ntrials')
                                diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).ntrials = [];
                            end
                        end
                    end
                end
            else
                continue;
            end
        end
    end
    % generate return variable
    if isfield(diff_tfr, 'difference')
        diff_tfr = diff_tfr.difference;
    end
end

