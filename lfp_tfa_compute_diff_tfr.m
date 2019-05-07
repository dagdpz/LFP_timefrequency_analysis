function [ diff_tfr ] = lfp_tfa_compute_diff_tfr( lfp_tfr, lfp_tfa_cfg )
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
    for dcn = 1:2:length(lfp_tfr.condition)-1
        % initially load the pre-injection data structure
        preinj_tfr = lfp_tfr.condition(dcn);
        diff_tfr.difference(dcn) = preinj_tfr;
        if isempty(preinj_tfr.hs_tuned_tfs) 
            continue;
        end
        % post injection
        postinj_tfr = lfp_tfr.condition(dcn + 1);
        if isempty(postinj_tfr.hs_tuned_tfs)
            diff_tfr.difference(dcn) = postinj_tfr;
            continue;
        end
        
        if ~isfield(postinj_tfr.hs_tuned_tfs, 'powspctrm') || ~isfield(preinj_tfr.hs_tuned_tfs, 'powspctrm')
            continue;
        end
        % change the condition label
        diff_tfr.difference(dcn).label = strrep(diff_tfr.difference(dcn).label, ...
            'Pre', 'Post - Pre');        

        % loop through handspace tunings
        for hs = 1:size(postinj_tfr.hs_tuned_tfs, 2)
            for st = 1:size(postinj_tfr.hs_tuned_tfs, 1)
                if isfield(preinj_tfr.hs_tuned_tfs(st, hs), 'powspctrm') && ...
                        isfield(postinj_tfr.hs_tuned_tfs(st, hs), 'powspctrm') && ...
                        ~isempty(preinj_tfr.hs_tuned_tfs(st, hs).powspctrm) && ...
                    ~isempty(postinj_tfr.hs_tuned_tfs(st, hs).powspctrm)
                    ntimebins = min(length(postinj_tfr.hs_tuned_tfs(st, hs).time), ...
                        length(diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).time));
                    % calculate the difference
                    diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).powspctrm = ...
                        postinj_tfr.hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins) - ...
                        preinj_tfr.hs_tuned_tfs(st, hs).powspctrm(1,:,1:ntimebins);
                    diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).time = ...
                        diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).time(1:ntimebins);
                    if isfield(diff_tfr.difference(dcn).hs_tuned_tfs(st, hs), 'ntrials')
                        diff_tfr.difference(dcn).hs_tuned_tfs(st, hs).ntrials = [];
                    end
                end
            end
        end
    end
    % generate return variable
    diff_tfr = diff_tfr.difference;
end

