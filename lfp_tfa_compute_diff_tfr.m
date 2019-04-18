function [ diff_tfr ] = lfp_tfa_compute_diff_tfr( lfp_tfr, lfp_tfa_cfg )
%lfp_tfa_compute_diff_tfr - function to compute the difference in time freq
%response between control and inactivation trials
%   Detailed explanation goes here
    diff_tfr = struct();
    %npostinj_conditions = npreinj_condtions;
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

