function [ site_lfp ] = lfp_tfa_compute_site_baseline( site_lfp, lfp_tfa_cfg )

% lfp_tfa_compute_baseline_power - Computes mean and stddev of baseline power according to
% the given configuration
%
% USAGE:
%	session_proc_lfp_out = lfp_tfa_compute_baseline_power( session_lfp_in, cfg_tfs )
%
% INPUTS:
%		sites_lfp_folder        - struct containing LFP data strcture for all
%		sites for one session, output of lfp_tfa_reject_noisy_lfp or
%		lfp_tfa_process_lfp
%       cfg_tfs                 - configuration structure for baseline computation
%           Required Fields: see lfp_tfa_settings
%           baseline_refstate           - reference state for baseline
%           baseline_period             - period around the reference state to be
%           considered for baseline computation
%           baseline_perturbation       - perturbation group to be considered for baseline
%           calculation
%           baseline_use_choice_trial   - 1/0 whether to consider choice trials (1)
%           or instructed trials (0) for baseline computation
% OUTPUTS:
%		session_proc_lfp_out      - same as input struct with additional
%		fields
%           baseline_mean   - mean spectral power in baseline period
%           baseline_std    - standard dev LFP spectral power in baseline
%           perios
%
% See also lfp_tfa_reject_noisy_lfp, lfp_tfa_process_lfp, lfp_tfa_settings,
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

    %fprintf('=============================================================\n');
    fprintf('Computing baseline ... \n');
    
    % loop through each site
    %for i = 1:length( session_lfp_in )
    % structure to store baseline
    baseline = struct();        
    % baseline conditions
    perturbation = unique([site_lfp.trials.perturbation]);
    choice = unique([site_lfp.trials.choice_trial]);
    % loop through each baseline condition
    cnd = 1;
    for p = perturbation
        for c = choice
            baseline(cnd).perturbation = p;
            baseline(cnd).choice = c;
            % struct to store lfp power during baseline period for each trial
            baseline_pow = cell(1, length(site_lfp.trials));
            % loop through trial
            for t = 1:length( site_lfp.trials )
                trial = site_lfp.trials(t);
                % whether this trial should be considered for baseline calculation
                consider_trial = ~(trial.noisy) & sum(trial.perturbation == p) ...
                    & sum(trial.choice_trial == c);
                if consider_trial 
                    if strcmp(lfp_tfa_cfg.baseline_ref_state, '') && strcmp(lfp_tfa_cfg.baseline_ref_period , 'trial')
                        baseline_pow{t} = trial.tfs.powspctrm(1, :, ...
                            trial.tfs.time >= trial.trialperiod(1) & trial.tfs.time <= trial.trialperiod(2));
                    elseif ~strcmp(lfp_tfa_cfg.baseline_ref_state, '')
                        ref = lfp_tfa_cfg.baseline_ref_state;
                        ref_onset = trial.states([trial.states.id] == ref) .onset_t;
                        ref_period = ref_onset + lfp_tfa_cfg.baseline_ref_period;
                        baseline_pow{t} = trial.tfs.powspctrm(1, :, trial.tfs.time >= ...
                            ref_period(1) & trial.tfs.time <= ref_period(2));
                    end
                end
            end
            % calculate baseline power mean and std
            arr_baseline_pow = cat(3, baseline_pow{:});
            baseline(cnd).pow_mean = nanmean(arr_baseline_pow, 3);
            baseline(cnd).pow_std = nanstd(arr_baseline_pow, 0, 3);
            site_lfp.baseline = baseline;
            cnd = cnd + 1;
            clear('arr_baseline_pow'); 
        end
    end

    clear('baseline_pow');
        
    %end    
    
    fprintf(' done\n');
    %fprintf('=============================================================\n');
end


