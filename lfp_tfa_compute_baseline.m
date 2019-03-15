function [ states_lfp, baseline ] = lfp_tfa_compute_baseline( states_lfp, cfg_tfs )

% computeBaseline - Computes mean and stddev of baseline power according to
% the given configuration
%
% USAGE:
%	[ states_lfp, baseline ] = lfp_tfa_compute_baseline(states_lfp, cfg_tfs)
%
% INPUTS:
%		states_lfp      - struct containing lfp power spectrogram per trial
%       cfg_tfs         - configuration structure for baseline computation
%           baseline_refstate   : reference state for baseline
%           baseline_period     : period around the reference state to be
%           considered for baseline computation
%           baseline_block      : block to be considered for baseline
%           calculation
%           choice_trial        : 1/0 whether to consider choice trials (1)
%           or instructed trials (0) for baseline computation
% OUTPUTS:
%		states_lfp      - input struct with added fields for baseline mean
%		and stddev for each site
%       baseline        - struct which saves the baseline mean and stddev for each site and 
%       configuration used for baseline calculation
% See also lfp_tfa_reject_noisy_trials

    % struct for storing results
    baseline = struct();
    baseline.cfg = cfg_tfs;
    fprintf('=============================================================\n');
    fprintf('Computing baseline ...\n');
    
    % loop through each site
    for i = 1:length( states_lfp )
        fprintf('Site %s \n', states_lfp(i).site_ID);
        % struct to store lfp power during baseline period for each trial
        baseline_pow = cell(1, length(states_lfp(i).trials));
        % loop through trial
        for t = 1:length( states_lfp(i).trials )
            trial = states_lfp(i).trials(t);
            % whether this trial should be considered for baseline calculation
            consider_trial = ~(trial.noisy) & (trial.block == cfg_tfs.baseline_block) ...
                & (trial.choice_trial == cfg_tfs.choice_trial);
            if consider_trial 
                if strcmp(cfg_tfs.baseline_ref_state, '') && strcmp(cfg_tfs.baseline_period , 'trial')
                    baseline_pow{t} = trial.tfs.powspctrm(1, :, ...
                        trial.tfs.time >= trial.trialperiod(1) & trial.tfs.time <= trial.trialperiod(2));
                elseif ~strcmp(cfg_tfs.baseline_ref_state, '')
                    ref = cfg_tfs.baseline_ref_state;
                    ref_onset = trial.states([trial.states.id] == ref) .onset_t;
                    ref_period = ref_onset + cfg_tfs.baseline_period;
                    baseline_pow{t} = trial.tfs.powspctrm(1, :, trial.tfs.time >= ...
                        ref_period(1) & trial.tfs.time <= ref_period(2));
                end
            end
        end
        % calculate baseline power mean and std
        arr_baseline_pow = cat(3, baseline_pow{:});
        baseline_pow_mean = nanmean(arr_baseline_pow, 3);
        baseline_pow_std = nanstd(arr_baseline_pow, 0, 3);
        states_lfp(i).baseline_mean = baseline_pow_mean;
        states_lfp(i).baseline_std = baseline_pow_std;
        
        baseline.sites(i).mean = baseline_pow_mean;
        baseline.sites(i).std = baseline_pow_std;
    end
    fprintf(' done\n');
    fprintf('=============================================================\n');
end

