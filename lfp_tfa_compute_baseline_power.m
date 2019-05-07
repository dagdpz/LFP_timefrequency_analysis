function [ session_proc_lfp_out ] = lfp_tfa_compute_baseline_power( session_lfp_in, cfg_tfs )

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

    % struct for storing results
    session_proc_lfp_out = session_lfp_in;
    %baseline.cfg = cfg_tfs;
    fprintf('=============================================================\n');
    fprintf('Computing baseline ...\n');
    
    % folder for saving results
    root_results_folder = cfg_tfs.results_folder;
    results_folder = fullfile(root_results_folder);
    if ~exist(results_folder, 'dir')
        mkdir(results_folder);
    end
    
    % loop through each site
    for i = 1:length( session_lfp_in )
        site_lfp = session_lfp_in(i);
        fprintf('Processing site (%g/%g): %s \n', i, length(session_lfp_in), site_lfp.site_ID);
        % struct to store lfp power during baseline period for each trial
        baseline_pow = cell(1, length(site_lfp.trials));
        % loop through trial
        for t = 1:length( site_lfp.trials )
            trial = site_lfp.trials(t);
            % whether this trial should be considered for baseline calculation
            consider_trial = ~(trial.noisy) & sum(trial.perturbation == cfg_tfs.baseline_perturbation) ...
                & (trial.choice_trial == cfg_tfs.baseline_use_choice_trial);
            if consider_trial 
                if strcmp(cfg_tfs.baseline_ref_state, '') && strcmp(cfg_tfs.baseline_ref_period , 'trial')
                    baseline_pow{t} = trial.tfs.powspctrm(1, :, ...
                        trial.tfs.time >= trial.trialperiod(1) & trial.tfs.time <= trial.trialperiod(2));
                elseif ~strcmp(cfg_tfs.baseline_ref_state, '')
                    ref = cfg_tfs.baseline_ref_state;
                    ref_onset = trial.states([trial.states.id] == ref) .onset_t;
                    ref_period = ref_onset + cfg_tfs.baseline_ref_period;
                    baseline_pow{t} = trial.tfs.powspctrm(1, :, trial.tfs.time >= ...
                        ref_period(1) & trial.tfs.time <= ref_period(2));
                end
            end
        end
        % calculate baseline power mean and std
        arr_baseline_pow = cat(3, baseline_pow{:});
        baseline_pow_mean = nanmean(arr_baseline_pow, 3);
        baseline_pow_std = nanstd(arr_baseline_pow, 0, 3);
        session_proc_lfp_out(i).baseline_mean = baseline_pow_mean;
        session_proc_lfp_out(i).baseline_std = baseline_pow_std;

        % store baseline power - commented for further inspection
%         baseline.sites(i).mean = baseline_pow_mean;
%         baseline.sites(i).std = baseline_pow_std;
        
        % save data
        site_lfp = session_proc_lfp_out(i);
        save(fullfile(cfg_tfs.session_results_fldr, ...
            [site_lfp.site_ID, '.mat']), 'site_lfp');

    end    
    
    fprintf(' done\n');
    fprintf('=============================================================\n');
end

