function [ session_proc_lfp_out ] = lfp_tfa_compute_baseline_power( session_lfp_in, cfg_tfs )

% computeBaseline - Computes mean and stddev of baseline power according to
% the given configuration
%
% USAGE:
%	[ sites_lfp_folder, baseline ] = lfp_tfa_compute_baseline_power( sites_lfp_folder, cfg_tfs )
%
% INPUTS:
%		sites_lfp_folder      - folder containing LFP data strcture for all
%		sites
%       cfg_tfs         - configuration structure for baseline computation
%           baseline_refstate   : reference state for baseline
%           baseline_period     : period around the reference state to be
%           considered for baseline computation
%           baseline_block      : block to be considered for baseline
%           calculation
%           use_choice_trial    : 1/0 whether to consider choice trials (1)
%           or instructed trials (0) for baseline computation
% OUTPUTS:
%		sites_lfp_folder      - folder to which the results are updated,
%		same as input data folder
%       baseline        - struct which saves the baseline mean and stddev for each site and 
%       configuration used for baseline calculation
% See also lfp_tfa_reject_noisy_trials

    % struct for storing results
    session_proc_lfp_out = session_lfp_in;
    %baseline.cfg = cfg_tfs;
    fprintf('=============================================================\n');
    fprintf('Computing baseline ...\n');
    
    % folder for saving results
    %sessionName = states_lfp(1).session;
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
            consider_trial = ~(trial.noisy) & (trial.block == cfg_tfs.baseline_block) ...
                & (trial.choice_trial == cfg_tfs.use_choice_trial);
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
        session_proc_lfp_out(i).baseline_mean = baseline_pow_mean;
        session_proc_lfp_out(i).baseline_std = baseline_pow_std;
        
%         baseline.sites(i).mean = baseline_pow_mean;
%         baseline.sites(i).std = baseline_pow_std;
        
        % save data
        site_lfp = session_proc_lfp_out(i);
        save(fullfile(cfg_tfs.session_results_fldr, ...
            [site_lfp.site_ID, '.mat']), 'site_lfp');
    end
    
    % save results
    %save(fullfile(results_folder, 'states_lfp.mat'), 'states_lfp');
    %save(fullfile(results_folder, 'baseline.mat'), 'baseline');
    session_proc_lfp_out(i) = site_lfp;
    fprintf(' done\n');
    fprintf('=============================================================\n');
end
