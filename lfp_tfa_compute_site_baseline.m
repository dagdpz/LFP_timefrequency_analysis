function [ site_lfp ] = lfp_tfa_compute_site_baseline( site_lfp, session_info, lfp_tfa_cfg )

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
    task_type = unique([lfp_tfa_cfg.compare.types]);
    task_effector = unique([lfp_tfa_cfg.compare.effectors]);
    perturbation = unique([lfp_tfa_cfg.compare.perturbations]);
    perturbation_blocks = unique([site_lfp.trials.perturbation]);
    %choice = unique([site_lfp.trials.choice_trial]);
    choice = unique([lfp_tfa_cfg.compare.choice_trials]);
    % loop through each baseline condition
    cnd = 0;
    
    for k = task_type
        for e = task_effector
            for p = perturbation
                for c = choice   
                    
                    if p == 0
                        blocks = session_info.Preinj_blocks;
                    elseif p == 1
                        blocks = session_info.Postinj_blocks;
                        if strcmp(blocks, 'all')
                            blocks = perturbation_blocks(perturbation_blocks ~= 0);
                        elseif strcmp(blocks, 'allbutfirst')
                            blocks = perturbation_blocks(perturbation_blocks ~= 0);
                            blocks = blocks(2:end);
                        end
                    end
                    
                    cnd = cnd+1;
                    baseline(cnd).type = k;
                    baseline(cnd).effector = e;
                    baseline(cnd).perturbation = p;
                    baseline(cnd).choice = c;
                    % struct to store lfp power during baseline period for each trial
                    baseline_pow = cell(1, length(site_lfp.trials));
                    % loop through trial
                    for t = 1:length( site_lfp.trials )
                        trial = site_lfp.trials(t);
                        
                        % whether this trial should be considered for baseline calculation
                        consider_trial = ~(trial.noisy) && trial.completed;
                        % based on task type
                        if ~isnan(k) && ~isinf(k)
                            consider_trial = consider_trial & sum(trial.type == k);
                        end
                        % based on task effector
                        if ~isnan(e) && ~isinf(e)
                            consider_trial = consider_trial & sum(trial.effector == e);
                        end
                        % based on perturbation
                        if ~isnan(p) && ~isinf(p)
                            consider_trial = consider_trial & sum(trial.perturbation == blocks);
                        end
                        % based on choice
                        if ~isnan(c) && ~isinf(c)
                            consider_trial = consider_trial & sum(trial.choice_trial == c);
                        end
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
                    clear('arr_baseline_pow'); 
                end
            end
        end
    end

    clear('baseline_pow');
    
    % baseline normalized power
    if isnan(lfp_tfa_cfg.baseline_perturbation)
        baseline_cnd_idx = isnan([site_lfp.baseline.perturbation]);
    end
    if isnan(lfp_tfa_cfg.baseline_use_choice_trial)
        baseline_cnd_idx = baseline_cnd_idx & isnan([site_lfp.baseline.choice]);
    end
    if ~isnan(lfp_tfa_cfg.baseline_perturbation) && ~isnan(lfp_tfa_cfg.baseline_use_choice_trial)
        baseline_cnd_idx = [site_lfp.baseline.perturbation] == ...
            lfp_tfa_cfg.baseline_perturbation & [site_lfp.baseline.choice] == ...
            lfp_tfa_cfg.baseline_use_choice_trial;
    end
    baseline_mean = site_lfp.baseline(baseline_cnd_idx).pow_mean;
    baseline_std = site_lfp.baseline(baseline_cnd_idx).pow_std;
    
%     trials_tfs = cat(1, site_lfp.trials(:).tfs);
%     site_lfp_pow = cat(3, trials_tfs.powspctrm);
%     site_lfp_pow_norm = (site_lfp_pow - repmat(baseline_mean, [1, 1, size(site_lfp_pow, 3)])) ./...
%         repmat(baseline_std, [1, 1, size(site_lfp_pow, 3)]);
%     
%     % plot
%     trial_pow_timestep = nanmean(diff(trials_tfs(1).time));
%     x = 1:trial_pow_timestep: trial_pow_timestep * (size(arr_concat_lfp_pow, 3) - 1);
%     y = trials_tfs(1).freq;
%     subplot(3,2,[5 6])
%     imagesc(x, y, squeeze(site_lfp_pow_norm), [-1 1]);
%     axis xy
%     title('LFP power spectrogram');
%     results_folder_noise = fullfile(lfp_tfa_cfg.noise.results_folder, 'LFP_Noise_Rejection');
%     saveas(hnoise, fullfile(results_folder_noise, [site_lfp.site_ID '_concat_raw_and_deriv_LFP.png']));
%             
%     clear('trials_tfs', 'site_lfp_pow', 'site_lfp_pow_norm');
    
    %end    
    
    fprintf(' done\n');
    %fprintf('=============================================================\n');
end


