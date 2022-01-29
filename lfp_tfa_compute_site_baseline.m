function [ site_lfp ] = lfp_tfa_compute_site_baseline( site_lfp, session_info, lfp_tfa_cfg )
% lfp_tfa_compute_site_baseline - Computes mean and standard deviation of
% LFP power for a site for each frequency during the defined baseline period
% for different perturbation and choice conditions
%
% USAGE:
%	site_lfp = lfp_tfa_compute_site_baseline( site_lfp, session_info, lfp_tfa_cfg )
%
% INPUTS:
%		site_lfp        - struct containing LFP time frequency spectrogram 
%       average for a site. i.e., the output of lfp_tfa_compute_site_tfr
%       session_info    - struct containing information about the session
%       being analysed, see settings/lfp_tfa_settings_example
%           Required fields:
%               Preinj_blocks - blocks to be considered as pre-injection,
%               can be an integer or integer array
%               Postinj_blocks - blocks to be considered as post-injection,
%               can be an integer, integer array, 'all' or 'allbutfirst'
%       lfp_tfa_cfg      - settings for baseline computation, see 
%       settings/lfp_tfa_settings_example
%           Required Fields: 
%               1. compare.perturbations    - perturbation conditions to be
%               compared, see settings/lfp_tfa_settings_example
%               2. baseline_refstate           - reference state for baseline
%               3. baseline_period             - period around the reference state to be
%               considered for baseline computation
% OUTPUTS:
%		site_lfp         - structure containing mean and standard deviation
%		of LFP power at each frequency of interest for different
%		perturbation and choice conditions
%           baseline.pow_mean   - mean LFP spectral power during the
%           baseline period
%           baseline.pow_std    - standard deviation of  LFP spectral power in baseline
%           period
%
% See also lfp_tfa_compute_site_tfr, lfp_tfa_process_lfp, 
% settings/lfp_tfa_settings_example
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
    task_type = inf;
    if isfield(lfp_tfa_cfg, 'baseline_use_type') && ...
            ~isempty(lfp_tfa_cfg.baseline_use_type)
        task_type = lfp_tfa_cfg.baseline_use_type;%unique([lfp_tfa_cfg.compare.types]);
    end
    task_effector = inf;
    if isfield(lfp_tfa_cfg, 'baseline_use_effector') && ...
            ~isempty(lfp_tfa_cfg.baseline_use_effector)
        task_effector = lfp_tfa_cfg.baseline_use_effector;%unique([lfp_tfa_cfg.compare.effectors]);
    end
    
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
                        consider_trial = ~(trial.noisy); % && trial.completed; MP debug
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
                        if strcmp(lfp_tfa_cfg.baseline_ref_state, '') && strcmp(lfp_tfa_cfg.baseline_ref_period , 'trial')
                            ref_period = trial.trialperiod;
                        elseif ~strcmp(lfp_tfa_cfg.baseline_ref_state, '')
                            ref_period = trial.states([trial.states.id] == lfp_tfa_cfg.baseline_ref_state) .onset_t + lfp_tfa_cfg.baseline_ref_period;
                        end
                        if consider_trial 
                            baseline_pow{t} = trial.tfs.powspctrm(1, :, trial.tfs.time >= ref_period(1) & trial.tfs.time <= ref_period(2));
                        else
                            baseline_pow{t} = NaN(size(trial.tfs.powspctrm,2),sum(trial.tfs.time >= ref_period(1) & trial.tfs.time <= ref_period(2));
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
    
    fprintf(' done\n');

end


