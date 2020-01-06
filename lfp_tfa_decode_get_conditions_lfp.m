function lfp_decode = lfp_tfa_decode_get_conditions_lfp( lfp_tfa_cfg )
%lfp_tfa_decode_get_conditions_lfp - Function to get raw LFP and LFP TFS
%for all the sessions for each condition specified in settings
%   lfp_decode = lfp_tfa_decode_get_conditions_lfp( lfp_tfa_cfg )
%   INPUTS:
%       lfp_tfa_cfg - struct containing the settings, see
%       settings/LFP_Decoding_Linus_8sessions/lfp_decoding_settings_Instr_control_IH_ISvsCS.m
%       for example
%       Required fields: 
%       analyse_epochs - windows to be analysed
%       session_info - info about sessions to be analysed
%   OUTPUTS: 
%       lfp_decode - struct containing session-wise LFP raw and TFS for all
%       trials from the specified conditions

% Average Evoked LFP across sites
lfp_decode = struct();
for i = 1:length(lfp_tfa_cfg.session_info)
    lfp_decode.raw_lfp.session(i).targets = {};
    lfp_decode.lfp_tfs.session(i).targets = {};
    lfp_decode.raw_lfp.session(i).trial = cell(1, size(lfp_tfa_cfg.analyse_epochs, 1));
    lfp_decode.raw_lfp.session(i).time = cell(1, size(lfp_tfa_cfg.analyse_epochs, 1));
    lfp_decode.raw_lfp.session(i).condition_idx = [];
    lfp_decode.lfp_tfs.session(i).trial = cell(1, size(lfp_tfa_cfg.analyse_epochs, 1));
    lfp_decode.lfp_tfs.session(i).time = cell(1, size(lfp_tfa_cfg.analyse_epochs, 1));
    lfp_decode.lfp_tfs.session(i).freq = cell(1, size(lfp_tfa_cfg.analyse_epochs, 1));
    lfp_decode.lfp_tfs.session(i).condition_idx = [];
    lfp_decode.lfp_pow.session(i).trial = cell(1, size(lfp_tfa_cfg.analyse_epochs, 1));
    lfp_decode.lfp_pow.session(i).freq = cell(1, size(lfp_tfa_cfg.analyse_epochs, 1));
    lfp_decode.lfp_pow.session(i).condition_idx = [];
end

results_folder = fullfile(lfp_tfa_cfg.root_results_fldr, 'LFP Decoding');
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end

% loop through each session to get lfp data from each site
for i = 1:length(lfp_tfa_cfg.session_info)
    session_info = lfp_tfa_cfg.session_info(i);
    proc_lfp_folder = session_info.proc_results_fldr;
    % read lfp data for each site
    site_lfp_files = dir(fullfile(proc_lfp_folder, 'site_lfp_*.mat'));
    nsites = 0;        
    for j = 1:length(site_lfp_files)
        filename = site_lfp_files(j).name;
        fprintf('Reading LFP data from %s\n', filename);
        load(fullfile(proc_lfp_folder, filename));
        % find the target for this site 
        nsites = nsites + 1;
        lfp_decode.raw_lfp.session(i).targets = ...
            [lfp_decode.raw_lfp.session(i).targets, site_lfp.target];
        lfp_decode.lfp_tfs.session(i).targets = ...
            [lfp_decode.lfp_tfs.session(i).targets, site_lfp.target];
              
    
        if nsites == 1
            % trial conditions to analyse
            conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg, ...
                {session_info.Preinj_blocks, session_info.Postinj_blocks});
            hs_labels = conditions(1).hs_labels;
            condition_trials_idx = cell(1, length(conditions) * length(hs_labels)); 
            cnd_idx = 0;
            class_conditions = conditions;
            class_conditions = rmfield(class_conditions, 'reach_hands');
            class_conditions = rmfield(class_conditions, 'reach_spaces');
            class_conditions = rmfield(class_conditions, 'hs_labels');
            for cnd = 1:length(conditions)
                for hs = 1:length(conditions(cnd).hs_labels)
                    cnd_idx = cnd_idx + 1;
                    % get trial indices for the given condition
                    cond_trials = lfp_tfa_get_condition_trials(site_lfp, conditions(cnd));                
                    % filter trials by hand-space labels
                    if ~strcmp(conditions(cnd).reach_hands{hs}, 'any')
                        cond_trials = cond_trials & ...
                            strcmp({site_lfp.trials.reach_hand}, ...
                            conditions(cnd).reach_hands{hs});
                    end
                    if ~strcmp(conditions(cnd).reach_spaces{hs}, 'any')
                        cond_trials = cond_trials & ...
                            strcmp({site_lfp.trials.reach_space}, ...
                            conditions(cnd).reach_spaces{hs});
                    end
                    % consider only non-noisy trials
                    cond_trials = cond_trials & ~[site_lfp.trials.noisy];
                    % find the condition index, state index, and hand-space index
                    % for this class
                    condition_trials_idx{cnd_idx} = find(cond_trials);
%                     % save a struct representing the condition
%                     if cnd_idx > 1
%                         class_conditions(cnd_idx) = class_conditions(cnd_idx - 1);
%                     end
                    class_conditions(cnd_idx).reach_hand = conditions(cnd).reach_hands{hs};
                    class_conditions(cnd_idx).reach_space = conditions(cnd).reach_spaces{hs};
                    class_conditions(cnd_idx).hs_label = conditions(cnd).hs_labels{hs};
                    class_conditions(cnd_idx).label = [conditions(cnd).label, ...
                        conditions(cnd).hs_labels{hs}];
                    class_conditions(cnd_idx).trials = condition_trials_idx{cnd_idx};
                    class_conditions(cnd_idx).type = conditions(cnd).type;
                    class_conditions(cnd_idx).effector = conditions(cnd).effector;
                    class_conditions(cnd_idx).choice = conditions(cnd).choice;
                    class_conditions(cnd_idx).perturbation = conditions(cnd).perturbation;
                    class_conditions(cnd_idx).perturbation_group = ...
                        conditions(cnd).perturbation_group;
                    class_conditions(cnd_idx).success = conditions(cnd).success;
                end                
            end
            % initialize trials data
            ntrials = length([condition_trials_idx{:}]);
            nepochs = size(lfp_tfa_cfg.analyse_epochs, 1);
            lfp_decode.raw_lfp.session(i).trial = cell(nepochs, ntrials);
            lfp_decode.lfp_tfs.session(i).trial = cell(nepochs, ntrials);
            lfp_decode.lfp_pow.session(i).trial = cell(nepochs, ntrials);
            lfp_decode.raw_lfp.session(i).conditions = class_conditions;
            lfp_decode.lfp_tfs.session(i).conditions = class_conditions;
        end
        
%         % for baseline normalization of tfs
%         cfg_baseline.method = lfp_tfa_cfg.baseline_method;
%         baseline_cnd_idx = [site_lfp.baseline.perturbation] == ...
%             lfp_tfa_cfg.baseline_perturbation & [site_lfp.baseline.choice] == ...
%             lfp_tfa_cfg.baseline_use_choice_trial;
%         cfg_baseline.mean = site_lfp.baseline(baseline_cnd_idx).pow_mean;
%         cfg_baseline.std = site_lfp.baseline(baseline_cnd_idx).pow_std;                
%         lfp_decode.lfp_tfs.session(i).baseline(nsites) = cfg_baseline;
        
        
        % loop through epochs to analyse
        for ep = 1:size(lfp_tfa_cfg.analyse_epochs, 1)
            epoch_refstate   = lfp_tfa_cfg.analyse_epochs{ep, 1};
            %epoch_name       = lfp_tfa_cfg.analyse_epochs{ep, 2};
            epoch_reftstart  = lfp_tfa_cfg.analyse_epochs{ep, 3};
            epoch_reftend    = lfp_tfa_cfg.analyse_epochs{ep, 4}; 
            
            condition_idx = zeros(1, length([condition_trials_idx{:}])); 
            
            for cnd = 1:length(condition_trials_idx)
                for t = condition_trials_idx{cnd}
                    trial_idx = find([condition_trials_idx{:}]==t);
                    condition_idx(trial_idx) = cnd;
                    % get timing information of epoch
                    states          = site_lfp.trials(t).states;
                    state_onset_t   = states([states(:).id] == ...
                        epoch_refstate).onset_t;
                    epoch_start_t   = states([states(:).id] == ...
                        epoch_refstate).onset_t + epoch_reftstart;
                    epoch_end_t     = states([states(:).id] == ...
                        epoch_refstate).onset_t + epoch_reftend;

                    % LFP data
                    trial_lfp = site_lfp.trials(t).lfp_data(:, ...
                        (site_lfp.trials(t).time >= epoch_start_t & ...
                        site_lfp.trials(t).time <= epoch_end_t));
                    trial_lfp_time = site_lfp.trials(t).time(1, ...
                        (site_lfp.trials(t).time >= epoch_start_t & ...
                        site_lfp.trials(t).time <= epoch_end_t)) - state_onset_t;
                    nsamples = size(trial_lfp, 2);
                    if nsites == 1
                        lfp_decode.raw_lfp.session(i).trial{ep, trial_idx} = ... 
                            cat(1, lfp_decode.raw_lfp.session(i).trial{ep, trial_idx}, ...
                            trial_lfp);
                        lfp_decode.raw_lfp.session(i).time{ep} = trial_lfp_time;
                    else
                        if size(lfp_decode.raw_lfp.session(i).trial{ep, trial_idx}, 2) < nsamples
                            nsamples = size(lfp_decode.raw_lfp.session(i).trial{ep, trial_idx}, 2);
                        end
                        trial_lfp = trial_lfp(:, 1:nsamples);
                        lfp_decode.raw_lfp.session(i).trial{ep, trial_idx} = ... 
                            cat(1, lfp_decode.raw_lfp.session(i).trial{ep, trial_idx}(:,1:nsamples), ...
                            trial_lfp);
                        lfp_decode.raw_lfp.session(i).time{ep} = trial_lfp_time(1:nsamples);
                    end                        
                    
                    % LFP time-freq spectrogram and power spectrum
                    trial_tfs = site_lfp.trials(t).tfs.powspctrm(1, :, ...
                        (site_lfp.trials(t).tfs.time >= epoch_start_t & ...
                        site_lfp.trials(t).tfs.time <= epoch_end_t));
                    % baseline normalization
%                     trial_tfs = lfp_tfa_baseline_normalization(...
%                         trial_tfs, cfg_baseline);
                    trial_tfs_time = site_lfp.trials(t).tfs.time(1, ...
                        (site_lfp.trials(t).tfs.time >= epoch_start_t & ...
                        site_lfp.trials(t).tfs.time <= epoch_end_t)) - state_onset_t;
                    trial_tfs_freq = site_lfp.trials(t).tfs.freq;
                    ntimebins = size(trial_tfs, 3);
                    if nsites == 1
                        lfp_decode.lfp_tfs.session(i).trial{ep, trial_idx} = ... 
                            cat(1, lfp_decode.lfp_tfs.session(i).trial{ep, trial_idx}, ...
                            trial_tfs);
                        lfp_decode.lfp_tfs.session(i).time{ep} = trial_tfs_time;
                        lfp_decode.lfp_tfs.session(i).freq{ep} = trial_tfs_freq;
                    else
                        if size(lfp_decode.lfp_tfs.session(i).trial{ep, trial_idx}, 3) < ntimebins
                            ntimebins = size(lfp_decode.lfp_tfs.session(i).trial{ep}{trial_idx}, 3);
                        end
                        trial_tfs = trial_tfs(:,:,1:ntimebins);
                        lfp_decode.lfp_tfs.session(i).trial{ep, trial_idx} = ... 
                            cat(1, ...
                            lfp_decode.lfp_tfs.session(i).trial{ep, trial_idx}(:,:,1:ntimebins), ...
                            trial_tfs);
                        lfp_decode.lfp_tfs.session(i).time{ep} = trial_tfs_time(1:ntimebins);
                    end                    
%                     lfp_decode.lfp_pow.session(i).trial{ep, trial_idx} = ... 
%                         cat(1, lfp_decode.lfp_pow.session(i).trial{ep, trial_idx}, ...
%                         sum(trial_tfs, 3));
                end
            end            
            
            if nsites == 1 && ep == 1
                lfp_decode.raw_lfp.session(i).condition_idx = condition_idx; 
                lfp_decode.lfp_tfs.session(i).condition_idx = ...
                    [lfp_decode.lfp_tfs.session(i).condition_idx, ...
                    condition_idx];
                lfp_decode.lfp_pow.session(i).condition_idx = ...
                    [lfp_decode.lfp_pow.session(i).condition_idx, ...
                    condition_idx];
            end
        end
    end
    % get same number of timebins for all trials
    for ep = 1:size(lfp_decode.raw_lfp.session(1).trial, 1)
        nsamples_rawlfp = min(cellfun(@length, lfp_decode.raw_lfp.session(i).trial(ep, :)));
        for t_idx = 1:size(lfp_decode.raw_lfp.session(i).trial, 2)
            lfp_decode.raw_lfp.session(i).trial{ep, t_idx} = ...
                lfp_decode.raw_lfp.session(i).trial{ep, t_idx}(:, 1:nsamples_rawlfp);
        end
        lfp_decode.raw_lfp.session(i).time{ep} = ...
            lfp_decode.raw_lfp.session(i).time{ep}(1:nsamples_rawlfp);
        lfp_decode.raw_lfp.session(i).epoch_name{ep} = ...
            lfp_tfa_cfg.analyse_epochs{ep, 2};
    end
    
    for ep = 1:size(lfp_decode.lfp_tfs.session(i).trial, 1)
        ndim_tfs = cellfun(@size, lfp_decode.lfp_tfs.session(i).trial(ep, :), 'uni', false);
        ndim_tfs = cat(1, ndim_tfs{:});
        if size(ndim_tfs, 2) > 2
            ntimebins_tfs = min(ndim_tfs(:,3));
        else
            ntimebins_tfs = 1;
        end
        for t_idx = 1:size(lfp_decode.lfp_tfs.session(i).trial, 2)
            lfp_decode.lfp_tfs.session(i).trial{ep, t_idx} = ...
                lfp_decode.lfp_tfs.session(i).trial{ep, t_idx}(:, :, 1:ntimebins_tfs);
        end
        lfp_decode.lfp_tfs.session(i).time{ep} = ...
            lfp_decode.lfp_tfs.session(i).time{ep}(1:ntimebins_tfs);
        lfp_decode.lfp_tfs.session(i).epoch_name{ep} = ...
            lfp_tfa_cfg.analyse_epochs{ep, 2};
    end
    
end



end

