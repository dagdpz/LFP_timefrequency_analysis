function lfp_decode = lfp_tfa_decode_get_sessions_lfp( lfp_tfa_cfg )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Average Evoked LFP across sites
lfp_decode = struct();
for i = 1:length(lfp_tfa_cfg.session_info)
    lfp_decode.raw_lfp.session(i).trial = cell(1, length(lfp_tfa_cfg.analyse_epochs));
    lfp_decode.raw_lfp.session(i).time = cell(1, length(lfp_tfa_cfg.analyse_epochs));
    lfp_decode.raw_lfp.session(i).classes = [];
    lfp_decode.lfp_tfs.session(i).trial = cell(1, length(lfp_tfa_cfg.analyse_epochs));
    lfp_decode.lfp_tfs.session(i).time = cell(1, length(lfp_tfa_cfg.analyse_epochs));
    lfp_decode.lfp_tfs.session(i).freq = cell(1, length(lfp_tfa_cfg.analyse_epochs));
    lfp_decode.lfp_tfs.session(i).classes = [];
    lfp_decode.lfp_pow.session(i).trial = cell(1, length(lfp_tfa_cfg.analyse_epochs));
    lfp_decode.lfp_pow.session(i).freq = cell(1, length(lfp_tfa_cfg.analyse_epochs));
    lfp_decode.lfp_pow.session(i).classes = [];
end

results_folder = fullfile(lfp_tfa_cfg.root_results_fldr, 'LFP Decoding');
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end

% loop through each session to get lfp data from each site
for i = 1:length(lfp_tfa_cfg.session_info)
    proc_lfp_folder = lfp_tfa_cfg.session_info(i).proc_results_fldr;
    % read lfp data for each site
    site_lfp_files = dir(fullfile(proc_lfp_folder, 'site_lfp_*.mat'));
    nsites = 0;        
    for j = 1:length(site_lfp_files)
        filename = site_lfp_files(j).name;
        fprintf('Reading LFP data from %s\n', filename);
        load(fullfile(proc_lfp_folder, filename));
        % find the class indices based on target hemisphere
        cl_idx = 1:length(lfp_tfa_cfg.decode.classes);
        if isfield(lfp_tfa_cfg.decode.classes, 'target')
            cl_idx = find(strcmp(...
                {lfp_tfa_cfg.decode.classes.target}, site_lfp.target));
        end 
        % find the target for this site 
        if isfield(lfp_tfa_cfg.decode.classes, 'target') 
            if ~any(strcmp({lfp_tfa_cfg.decode.classes.target}, site_lfp.target))
                continue;
            else
                nsites = nsites + 1;
            end
        end    
    
        if nsites == 1     
            class_trials_idx = cell(1, length(lfp_tfa_cfg.decode.classes));         
            for cl = cl_idx
                % find the condition index, state index, and hand-space index
                % for this class
                class_trials_idx{cl} = lfp_tfa_decode_get_class_trials(...
                    site_lfp, lfp_tfa_cfg.decode.classes(cl));
            end
            % initialize trials data
            ntrials = length([class_trials_idx{:}]);
            nepochs = size(lfp_tfa_cfg.analyse_epochs, 1);
            lfp_decode.raw_lfp.session(i).trial = cell(nepochs, ntrials);
            lfp_decode.lfp_tfs.session(i).trial = cell(nepochs, ntrials);
            lfp_decode.lfp_pow.session(i).trial = cell(nepochs, ntrials);
        end
        
        % for baseline normalization of tfs
        cfg_baseline.method = lfp_tfa_cfg.baseline_method;
        baseline_cnd_idx = [site_lfp.baseline.perturbation] == ...
            lfp_tfa_cfg.baseline_perturbation & [site_lfp.baseline.choice] == ...
            lfp_tfa_cfg.baseline_use_choice_trial;
        cfg_baseline.mean = site_lfp.baseline(baseline_cnd_idx).pow_mean;
        cfg_baseline.std = site_lfp.baseline(baseline_cnd_idx).pow_std;                
        lfp_decode.lfp_tfs.session(i).baseline(nsites) = cfg_baseline;
        
        
        % loop through epochs to analyse
        for ep = 1:size(lfp_tfa_cfg.analyse_epochs, 1)
            epoch_refstate   = lfp_tfa_cfg.analyse_epochs{ep, 1};
            %epoch_name       = lfp_tfa_cfg.analyse_epochs{ep, 2};
            epoch_reftstart  = lfp_tfa_cfg.analyse_epochs{ep, 3};
            epoch_reftend    = lfp_tfa_cfg.analyse_epochs{ep, 4}; 
            
            class_labels = zeros(1, length([class_trials_idx{:}])); 
            
            for cl = cl_idx
                for t = class_trials_idx{cl}
                    trial_idx = find([class_trials_idx{:}]==t);
                    class_labels(trial_idx) = cl;
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
                    trial_tfs = lfp_tfa_baseline_normalization(...
                        trial_tfs, cfg_baseline);
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
            
%             epoch_lfp_time = site_lfp.trials(t).time(1, ...
%                 (site_lfp.trials(t).time >= epoch_start_t & ...
%                 site_lfp.trials(t).time <= epoch_end_t)) - state_onset_t;
%             epoch_lfp_time = epoch_lfp_time(1, 1:nsamples);
%             epoch_tfs_time = site_lfp.trials(t).tfs.time(1, ...
%                 (site_lfp.trials(t).tfs.time >= epoch_start_t & ...
%                 site_lfp.trials(t).tfs.time <= epoch_end_t)) - state_onset_t;
%             epoch_tfs_time = epoch_tfs_time(1, 1:ntimebins);
%             epoch_freq = site_lfp.trials(t).tfs.freq;
%             % freq normalization
%             epoch_psd = 10*log10(epoch_psd);
%             % baseline normalization
%             epoch_tfs = lfp_tfa_baseline_normalization(...
%                 epoch_tfs, cfg_baseline);
% 
%             % save this into the lfp_decode struct
%             % same number of lfp samples
%             if nsites == 1
%                 lfp_decode.raw_lfp.session(i).trial{ep} = {epoch_lfp};                        
%             else
%                 nsamples = size(lfp_decode.raw_lfp.session(i).trial{ep}{nsites}, 2);
%                 if nsamples > size(epoch_lfp, 2)
%                     nsamples = size(epoch_lfp, 2);
%                 end
%                 lfp_decode.raw_lfp.session(i).trial{ep} = ...
%                     [lfp_decode.raw_lfp.session(i).trial{ep}(:,1:nsamples), ...
%                     epoch_lfp(:,1:nsamples)];                
%             end
%             lfp_decode.raw_lfp.session(i).time{ep} = epoch_lfp_time(:,1:nsamples); 
            if nsites == 1 && ep == 1
                lfp_decode.raw_lfp.session(i).classes = class_labels; 
                lfp_decode.lfp_tfs.session(i).classes = ...
                    [lfp_decode.lfp_tfs.session(i).classes, ...
                    class_labels];
                lfp_decode.lfp_pow.session(i).classes = ...
                    [lfp_decode.lfp_pow.session(i).classes, ...
                    class_labels];
            end
%             % same number of time bins for all sites                    
%             if nsites == 1
%                 lfp_decode.lfp_tfs.session(i).trial{ep} = {epoch_tfs};
%             else
%                 ntimebins = size(lfp_decode.lfp_tfs.session(i).trial{ep}, 3);
%                 if ntimebins > length(epoch_tfs_time)
%                     ntimebins = length(epoch_tfs_time);
%                 end
%                 lfp_decode.lfp_tfs.session(i).trial{ep} = ...
%                     [lfp_decode.lfp_tfs.session(i).trial{ep}(:,:,1:ntimebins), ...
%                     epoch_tfs(:,:,1:ntimebins)];                    
%             end
%             lfp_decode.lfp_tfs.session(i).time{ep} = epoch_tfs_time(:,1:ntimebins);
%             lfp_decode.lfp_tfs.session(i).freq{ep} = epoch_freq;
%             if nsites == 1 && ep == 1
%                 
%             end
%             % power spectra
%             if nsites == 1
%                 lfp_decode.lfp_pow.session(i).trial{ep} = {epoch_psd};
%                 lfp_decode.lfp_pow.session(i).freq{ep} = epoch_freq;
%             else
%                 lfp_decode.lfp_pow.session(i).trial{ep} = ...
%                     [lfp_decode.lfp_pow.session(i).trial{ep}, ...
%                     epoch_psd];
%             end
%             if nsites == 1 && ep == 1
%                 
%             end            
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

% assign class labels
% for i = 1:length(lfp_decode.lfp_pow.session)
% 
%     lfp_decode.raw_lfp.sessions(i).pred_classes = ...
%         nan(size(lfp_decode.raw_lfp.session(i).classes));
%     lfp_decode.raw_lfp.sessions(i).actual_labels = ...
%         {lfp_tfa_cfg.decode.classes(lfp_decode.raw_lfp.session(i).classes).label};
%     lfp_decode.raw_lfp.sessions(i).predicted_labels = repmat(...
%         {''}, size(lfp_decode.raw_lfp.session(i).classes));
%     lfp_decode.lfp_tfs.sessions(i).pred_classes = ...
%         nan(size(lfp_decode.lfp_tfs.session(i).classes));
%     lfp_decode.lfp_tfs.sessions(i).actual_labels = ...
%         {lfp_tfa_cfg.decode.classes(lfp_decode.lfp_tfs.sessions(i).classes).label};
%     lfp_decode.lfp_tfs.session(i).predicted_labels = repmat(...
%         {''}, size(lfp_decode.raw_lfp.session(i).classes));
%     lfp_decode.lfp_pow.sessions(i).pred_classes = ...
%         nan(size(lfp_decode.lfp_pow.session(i).classes));
%     lfp_decode.lfp_pow.session(i).actual_labels = ...
%         {lfp_tfa_cfg.decode.classes(lfp_decode.lfp_pow.session(i).classes).label};
%     lfp_decode.lfp_pow.session(i).predicted_labels = repmat(...
%         {''}, size(lfp_decode.raw_lfp.session(i).classes));
% end

% save file
save(fullfile(results_folder, 'lfp_decode_data.mat'), 'lfp_decode');

%             for cn = class_trials_idx
%                 lfp_decode.classes(cl).raw_lfp.nsites = 0;
%                 for i = 1:length(lfp_evoked.session) 
%                     for j = 1:length(lfp_evoked.session(i).sites)
%                         if isfield(lfp_tfa_cfg.decode.classes(cl), 'target') && ...
%                                 ~strcmp(lfp_evoked.session(i).sites(j).target, lfp_tfa_cfg.decode.classes(cl).target)
%                             continue;
%                         end
%                         if ~lfp_evoked.session(i).sites(j).use_for_avg
%                             continue;
%                         end
% 
%                         if ~isempty(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked) && ... 
%                             isfield(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked, 'lfp')
%                             if isfield(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st_idx, hs_idx), 'lfp') ...
%                                     && ~isempty(lfp_evoked.session(i).sites(j).condition(cn).hs_tuned_evoked(st_idx, hs_idx).lfp)
%                                 lfp_decode.classes(cl).raw_lfp.nsites = ...
%                                     lfp_decode.classes(cl).raw_lfp.nsites + 1;
%                                 if lfp_decode.classes(cl).raw_lfp.nsites == 1%~isfield(sessions_avg.condition(cn).hs_tuned_evoked, 'mean')
% 
%                                     if isfield(lfp_evoked.session(i).sites(j).condition(cnd_idx).hs_tuned_evoked(st_idx, hs_idx(1)), 'state') && ...
%                                             isfield(lfp_evoked.session(i).sites(j).condition(cnd_idx).hs_tuned_evoked(st_idx, hs_idx(1)), 'state_name')
%                                         lfp_decode.classes(cl).raw_lfp.state ...
%                                             = lfp_evoked.session(i).sites(j).condition(cnd_idx).hs_tuned_evoked(st_idx, hs_idx(1)).state;
%                                         lfp_decode.classes(cl).raw_lfp.state_name ...
%                                             = lfp_evoked.session(i).sites(j).condition(cnd_idx).hs_tuned_evoked(st_idx, hs_idx(1)).state_name;
%                                     end
% 
%                                     lfp_decode.classes(cl).raw_lfp.time ...
%                                         = lfp_evoked.session(i).sites(j).condition(cnd_idx).hs_tuned_evoked(st_idx, hs_idx(1)).time;
%                                     lfp_decode.classes(cl).raw_lfp.trial ...
%                                         = cat(1, lfp_evoked.session(i).sites(j).condition(cnd_idx).hs_tuned_evoked(st_idx, hs_idx).lfp);
% 
%                                 else
% 
%                                     nsamples = length(lfp_decode.classes(cn).raw_lfp.time);
%                                     % average same number of time bins
%                                     if nsamples > length(lfp_evoked.session(i).sites(j).condition(cnd_idx).hs_tuned_evoked(st_idx, hs_idx).time)
%                                         nsamples = length(lfp_evoked.session(i).sites(j).condition(cnd_idx).hs_tuned_evoked(st_idx, hs_idx).time);
%                                     end
%                                     lfp_decode.classes(cl).raw_lfp.time = ...
%                                         lfp_evoked.session(i).sites(j).condition(cnd_idx).hs_tuned_evoked(st_idx, hs_idx(1)).time(1:nsamples);
%                                     for hs = hs_idx
%                                         lfp_decode.classes(cl).raw_lfp.trial ...
%                                             = cat(1, ...
%                                             lfp_decode.classes(cl).raw_lfp.trial(:, 1:nsamples), ...
%                                             (lfp_evoked.session(i).sites(j).condition(cnd_idx).hs_tuned_evoked(st_idx, hs).lfp(:,1:nsamples)));
%                                     end
% 
%                                 end
%                             end
% 
%                         end
%                     end
%                 end
%             end
% 
%             lfp_decode.classes(cl).raw_lfp.ntrials = size(lfp_decode.classes(cl).raw_lfp.trial, 1);
% 
% 
%         end
%     end

end

