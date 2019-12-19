function lfp_decode = lfp_tfa_decode_get_classes_lfp( lfp_tfa_cfg )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Average Evoked LFP across sites
lfp_decode = struct();
for cl = 1:length(lfp_tfa_cfg.decode.classes)
    lfp_decode.classes(cl).cfg_condition = lfp_tfa_cfg.decode.classes(cl);
    lfp_decode.classes(cl).label = lfp_tfa_cfg.decode.classes(cl).label;
    lfp_decode.classes(cl).raw_lfp = struct();
    lfp_decode.classes(cl).lfp_tfs = struct();
    lfp_decode.classes(cl).lfp_pow = struct();
    lfp_decode.classes(cl).raw_lfp.trial = cell(1, length(lfp_tfa_cfg.analyse_epochs));
    lfp_decode.classes(cl).raw_lfp.time = cell(1, length(lfp_tfa_cfg.analyse_epochs));
    lfp_decode.classes(cl).lfp_tfs.trial = cell(1, length(lfp_tfa_cfg.analyse_epochs));
    lfp_decode.classes(cl).lfp_tfs.time = cell(1, length(lfp_tfa_cfg.analyse_epochs));
    lfp_decode.classes(cl).lfp_tfs.freq = cell(1, length(lfp_tfa_cfg.analyse_epochs));
    lfp_decode.classes(cl).lfp_pow.trial = cell(1, length(lfp_tfa_cfg.analyse_epochs));
    lfp_decode.classes(cl).lfp_pow.freq = cell(1, length(lfp_tfa_cfg.analyse_epochs));
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
    first_site_found = false;        
    for j = 1:length(site_lfp_files)
        % find data from this site need to be used
        if isfield(lfp_tfa_cfg.session_info(i), 'sites') && ...
                ~empty(lfp_tfa_cfg.session_info(i).sites)
            if ~any(lfp_tfa_cfg.session_info(i).sites == j)
                continue;
            end
        end
        filename = site_lfp_files(j).name;
        fprintf('Reading LFP data from %s\n', filename);
        load(fullfile(proc_lfp_folder, filename));
        % find the class indices based on target hemisphere
        cl_idx = 1:length(lfp_tfa_cfg.decode.classes);
        if isfield(lfp_tfa_cfg.decode.classes, 'target')
            cl_idx = find(strcmp(...
                {lfp_tfa_cfg.decode.classes.target}, site_lfp.target));
        end    
    
        if ~isempty(cl_idx) && ~first_site_found            
            class_trials_idx = cell(1, length(cl_idx));         
            for cl = cl_idx
                % find the condition index, state index, and hand-space index
                % for this class
                class_trials_idx{cl} = lfp_tfa_decode_get_class_trials(...
                    site_lfp, lfp_tfa_cfg.decode.classes(cl));
            end
            first_site_found = true;
        end
        
        % for baseline normalization of tfs
        cfg_baseline.method = lfp_tfa_cfg.baseline_method;
        baseline_cnd_idx = [site_lfp.baseline.perturbation] == ...
            lfp_tfa_cfg.baseline_perturbation & [site_lfp.baseline.choice] == ...
            lfp_tfa_cfg.baseline_use_choice_trial;
        cfg_baseline.mean = site_lfp.baseline(baseline_cnd_idx).pow_mean;
        cfg_baseline.std = site_lfp.baseline(baseline_cnd_idx).pow_std;                
        
        for cl = cl_idx
            
            % loop through epochs to analyse
            for ep = 1:size(lfp_tfa_cfg.analyse_epochs, 1)
                epoch_refstate   = lfp_tfa_cfg.analyse_epochs{ep, 1};
                epoch_name       = lfp_tfa_cfg.analyse_epochs{ep, 2};
                epoch_reftstart  = lfp_tfa_cfg.analyse_epochs{ep, 3};
                epoch_reftend    = lfp_tfa_cfg.analyse_epochs{ep, 4};

                epoch_lfp = []; % raw lfp
                epoch_tfs = []; % lfp tfr
                epoch_psd = []; % power spectrum

                for t = class_trials_idx{cl}

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
                    nsamples = size(trial_lfp, 2);
                    if isempty(epoch_lfp)
                        epoch_lfp = trial_lfp;
                    else
                        if size(epoch_lfp, 2) < nsamples
                            nsamples = size(epoch_lfp, 2);
                        end
                        epoch_lfp = cat(1, epoch_lfp(:, 1:nsamples), ...
                            trial_lfp(:, 1:nsamples));
                    end
                    
                    % LFP time-freq spectrogram and power spectrum
                    trial_tfs = site_lfp.trials(t).tfs.powspctrm(1, :, ...
                        (site_lfp.trials(t).tfs.time >= epoch_start_t & ...
                        site_lfp.trials(t).tfs.time <= epoch_end_t));                     
                    ntimebins = size(trial_tfs, 3);
                    if isempty(epoch_tfs)
                        epoch_tfs = trial_tfs;
                    else
                        if size(epoch_tfs, 3) < ntimebins
                            ntimebins = size(epoch_tfs, 3);
                        end
                        epoch_tfs = cat(1, epoch_tfs(:,:,1:ntimebins), trial_tfs(:,:,1:ntimebins));                      
                    end
                    epoch_psd = cat(1, epoch_psd, sum(trial_tfs, 3));                    
                end
                epoch_lfp_time = site_lfp.trials(t).time(1, ...
                    (site_lfp.trials(t).time >= epoch_start_t & ...
                    site_lfp.trials(t).time <= epoch_end_t)) - state_onset_t;
                epoch_lfp_time = epoch_lfp_time(1, 1:nsamples);
                epoch_tfs_time = site_lfp.trials(t).tfs.time(1, ...
                    (site_lfp.trials(t).tfs.time >= epoch_start_t & ...
                    site_lfp.trials(t).tfs.time <= epoch_end_t)) - state_onset_t;
                epoch_tfs_time = epoch_tfs_time(1, 1:ntimebins);
                epoch_freq = site_lfp.trials(t).tfs.freq;
                % freq normalization
                epoch_psd = 10*log10(epoch_psd);
                % baseline normalization
                epoch_tfs = lfp_tfa_baseline_normalization(...
                   epoch_tfs, cfg_baseline);
                
                % save this into the lfp_decode struct
                % same number of lfp samples
                if isempty(lfp_decode.classes(cl).raw_lfp.time{ep})
                    lfp_decode.classes(cl).raw_lfp.trial{ep} = epoch_lfp;                        
                else
                    nsamples = size(lfp_decode.classes(cl).raw_lfp.trial{ep}, 2);
                    if nsamples > size(epoch_lfp, 2)
                        nsamples = size(epoch_lfp, 2);
                    end
                    lfp_decode.classes(cl).raw_lfp.trial{ep} = ...
                        cat(1, lfp_decode.classes(cl).raw_lfp.trial{ep}(:,1:nsamples), ...
                        epoch_lfp(:,1:nsamples));                
                end
                lfp_decode.classes(cl).raw_lfp.time{ep} = epoch_lfp_time(:,1:nsamples);                
                % same number of time bins of spectrogram                    
                if isempty(lfp_decode.classes(cl).lfp_tfs.time{ep})
                    lfp_decode.classes(cl).lfp_tfs.trial{ep} = epoch_tfs;
                else
                    ntimebins = size(lfp_decode.classes(cl).lfp_tfs.trial{ep}, 3);
                    if ntimebins > length(epoch_tfs_time)
                        ntimebins = length(epoch_tfs_time);
                    end
                    lfp_decode.classes(cl).lfp_tfs.trial{ep} = ...
                        cat(1, lfp_decode.classes(cl).lfp_tfs.trial{ep}(:,:,1:ntimebins), ...
                        epoch_tfs(:,:,1:ntimebins));                    
                end
                lfp_decode.classes(cl).lfp_tfs.time{ep} = epoch_tfs_time(:,1:ntimebins);
                lfp_decode.classes(cl).lfp_tfs.freq{ep} = epoch_freq;
                % power spectra
                if ~isempty(lfp_decode.classes(cl).lfp_pow.freq{ep})
                    lfp_decode.classes(cl).lfp_pow.freq{ep} = epoch_freq;
                end
                lfp_decode.classes(cl).lfp_pow.trial{ep} = ...
                    cat(1, lfp_decode.classes(cl).lfp_pow.trial{ep}, ...
                    epoch_psd);
                
            end            
        end
    end
end

% assign class labels
for cl = 1:length(lfp_decode.classes)
    lfp_decode.classes(cl).ntrials = ...
        size(lfp_decode.classes(cl).lfp_pow.trial{1}, 1);
    % train and test trials
    % Cross varidation (train: 50%, test: 50%)
    cv = cvpartition(lfp_decode.classes(cl).ntrials,'HoldOut',0.5);
    idx = cv.test;
    % Separate to training and test data
    lfp_decode.classes(cl).train_trials = ~idx;
    lfp_decode.classes(cl).test_trials = idx;
    lfp_decode.classes(cl).raw_lfp.actual_labels = repmat(...
        {lfp_decode.classes(cl).label}, 1, ...
        lfp_decode.classes(cl).ntrials);
    lfp_decode.classes(cl).raw_lfp.predicted_labels = repmat(...
        {''}, 1, size(lfp_decode.classes(cl).raw_lfp.trial, 2), ...
        lfp_decode.classes(cl).ntrials);
    lfp_decode.classes(cl).lfp_tfs.actual_labels = repmat(...
        {lfp_decode.classes(cl).label}, 1, ...
        lfp_decode.classes(cl).ntrials);
    lfp_decode.classes(cl).lfp_tfs.predicted_labels = repmat(...
        {''}, size(lfp_decode.classes(cl).lfp_tfs.trial, 2), ...
        lfp_decode.classes(cl).ntrials);
    lfp_decode.classes(cl).lfp_pow.actual_labels = repmat(...
        {lfp_decode.classes(cl).label}, 1, ...
        lfp_decode.classes(cl).ntrials);
    lfp_decode.classes(cl).lfp_pow.predicted_labels = repmat(...
        {''}, size(lfp_decode.classes(cl).lfp_pow.trial, 2), ...
        lfp_decode.classes(cl).ntrials);
end

% save file
save(fullfile(results_folder, 'lfp_decode_classes.mat'), 'lfp_decode');

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

