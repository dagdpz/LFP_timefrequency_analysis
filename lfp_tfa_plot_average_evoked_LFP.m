function [ cond_based_evoked ] = lfp_tfa_plot_average_evoked_LFP( sites_lfp_folder, analyse_states, lfp_tfa_cfg ) 

% lfp_tfa_plot_average_evoked_LFP  - plots average evoked LFP for
% different hand-space tuning conditions for each site and across sites
%
% USAGE:
%	[ cond_based_evoked ] = lfp_tfa_plot_average_evoked_LFP( sites_lfp_folder, analyse_states, lfp_tfa_cfg )
%
% INPUTS:
%		sites_lfp_folder  	- folder containing lfp data, output from
%		lfp_tfa_compute_baseline
%       analyse_states  - cell array containing ids of states to be
%       analysed
%       lfp_tfa_cfg     - struct containing configuration for TFR 
%           Required fields:
%           trial_condition.blocks              - blocks to be analysed, 
%           leave empty to analyse each block separately
%           results_folder                      - folder to save results
%
% OUTPUTS:
%		cond_based_evoked	- output structure which saves the average evoked LFP for  
%                         trials of a given condition for different handspace 
%                         tunings and periods around the states analysed
%                         same datastructure as input ft_data_sites, but
%                         with additional field to store condition-wise lfp
%                         time freq response		
%
%
%
% See also lfp_tfa_compute_baseline, lfp_tfa_define_states
    
    % suppress warning for xticklabel
    warning ('off', 'MATLAB:hg:willberemoved');

    % make a folder to save figures
    %sessionName = states_lfp(1).session;
    results_folder_tfr = fullfile(lfp_tfa_cfg.results_folder, 'Condition_based_Evoked_LFP');
    if ~exist(results_folder_tfr, 'dir')
        mkdir(results_folder_tfr);
    end
    
    states_lfp = struct();
    % load site lfps
    site_lfp_files = dir([sites_lfp_folder, '\*.mat']);
    for i = 1:length( site_lfp_files )
        load(fullfile(sites_lfp_folder, site_lfp_files(i).name));
        if i == 1
            states_lfp = site_lfp;
        else
            states_lfp = [states_lfp, site_lfp];
        end
    end
    
    recorded_hemispace = unique([states_lfp.recorded_hemispace]);
    choice = unique([states_lfp(1).trials.choice_trial]);
    perturbation = unique([states_lfp(1).trials.perturbation]);
    blocks = unique([states_lfp(1).trials.block]);
    
    % create conditions
    cfg_conditions = struct();
       
    i = 0;
    for rec_hem = recorded_hemispace        
        for c = choice
            for b = blocks
                i = i + 1;
                cfg_conditions(i).recorded_hemispace = rec_hem;
                cfg_conditions(i).choice = c;
                cfg_conditions(i).block = b;
                cfg_conditions(i).perturbation = ...
                    unique([states_lfp(1).trials([states_lfp(1).trials.block] == b).perturbation]);
                cond_label = [];
                if cfg_conditions(i).recorded_hemispace == 'L'
                    cond_label = [cond_label 'Left_hemisphere_'];
                else
                    cond_label = [cond_label 'Right_hemisphere_'];
                end
                if cfg_conditions(i).choice == 0
                    cond_label = [cond_label 'Instructed_'];
                else
                    cond_label = [cond_label 'Choice_'];
                end
                if cfg_conditions(i).perturbation == 0
                    cond_label = [cond_label 'Control_'];
                else
                    cond_label = [cond_label 'Inactivation_'];
                end
                cond_label = [cond_label, 'Block_', num2str(cfg_conditions(i).block)];
                cfg_conditions(i).label = cond_label;
                
                % create a folder for storing results for this condition
                cfg_conditions(i).results_folder = fullfile(results_folder_tfr); %, cfg_conditions(i).label
                if ~exist(cfg_conditions(i).results_folder, 'dir')
                    mkdir(cfg_conditions(i).results_folder)
                end
                                
            end
        end
    end
    
    
    if isfield(lfp_tfa_cfg, 'add_conditions')
        for c = 1:length(lfp_tfa_cfg.add_conditions)
            if ~isempty(lfp_tfa_cfg.add_conditions(c))
                if ~isempty(lfp_tfa_cfg.add_conditions(c).blocks)
                    for rec_hem = recorded_hemispace        
                        for ch = choice
                            i = i + 1;
                            if strcmp(lfp_tfa_cfg.add_conditions(c).blocks, 'inactivation')
                                cfg_conditions(i).block = blocks(perturbation ~= 0);
                                cfg_conditions(i).perturbation = 1;
                            else
                                cfg_conditions(i).block = lfp_tfa_cfg.add_conditions(c).blocks;
                                cfg_conditions(i).perturbation = perturbation(blocks == lfp_tfa_cfg.add_conditions(c).blocks(1));
                                
                            end                    
                            cfg_conditions(i).choice = ch;
                            if isfield(lfp_tfa_cfg.add_conditions(c), 'perturbation')
                                cfg_conditions(i).perturbation = lfp_tfa_cfg.add_conditions(c).perturbation;
                            end
                            cfg_conditions(i).recorded_hemispace = rec_hem;
                            cond_label = [];
                            if cfg_conditions(i).recorded_hemispace == 'L'
                                cond_label = [cond_label 'Left_hemispace_'];
                            else
                                cond_label = [cond_label 'Right_hemispace_'];
                            end
                            if cfg_conditions(i).choice == 0
                                cond_label = [cond_label 'Instructed_'];
                            else
                                cond_label = [cond_label 'Choice_'];
                            end
                            if cfg_conditions(i).perturbation == 0
                                cond_label = [cond_label 'Control_'];
                            else
                                cond_label = [cond_label 'Inactivation_'];
                            end
                            cond_label = [cond_label, 'Block_', num2str(cfg_conditions(i).block)];
                            cfg_conditions(i).label = cond_label;
                            
                            % create a folder for storing results for this condition
                            cfg_conditions(i).results_folder = fullfile(results_folder_tfr); %, cfg_conditions(i).label
                            if ~exist(cfg_conditions(i).results_folder, 'dir')
                                mkdir(cfg_conditions(i).results_folder)
                            end
                        end
                    end
                end
            end
        end
    end
    
       
    % condition based TFS
    cond_based_evoked = struct();
    for cn = 1:length(cfg_conditions)
        
        cond_based_evoked(cn).label = cfg_conditions(cn).label;
        cond_based_evoked(cn).cfg_cond = cfg_conditions(cn);
        
        % states to be analysed
    %     analyse_states = {states_lfp(1).states.name};

        % hand-space tuning of LFP
        hs_labels = unique({states_lfp(1).trials.hndspc_lbl});
        
        % current trial condition analysed
        cfg_condition = cfg_conditions(cn);
        
        
        % cell array to store time frequency average across sites
        TFR_avg = cell(length(analyse_states),length(hs_labels));

        % loop through each site
        nsites = length(states_lfp);
        
        cond_based_evoked(cn).sites = struct();

        % loop through each site
        for i = 1:length(states_lfp)            
            
            %cond_based_tfs(cn).sites(i) = struct();
            % consider site based on recorded hemispace
            if strcmp(states_lfp(i).recorded_hemispace, cfg_conditions(cn).recorded_hemispace) 
                cond_based_evoked(cn).sites(i).site_ID = states_lfp(i).site_ID;
                cond_based_evoked(cn).sites(i).session = states_lfp(i).session;
                cond_based_evoked(cn).sites(i).target = states_lfp(i).target;
                % make a struct for concatenating TFR for all states
                cond_based_evoked(cn).sites(i).tfs_avg_site = struct(); 
                % struct to store evoked LFP
                cond_based_evoked(cn).sites(i).site_evoked_lfp = struct();
                % struct to store LFP power spectrum
                cond_based_evoked(cn).sites(i).site_lfp_psd = struct();
                %cond_based_tfs(i).tfs_avg_site.powspctrm = cell(length(analyse_states),length(hs_labels));
                cond_based_evoked(cn).sites(i).ntrials = zeros(1,length(hs_labels));

                for hs = 1:length(hs_labels)
                    cond_trials = zeros(1, length(states_lfp(i).trials));
                    % get the trials for given condition and this hs label
%                     if ~isnan(cfg_conditions(cn).perturbation)
%                         cond_trials = cond_trials & ...
%                             ([states_lfp(i).trials.perturbation]) == ...
%                             (cfg_conditions(cn).perturbation);
%                     end
                    if ~isnan(cfg_conditions(cn).block)
                        for b = cfg_conditions(cn).block
                            cond_trials = cond_trials | ...
                                ([states_lfp(i).trials.block] == b);
                        end
                    end
                    if ~isnan(cfg_conditions(cn).choice)
                        cond_trials = cond_trials & ...
                            ([states_lfp(i).trials.choice_trial] == ...
                            cfg_conditions(cn).choice);
                    end
                    cond_trials = cond_trials & ...
                        strcmp({states_lfp(i).trials.hndspc_lbl}, hs_labels(hs));
                    cond_based_evoked(cn).sites(i).ntrials(hs) = sum(cond_trials);
                                        
                    fprintf('Condition %s - %s\n', cfg_conditions(cn).label, hs_labels{hs});
                    fprintf('Total number of trials %g\n', sum(cond_trials));
                    
                    cond_based_evoked(cn).sites(i).noisytrials(hs) = ...
                        sum(cond_trials & [states_lfp(i).trials.noisy]); 
                                        
                    % consider only non noisy trials
                    fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                        & [states_lfp(i).trials.noisy]));
                    cond_trials = cond_trials & ~[states_lfp(i).trials.noisy];
                    
                                 
                    % loop through trials 

                    for st = 1:length(analyse_states)
                        state_tfs.time = {}; % timebins fo spectrogram
                        state_tfs.freq = {}; % freq bins
                        state_tfs.lfp = {}; % evoked LFP response
                        
                        for t = find(cond_trials)

                            states          = states_lfp(i).trials(t).states;
                            state_onset_t   = states([states(:).id] == ...
                                analyse_states{st}).onset_t;
                            state_start_t   = states([states(:).id] == ...
                                analyse_states{st}).start_t;
                            state_end_t     = states([states(:).id] == ...
                                analyse_states{st}).end_t;
                            
                            % evoked LFP for this state
                            state_tfs.lfp = [state_tfs.lfp, ...
                                states_lfp(i).trials(t).lfp_data(...
                                (states_lfp(i).trials(t).time >= state_start_t & ...
                                states_lfp(i).trials(t).time <= state_end_t))];
                            state_tfs.lfp_time = states_lfp(i).trials(t).time(...
                                (states_lfp(i).trials(t).time >= state_start_t & ...
                                states_lfp(i).trials(t).time <= state_end_t)) - state_onset_t;
                            
                        end

                        
                        if ~isempty(state_tfs.lfp)
                        
                            % evoked LFP average
                            % crop each lfp to same number of samples
                            nsamples = min(cellfun('length', state_tfs.lfp));
                            for k = 1:length(state_tfs.lfp)
                                state_tfs.lfp{k} = state_tfs.lfp{k}(1:nsamples);
                            end
                            state_tfs.lfp_time = state_tfs.lfp_time(1:nsamples);
                            arr_state_lfp = vertcat(state_tfs.lfp{:});
                            evoked_state_lfp_mean = nanmean(arr_state_lfp, 1);
                            evoked_state_lfp_std = nanstd(arr_state_lfp, 0, 1);

                            % save evoked LFP
                            cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).mean = evoked_state_lfp_mean;
                            cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).std = evoked_state_lfp_std; 
                            cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).time = state_tfs.lfp_time;
                            cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).trials = find(cond_trials);
                            cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).hs_label = hs_labels(hs);
                            cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).state = analyse_states{st};
                             
                            
                        end

                    end

                end
            else
                continue;

            end

        end
        
      
        % calculate the average evoked LFP across sites for this condition 
        if isfield(cond_based_evoked(cn).sites, 'site_evoked_lfp')
            % struct to store average evoked LFP across sites
            cond_based_evoked(cn).evoked_lfp_session = struct();%cell(length(analyse_states),length(hs_labels));
            
            
            for i = 1:length(cond_based_evoked(cn).sites)

                for hs = 1:length(hs_labels)
                    for st = 1:length(analyse_states)
                        if ~isempty(cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs))
                            
                            if i == 1
                                cond_based_evoked(cn).evoked_lfp_session(st, hs).mean = ...
                                    (1/length(cond_based_evoked(cn).sites))*cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).mean ;
                                cond_based_evoked(cn).evoked_lfp_session(st, hs).std = ...
                                    (1/length(cond_based_evoked(cn).sites))*cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).std ;
                                cond_based_evoked(cn).evoked_lfp_session(st, hs).time = ...
                                    cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).time;
                            else
                                nsamples = length(cond_based_evoked(cn).evoked_lfp_session(st, hs).time);
                                if nsamples > length(cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).time)
                                    nsamples = length(cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).time);
                                end
                                cond_based_evoked(cn).evoked_lfp_session(st, hs).mean = ...
                                    cond_based_evoked(cn).evoked_lfp_session(st, hs).mean(1:nsamples) + ...
                                    (1/length(cond_based_evoked(cn).sites)) * cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).mean(1:nsamples) ;
                                cond_based_evoked(cn).evoked_lfp_session(st, hs).std = ...
                                    cond_based_evoked(cn).evoked_lfp_session(st, hs).std(1:nsamples) + ...
                                    (1/length(cond_based_evoked(cn).sites)) * cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).std(1:nsamples) ;
                                
                            end
                            cond_based_evoked(cn).evoked_lfp_session(st, hs).hs_label = ...
                                cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).hs_label;
                            cond_based_evoked(cn).evoked_lfp_session(st, hs).state = ...
                                cond_based_evoked(cn).sites(i).site_evoked_lfp(st, hs).state;
                            cond_based_evoked(cn).evoked_lfp_session(st, hs).nsites = ...
                                length(cond_based_evoked(cn).sites);
                        end
                        
                    end
                end
            end
        end
               
    
        % plots

        % site averages
        for i = 1:length(cond_based_evoked(cn).sites)            
            % Evoked LFP
            if isfield(cond_based_evoked(cn).sites, 'site_evoked_lfp')
                plottitle = ['Site ID: ', cond_based_evoked(cn).sites(i).site_ID ', Target = ' ...
                    cond_based_evoked(cn).sites(i).target ', '  ...
                    '(block ' num2str(cfg_condition.block) '), '];
                if cfg_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                else
                    plottitle = [plottitle 'Choice trials'];
                end
                result_folder = fullfile(cfg_conditions(cn).results_folder, cond_based_evoked(cn).sites(i).site_ID);
                if ~exist(result_folder, 'dir')
                    mkdir(result_folder);
                end
                result_file = fullfile(cfg_conditions(cn).results_folder, cond_based_evoked(cn).sites(i).site_ID, ...
                    ['Evoked_LFP_' cond_based_evoked(cn).sites(i).site_ID '_' cfg_conditions(cn).label '.png']);
                
                lfp_tfa_plot_evoked_lfp (cond_based_evoked(cn).sites(i).site_evoked_lfp, lfp_tfa_cfg, ...
                    plottitle, result_file);
            end
        end        
        
        % plot average evoked LFP across sites for this session
        plottitle = ['Session: ', cond_based_evoked(cn).sites(1).session ', Target = ' cond_based_evoked(cn).sites(1).target ', '  ...
        'Block ' num2str(cfg_conditions(cn).block) ', '];
        if cfg_conditions(cn).choice == 0
            plottitle = [plottitle 'Instructed trials'];
        else
            plottitle = [plottitle 'Choice trials'];
        end
        result_file = fullfile(cfg_conditions(cn).results_folder, ['Evoked_LFP_' ...
            cond_based_evoked(cn).sites(1).session '_' cfg_conditions(cn).label '.png']);
        lfp_tfa_plot_evoked_lfp (cond_based_evoked(cn).evoked_lfp_session, lfp_tfa_cfg, ...
            plottitle, result_file);
        
    end
    
    close all;
    % save mat files
    save(fullfile(results_folder_tfr, 'cfg_conditions.mat'), 'cfg_conditions');
    save(fullfile(results_folder_tfr, 'cond_based_evoked.mat'), 'cond_based_evoked');
    save(fullfile(results_folder_tfr, 'lfp_tfa_cfg.mat'), 'lfp_tfa_cfg');
end
        
%         % read baseline tfs from all trials
%         baseline_pow = cell(1, length(states_lfp(i).trials));
%         for t = 1:length(states_lfp(i).trials)
%             trial = states_lfp(i).trials(t);
%             % trials to be considered for baseline calculation
%             consider_trial = trial.choice_trial == choice & trial.ipsilateral == hemisphere & ...
%                 trial.perturbation == 0 & trial.noisy == 0;
%             if consider_trial
% %             if t == 1
% %                 baseline_pow{t} = trial.tfs.powspctrm;
% %             else
%                 baseline_pow{t} = trial.tfs.powspctrm(1,:,...
%                 trial.tfs.time >= trial.baseline.start_t & ...
%                 trial.tfs.time <= trial.baseline.end_t);
%                 %baseline_pow{t} = tfs_baseline.powspctrm + baseline_pow;
%             end
%         end
%         
%         % find mean and std of baseline power for each frequency bin
%         concat_baseline_pow = cat(3, baseline_pow{:});
%         baseline_mean = nanmean(concat_baseline_pow, 3);
%         baseline_std  = nanstd(concat_baseline_pow, 0, 3);
%         
%         % handspace tuned LFP TF analysis
%         tfs_sites = struct();
%         % initialize 
%         for h = 1:length(hs_label)
%             tfs_sites(i).handspace(h) = struct();
%             tfs_sites(i).handspace(h).name = hs_label(h);
%             for s = 1:length(analyse_states)
%                 tfs_sites(i).handspace(h).state(s).name = analyse_states(s);
%                 tfs_sites(i).handspace(h).state(s).tfr_avg = struct();
%                 %tfs_sites(i).handspace(h).state(s).pow_norm = [];
%             end
%         end
%         
%         powspctrm = cell(length(hs_label), length(analyse_states));
%         %powspctrm_norm = length(states_lfp(i).trials);
%         for t = 1:length(states_lfp(i).trials)
%             trial = states_lfp(i).trials(t);
%             % trial to be considered for analysis
%             consider_trial = trial.choice_trial == choice & trial.ipsilateral == hemisphere & ...
%                 isany(trial.block == blocks) & trial.noisy == 0;
%             if consider_trial
%                 % get stored info about all states
%                 states = trial.states;
%                 % handspace tuning for the trial
%                 hndspc = trial.hndspc_lbl;
%                 hs_idx = find(hs_label == hndspc);
%                 for s = 1:length(analyse_states)
%                     st_idx = find(analyse_states(s) == states);
%                     % get the tfs for state and handspace tuning
%                     tfs_sites(i).handspace(hs_idx).state(s).tfr_avg = trial.tfs;
%                     powspctrm{t} = trial.tfs.powspctrm(1,:,...
%                         trial.tfs.time >= trial.states(st_idx).start_t & ...
%                         trial.tfs.time <= trial.states(st_idx).end_t);
%                     tfs_sites(i).handspace(hs_idx).state(s).tfr_avg.time = trial.tfs.time(1, ...
%                         trial.tfs.time >= trial.states(st_idx).start_t & ...
%                         trial.tfs.time <= trial.states(st_idx).end_t);
%                     tfs_sites(i).handspace(hs_idx).state(s).tfr_avg.time = ...
%                         tfs_sites(i).handspace(hs_idx).state(s).tfr_avg.time - ...
%                         trial.states(st_idx).onset_t;
%                     
%                     % baseline nomrlaization
%                     tfs_sites(i).handspace(hs_idx).state(s).tfr_avg.powspctrm = baselineNormalization(...
%                         tfs_sites.handspace(hs_idx).state(s).pow, baseline_mean, baseline_std, 'division');
%                     % combined TFR for all states
%                     tfs_sites(i).handspace(hs_idx).tfr_avg = ...
%                         tfs_sites(i).handspace(hs_idx).state(s).tfr_avg;
%                     
%                 end
%             end
%             
%         end
%         
%         % save tfs for the site
%         states_lfp(i).tfs = tfs_sites(i);
%         
%         % plot TFR for each site and calculate average across sites
%         figure; nsubplot = 0;
%         for h = 1:length(tfs_sites(i).handspace)
%             for s = 1:length(tfs_sites(i).handspace(h).state)
%                 tfr_avg = tfs_sites(i).handspace(h).state(s).tfr_avg;
%                 if ~isempty(tfr_avg)
%                     % configuration for ft_singleplotTFR
%                     cfg = [];
%                     cfg.baseline     = 'no'; %baseline_shift;                 % -400ms to -100ms before the onset of first state
%                     cfg.baselinetype = 'absolute';  
%                     cfg.maskstyle    = 'saturation';
%                     %cfg.interactive  = 'no';
%                     cfg.zlim         = [-2 2];
%                     %cfg.masknans     = 'yes';
%                     cfg.channel  = concat_states_TFR.label;
%                     nsubplot = nsubplot + 1;
%                     subplot(2,4,nsubplot)
%                     ft_singleplotTFR(cfg, concat_states_TFR);
%                     line([0 0], ylim);
%                     title(tfs_sites(i).handspace(h).name + ': ' + ...
%                         tfs_sites(i).handspace(h).state(s).name);
%                 end
%             end
%         end
%     end
%     % calculate the average TFR across sites
%     tfs_avg_sites = struct();
%     % initialize
%     tfs_avg_sites = tfs_sites(1);
%     tfs_powspctrm = cat(3, tfs_sites(:).handspace(h).state(s).tfr_avg.powspctrm);
%     for i = 1:length(tfs_sites)
%         for h = 1:length(tfs_sites(i).handspace)
%             for s = 1:length(tfs_sites(i).handspace(h).state)
%                 tfr_avg = tfs_sites(i).handspace(h).state(s).tfr_avg;
% end
%         
%         
%         
%         
%                     
%         
% %         % baseline for this site
% %         data_lfp_baseline = ft_data_sites(i).baseline;
% %         % get baseline power for normalization
% %         if ~isempty(data_lfp_baseline.trial)
% %             nsample = length(data_lfp_baseline.time{1});
% %             fs = data_lfp_baseline.fsample;
% %             ts = 1/fs;
% %             tstart = data_lfp_baseline.time{1}(1);
% %             tend = data_lfp_baseline.time{1}(end);
% %             % configuration for TFR of baseline
% %             cfg              = [];
% %             cfg.output       = 'pow';
% %             cfg.method       = 'mtmconvol';
% %             cfg.taper        = 'hanning';
% %             cfg.pad          = 'nextpow2';
% %             cfg.foi          = fs/500:fs/500:fs/10;             % analysis 2 to 100 Hz in steps of 2 Hz
% %             cfg.t_ftimwin    = ones(length(cfg.foi),1).*200*ts;    % length of time window = 0.2 sec
% %             %cfg.t_ftimwin    = 4./cfg.foi;    % 7 cycles per time window
% %             cfg.toi          = linspace(tstart, tend, nsample);% time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
% %             cfg.channel      = data_lfp_baseline.label;
% %             disp(['Processing ' cfg.channel]);
% %             TFR_baseline     = ft_freqanalysis(cfg, data_lfp_baseline);
% %             % now save the baseline TFR to corresponding hand space tuning
% %             ft_data_sites(i).TFR_baseline = TFR_baseline;
% %             % get mean value of baseline power for each frequency
% %             baseline_mean = nanmean(TFR_baseline.powspctrm, 3);
% %             % get std value of baseline power for each frequency
% %             baseline_std = nanstd(TFR_baseline.powspctrm, 0, 3);
% %         else
% %             ft_data_sites(i).TFR_baseline = [];
% %         end
% %         
% %         % loop through each hand-space tuning
% %         for hs = 1:length(hs_label)            
% %             % loop through each state
% %             for st = 1:length(ft_data_sites(i).states)
% %                 % now this is a FT_DATATYPE_RAW on which FT time frequency
% %                 % analzsis can be done
% %                 data_lfp_st_hs = ft_data_sites(i).states(st).hndspc(hs);
% %                 nsample = length(ft_data_sites(i).states(st).hndspc(hs).time{1});
% %                 if ~isempty(data_lfp_st_hs.trial)
% %                     fs = data_lfp_st_hs.fsample;
% %                     ts = 1/fs;
% %                     tstart = data_lfp_st_hs.time{1}(1);
% %                     tend = data_lfp_st_hs.time{1}(end);
% %                     % configuration for TFR
% %                     cfg              = [];
% %                     cfg.output       = 'pow';
% %                     cfg.method       = 'mtmconvol';
% %                     cfg.taper        = 'hanning';
% %                     cfg.pad          = 'nextpow2';
% %                     cfg.foi          = fs/500:fs/500:fs/10;             % analysis 2 to 100 Hz in steps of 2 Hz
% %                     cfg.t_ftimwin    = ones(length(cfg.foi),1).*200*ts;    % length of time window = 0.2 sec
% %                     % 7 cycles per time window
% %                     cfg.toi          = linspace(tstart, tend, nsample);% time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
% %                     cfg.channel      = data_lfp_st_hs.label;
% %                     disp(['Processing ' cfg.channel]);
% %                     TFRhann_fixed = ft_freqanalysis(cfg, data_lfp_st_hs); 
% %                     % now save the TFR to corresponding hand space tuning
% %                     ft_data_sites(i).TFR{st, hs} = TFRhann_fixed;
% %                     % baseline normalization
% % %                     TFRhann_fixed_norm = baselinePowerNormalize(TFRhann_fixed, ...
% % %                         TFR_baseline, 'subtraction');
% %                     % subtract baseline mean for each frequency
% %                     TFRhann_fixed_norm = TFRhann_fixed;
% %                     TFRhann_fixed_norm.powspctrm = (TFRhann_fixed.powspctrm - ...
% %                         repmat(baseline_mean, [1 1 size(TFRhann_fixed.powspctrm, 3)])) ./ ...
% %                         repmat(baseline_std, [1 1 size(TFRhann_fixed.powspctrm, 3)]);
% %                     % concatenate TFRs for each state
% %                     if st == 1
% %                         concat_states_TFR = TFRhann_fixed_norm;
% %                     else
% %                         % concatenate times
% %                         concat_states_TFR.time = [concat_states_TFR.time, TFRhann_fixed_norm.time];
% %                         concat_states_TFR.powspctrm = cat(3, concat_states_TFR.powspctrm, TFRhann_fixed_norm.powspctrm);
% %                     end
% % 
% %                     % accumulate TFR per site for averaging across sites
% %                     if i == 1
% %                         TFR_avg{st, hs} = TFRhann_fixed_norm;
% %                         TFR_avg{st, hs}.powspctrm = (1/nsites) * TFRhann_fixed_norm.powspctrm;
% %                         TFR_avg{st, hs}.ntrials = length(data_lfp_st_hs.trial);
% %                     else
% %                         TFR_avg{st, hs}.powspctrm = TFR_avg{st, hs}.powspctrm + (1/nsites) * TFRhann_fixed_norm.powspctrm;
% %                         TFR_avg{st, hs}.ntrials = TFR_avg{st, hs}.ntrials + length(data_lfp_st_hs.trial);
% %                     end
% %                     TFR_avg{st, hs}.label = strcat(analyse_states(st),hs_label(hs));
% %                 else
% %                     ft_data_sites(i).TFR{st, hs} = [];
% %                 end
% %             
% %             end 
% % 
% %             % Visualization
% %             if ~isempty(concat_states_TFR) > 0
% %                 % normalize power spectrogram by z-score
% % %                 concat_states_TFR.powspctrm = (concat_states_TFR.powspctrm - ...
% % %                     nanmean(concat_states_TFR.powspctrm(:))) / nanstd(concat_states_TFR.powspctrm(:));
% %                 time_to_event = concat_states_TFR.time;
% %                 tshift = time_to_event(1);
% %                 % first rearrange the time axis of TFR
% %                 bin_t = concat_states_TFR.time(2) - concat_states_TFR.time(1);
% %                 concat_states_TFR.time = 0:bin_t:bin_t*(length(concat_states_TFR.time)-1);
% %                 % now find the state onsets in the new time axis
% %                 rearr_state_onset_t = concat_states_TFR.time(time_to_event == 0);
% %                 %600*ts%0:1:length(concat_states_TFR.time)-1;
% %                 % rearrange the baseline
% %                 baseline = [-0.4 -0.1];
% %                 baseline_shift = baseline - tshift;
% % 
% %                 cfg = [];
% %                 cfg.baseline     = 'no'; %baseline_shift;                 % -400ms to -100ms before the onset of first state
% %                 cfg.baselinetype = 'absolute';  
% %                 cfg.maskstyle    = 'saturation';
% %                 %cfg.interactive  = 'no';
% %                 cfg.zlim         = [-2 2];
% %                 %cfg.masknans     = 'yes';
% %                 cfg.channel  = concat_states_TFR.label;
% %                 subplot(2,2,hs)
% %                 ft_singleplotTFR(cfg, concat_states_TFR);
% %                 % mark state onsets
% %                 set(gca,'xtick',rearr_state_onset_t)
% %                 set(gca,'xticklabels',strrep(analyse_states, '_', '\_'))
% %                 for t = rearr_state_onset_t
% %                     line([t t], ylim);
% %                 end
% %                 xlabel('Time');
% %                 ylabel('Frequency');
% %                 title(sprintf('%s (ntrials = %g)', cfg.channel{1}, length(data_lfp_st_hs.trial)));
% %                 line([0 0], ylim, 'color', 'k');
% %             end
% %         end
% %         % plot title
% %     %     sgtitle(sprintf('Session: %s, Target = %s, Site: %s (nblocks = %g)', ft_data_sites(i).session, ...
% %     %         ft_data_sites(i).target, ft_data_sites(i).siteID, ft_data_sites(i).nblocks));
% %         plottitle = ['Session: ' ft_data_sites(i).session ', Target = ' ft_data_sites(i).target ', Site ID: ' ft_data_sites(i).siteID ...
% %             '(nblocks = ' num2str(ft_data_sites(i).nblocks) ')'];
% %         ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
% %             , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% %         saveas(gca, fullfile(fig_folder_tfr, [ft_data_sites(i).siteID '.png']));
% %     end
% % 
% %     % save variables
% %     save ft_data_sites ft_data_sites;
% %     %save TFR_avg TFR_avg;
% % 
% %     %%% Visualiyation of time frequency average across site %%%
% %     figure; colormap jet;
% %     for hs = 1:length(hs_label)
% %         for st = 1:length(analyse_states)
% %             % first concatenate TFR for different states for one hs tuning for
% %             % visualization
% %             if st == 1
% %                 concat_TFR_avg = TFR_avg{st,hs};
% %                 concat_TFR_avg.label = hs_label{hs};
% %             else
% %                 concat_TFR_avg.time = [concat_TFR_avg.time TFR_avg{st,hs}.time];
% %                 concat_TFR_avg.powspctrm = cat(3, concat_TFR_avg.powspctrm, TFR_avg{st,hs}.powspctrm);
% %             end
% % 
% %         end
% % 
% %         % Visualization
% % 
% %         time_to_event = concat_TFR_avg.time;
% %         tshift = time_to_event(1);
% %         % first rearrange the time axis of TFR
% %         bin_t = concat_TFR_avg.time(2) - concat_TFR_avg.time(1);
% %         concat_TFR_avg.time = 0:bin_t:bin_t*(length(concat_TFR_avg.time)-1);
% %         % now find the state onsets in the new time axis
% %         rearr_state_onset_t = concat_TFR_avg.time(time_to_event == 0);
% %         % rearrange the baseline
% %         baseline = [-0.4 -0.1];
% %         baseline_shift = baseline - tshift;
% % 
% %         cfg = [];
% %         cfg.baseline     = 'no';%baseline_shift;                 % -400ms to -100ms before the state onset
% %         cfg.baselinetype = 'absolute';  
% %         cfg.maskstyle    = 'saturation';
% %         cfg.zlim         = [-2 2];
% %         cfg.masknans     = 'yes';
% %         %concat_TFR_avg.label = TFR_avg{st,hs}.label{1};
% %         cfg.channel  = concat_TFR_avg.label;
% %         subplot(2,2,hs)
% %         ft_singleplotTFR(cfg, concat_TFR_avg);
% %         % mark state onsets
% %         set(gca,'xtick',rearr_state_onset_t)
% %         set(gca,'xticklabels',strrep(analyse_states, '_', '\_'))
% %         for t = rearr_state_onset_t
% %             line([t t], ylim);
% %         end
% %         xlabel('Time');
% %         ylabel('Frequency');
% %         title(sprintf('%s', cfg.channel{1}));
% %         line([0 0], ylim, 'color', 'k');
% %     end
% % 
% %     % plot title
% %     plottitle = ['Session: ' ft_data_sites(i).session '(nsites = ', ...
% %         num2str(length(ft_data_sites)), 'nblocks = ' num2str(ft_data_sites(i).nblocks) ')'];
% %     ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
% %         , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
% %     % sgtitle(sprintf('Session: %s, Target = %s, (nsites = %g)', ft_data_sites(1).session, ...
% %     %     ft_data_sites(1).target(1:end-2), length(ft_data_sites)));
% %     saveas(gca, fullfile(fig_folder_tfr, [ft_data_sites(1).session '.png']));
% % 
% % end
% % 
