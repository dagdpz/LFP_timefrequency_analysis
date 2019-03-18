function [ cond_based_tfs ] = lfp_tfa_compute_plot_tfr( states_lfp, analyse_states, cfg_condition, cfg_baseline, root_results_folder ) 

% lfp_tfa_compute_plot_tfr  - plots lfp time freq response for
% different hand-space tuning conditions for each site and across sites
%
% USAGE:
%	[ cond_based_tfs ] = computeAndPlotTFR( states_lfp, analyse_states, cfg_condition, cfg_baseline, root_results_folder )
%
% INPUTS:
%		states_lfp  	- struct containing lfp data per trial, output from
%		lfp_tfa_compute_baseline
%       analyse_states  - ids of states to be analysed, see 
%       lfp_tfa_define_states
%		cfg_condition   - configuration structure to specify the conditions
%		for selecting trials for analysis
%           choice          : whether instructed (0) or choice (1) trials to be
%           selected
%           perturbation    : whether control or inactivation trials to be
%           selected
%           blocks          : blocks to be selected
%           recorded_hemispace : hemispace from which LFP is recorded that
%           has to be analysed
%       cfg_baseline    - configuration structure for baseline corection
%           method          : method to be used for baseline correction
%           ('subtraction', 'division', 'relchange', 'zscore')
%       root_results_folder     - folder to save results
% 
% OUTPUTS:
%		cond_based_tfs	- output structure which saves the average tfs for  
%                         trials of a given condition for different handspace 
%                         tunings and periods around the states analysed
%                         same datastructure as input ft_data_sites, but
%                         with additional field to store condition-wise lfp
%                         time freq response		
%       TFR_avg         - condition-wise lfp time freq response across all
%                         sites
%
%
% See also lfp_tfa_compute_baseline, lfp_tfa_define_states
    
    
    % make a folder to save figures
    sessionName = states_lfp(1).session;
    results_folder_tfr = fullfile(root_results_folder, date, sessionName, 'Condition_based_TFS');
    if ~exist(results_folder_tfr, 'dir')
        mkdir(results_folder_tfr);
    end
    
    % create a label for the given condition
    cond_label = [];
    if cfg_condition.recorded_hemispace == 'L'
        cond_label = [cond_label 'Left_hemispace_'];
    else
        cond_label = [cond_label 'Right hemispace_'];
    end
    if cfg_condition.choice == 0
        cond_label = [cond_label 'Instructed_'];
    else
        cond_label = [cond_label 'Choice_'];
    end
    if cfg_condition.perturbation == 0
        cond_label = [cond_label 'Control_'];
    else
        cond_label = [cond_label 'Inactivation_'];
    end
    cond_label = [cond_label, 'Block_', num2str(cfg_condition.blocks)];
    
    % create a folder for storing results for this condition
    cond_tfs_folder = fullfile(results_folder_tfr, cond_label);
    if ~exist(cond_tfs_folder, 'dir')
        mkdir(cond_tfs_folder)
    end
    
    % condition based TFS
    cond_based_tfs = struct();
    cond_based_tfs.cfg = cfg_condition;

    % states to be analysed
%     analyse_states = {states_lfp(1).states.name};

    % hand-space tuning of LFP
    hs_labels = unique({states_lfp(1).trials.hndspc_lbl});
    
    % cell array to store time frequency average across sites
    TFR_avg = cell(length(analyse_states),length(hs_labels));
    
    % loop through each site
    nsites = length(states_lfp);
    
    % loop through each site
    for i = 1:length(states_lfp)
        
        % consider site based on recorded hemispace
        if strcmp(states_lfp(i).recorded_hemispace, cfg_condition.recorded_hemispace) 
            cond_based_tfs(i).site_ID = states_lfp(i).site_ID;
            % make a struct for concatenating TFR for all states
            cond_based_tfs(i).tfs_avg_site = struct(); 
            %cond_based_tfs(i).tfs_avg_site.powspctrm = cell(length(analyse_states),length(hs_labels));
            cond_based_tfs(i).ntrials = zeros(1,length(hs_labels));
            
            for hs = 1:length(hs_labels)
                cond_trials = ones(1, length(states_lfp(i).trials));
                % consider only non noisy trials
                cond_trials = cond_trials & ~[states_lfp(i).trials.noisy];
                % get the trials for given condition and this hs label
                if ~isnan(cfg_condition.perturbation)
                    cond_trials = cond_trials & [states_lfp(i).trials.perturbation] == cfg_condition.perturbation;
                end
                if ~isnan(cfg_condition.blocks)
                    cond_trials = cond_trials & ([states_lfp(i).trials.block] == cfg_condition.blocks);
                end
                if ~isnan(cfg_condition.choice)
                    cond_trials = cond_trials & ([states_lfp(i).trials.choice_trial] == cfg_condition.choice);
                end
                cond_trials = cond_trials & strcmp({states_lfp(i).trials.hndspc_lbl}, hs_labels(hs));
                cond_based_tfs(i).ntrials(hs) = sum(cond_trials);

                % loop through trials 

                for st = 1:length(analyse_states)
                    state_tfs.powscptrm = {};
                    state_tfs.time = {};
                    state_tfs.freq = {};
                    state_tfs.lfp = {};

                    for t = find(cond_trials)

                        states          = states_lfp(i).trials(t).states;
                        state_onset_t   = states([states(:).id] == analyse_states{st}).onset_t;
                        state_start_t   = states([states(:).id] == analyse_states{st}).start_t;
                        state_end_t     = states([states(:).id] == analyse_states{st}).end_t;

                        % crop the tfs for this state
                        state_tfs.powscptrm = [state_tfs.powscptrm, ...
                            states_lfp(i).trials(t).tfs.powspctrm(1, :, ...
                            (states_lfp(i).trials(t).tfs.time >= state_start_t & ...
                            states_lfp(i).trials(t).tfs.time <= state_end_t))];
    %                     if sum (states_lfp(i).trials(t).tfs.time >= state_start_t & ...
    %                         states_lfp(i).trials(t).tfs.time <= state_end_t) < ntimebins
                        state_tfs.time = [state_tfs.time, states_lfp(i).trials(t).tfs.time(1, ...
                            (states_lfp(i).trials(t).tfs.time >= state_start_t & ...
                            states_lfp(i).trials(t).tfs.time <= state_end_t)) - state_onset_t];
                        state_tfs.freq = [state_tfs.freq states_lfp(i).trials(t).tfs.freq];  
    %                     end
                        % evoked LFP for this state
    %                     state_tfs.lfp = [state_tfs.lfp, states_lfp(i).trials(t).lfp_data(1, ...
    %                         (states_lfp(i).trials(t).time >= state_start_t & ...
    %                         states_lfp(i).trials(t).time <= state_end_t))];
                        % LFP power spectrum
    %                     state_tfs.raw_psd =                     

                    end

                    % find number of time bins in power
                    % spectrogram
                    ntimebins = min(cellfun('length', state_tfs.time));
                    nfreqbins = numel(state_tfs.freq);
                    % crop each tfs to the ntimebins
                    for k = 1:length(state_tfs.powscptrm)
                        state_tfs.powscptrm{k} = state_tfs.powscptrm{k}(1,:,1:ntimebins);
                    end
                    % average power spectrum for each state
                    arr_state_pow = zeros(1, nfreqbins, ntimebins);

                    % find the average TFS for each state
                    arr_state_pow = cat(1, state_tfs.powscptrm{:});
                    state_tfs.powspctrm_rawmean = nanmean(arr_state_pow, 1);
                    % baseline normalization
                    cfg_baseline.mean = states_lfp(i).baseline_mean;
                    cfg_baseline.std = states_lfp(i).baseline_std;
                    state_tfs.powspctrm_normmean = lfp_tfa_baseline_normalization(...
                        state_tfs.powspctrm_rawmean, cfg_baseline); 
                    % save average tfs
                    cond_based_tfs(i).tfs_avg_site(st, hs).powspctrm = state_tfs.powspctrm_normmean;
                    cond_based_tfs(i).tfs_avg_site(st, hs).time = state_tfs.time{1};
                    cond_based_tfs(i).tfs_avg_site(st, hs).freq = state_tfs.freq{1};                
                    cond_based_tfs(i).tfs_avg_site(st, hs).label = hs_labels(hs);

                end
                
            end
        
        end
                        
    end
    
    % calculate the average across sites for this condition 
    if isfield(cond_based_tfs, 'tfs_avg_site')
        % struct to store average tfr across sites
        cond_based_tfs(1).tfs_avg_session = cell(length(analyse_states),length(hs_labels));
        
        for i = 1:length(cond_based_tfs)
            
            for hs = 1:length(hs_labels)
                for st = 1:length(analyse_states)
                    if i == 1
                        cond_based_tfs(1).tfs_avg_session{st, hs}.powspctrm = ...
                            (1/length(cond_based_tfs))*cond_based_tfs(i).tfs_avg_site(st, hs).powspctrm ;

                    else
                        cond_based_tfs(1).tfs_avg_session{st, hs}.powspctrm = ...
                            cond_based_tfs(1).tfs_avg_session{st, hs}.powspctrm + ...
                            (1/length(cond_based_tfs)) * cond_based_tfs(i).tfs_avg_site(st, hs).powspctrm ;
                    end
                    cond_based_tfs(1).tfs_avg_session{st, hs}.time = cond_based_tfs(i).tfs_avg_site(st, hs).time;
                    cond_based_tfs(1).tfs_avg_session{st, hs}.freq = cond_based_tfs(i).tfs_avg_site(st, hs).freq;
                    cond_based_tfs(1).tfs_avg_session{st, hs}.label = cond_based_tfs(i).tfs_avg_site(st, hs).label;                

                end
            end
        end
    end
    
    % plots
    % site averages
    for i = 1:length(cond_based_tfs)
        if isfield(cond_based_tfs, 'tfs_avg_site')
            figure;
            cm = colormap('jet');
            cm(1,:,:) = [1,1,1];
            colormap(cm);
            
            % loop through handspace
            for hs = 1:size(cond_based_tfs(i).tfs_avg_site, 2)
                % concatenate states
                concat_states_tfs = struct();
                concat_states_tfs.powspctrm = [];
                concat_states_tfs.state_time = [];
                concat_states_tfs.freq = cond_based_tfs(i).tfs_avg_site(1, hs).freq;
                concat_states_tfs.label = cond_based_tfs(i).tfs_avg_site(1, hs).label;

                state_info = struct();
                for st = 1:size(cond_based_tfs(i).tfs_avg_site, 1)
                    
                    
                    concat_states_tfs.powspctrm = cat(3, concat_states_tfs.powspctrm, cond_based_tfs(i).tfs_avg_site(st, hs).powspctrm);
                    concat_states_tfs.state_time = [concat_states_tfs.state_time, ...
                        cond_based_tfs(i).tfs_avg_site(st, hs).time];
                    % append nans to separate the states
                    concat_states_tfs.powspctrm = cat(3, concat_states_tfs.powspctrm, ...
                        nan(1, length(concat_states_tfs.freq), 100));
                    concat_states_tfs.state_time = [concat_states_tfs.state_time, ...
                        nan(1, 100)];
                    
                    state_info(st).onset_s = find(cond_based_tfs(i).tfs_avg_site(st, hs).time <= 0, 1, 'last');
                    state_info(st).onset_t = 0;
                    state_info(st).start_s = 1;
                    state_info(st).start_t = cond_based_tfs(i).tfs_avg_site(st, hs).time(1);
                    state_info(st).finish_s = length(cond_based_tfs(i).tfs_avg_site(st, hs).time);
                    state_info(st).finish_t = cond_based_tfs(i).tfs_avg_site(st, hs).time(end);                    
                    
                    if st > 1
                        state_info(st).start_s = length(cond_based_tfs(i).tfs_avg_site(st-1, hs).time) + ...
                            state_info(st).start_s + (st-1)*100;
                        state_info(st).finish_s = length(cond_based_tfs(i).tfs_avg_site(st-1, hs).time) + ...
                            state_info(st).finish_s + (st-1)*100;
                        state_info(st).onset_s = length(cond_based_tfs(i).tfs_avg_site(st-1, hs).time) + ...
                            state_info(st).onset_s + (st-1)*100;
                    end
                    
                end
                concat_states_tfs.time = 1:1:size(concat_states_tfs.powspctrm, 3);
                state_onsets = find(concat_states_tfs.state_time(1:end-1) .* ...
                    concat_states_tfs.state_time(2:end) <= 0);
                state_samples = sort([state_info.start_s, state_info.onset_s, ...
                    state_info.finish_s]);

                % now plot
                subplot(2,2,hs)
                cfg = [];
                cfg.baseline     = 'no'; %baseline_shift;                 % -400ms to -100ms before the onset of first state
                cfg.maskstyle    = 'saturation';
                cfg.interactive  = 'no';
                cfg.zlim         = [-1 1];
                cfg.channel  = concat_states_tfs.label;
                subplot(2,2,hs)
                ft_singleplotTFR(cfg, concat_states_tfs);
                % mark state onsets
                set(gca,'xtick',state_samples)
%                 xticklabels = [];
                for so = state_onsets
                    line([so so], ylim); 
                end
                set(gca,'xticklabels', round(concat_states_tfs.state_time(state_samples), 2))
                xlabel('Time');
                ylabel('Frequency');
                title(sprintf('%s (ntrials = %g)', cfg.channel{1}, cond_based_tfs(i).ntrials(hs)));
                line([0 0], ylim, 'color', 'k');
            end
            plottitle = ['Site ID: ', states_lfp(i).site_ID ', Target = ' states_lfp(i).target ', '  ...
            '(block ' num2str(cfg_condition.blocks) '), '];
            if cfg_condition.choice == 0
                plottitle = [plottitle 'Instructed trials'];
            else
                plottitle = [plottitle 'Choice trials'];
            end
            ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
                , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
            saveas(gca, fullfile(cond_tfs_folder, [states_lfp(i).site_ID '.png']));
        end
    end
    % plot average across sites for this session
    figure;
    cm = colormap('jet');
    cm(1,:,:) = [1,1,1];
    colormap(cm);
    % loop through handspace
    for hs = 1:size(cond_based_tfs(i).tfs_avg_site, 2)
        % concatenate states
        concat_states_tfs = struct();
        concat_states_tfs.powspctrm = [];
        concat_states_tfs.state_time = [];
        concat_states_tfs.freq = cond_based_tfs(i).tfs_avg_site(1, hs).freq;
        concat_states_tfs.label = cond_based_tfs(i).tfs_avg_site(1, hs).label;

        state_info = struct();
        for st = 1:size(cond_based_tfs(i).tfs_avg_site, 1)


            concat_states_tfs.powspctrm = cat(3, concat_states_tfs.powspctrm, cond_based_tfs(i).tfs_avg_site(st, hs).powspctrm);
            concat_states_tfs.state_time = [concat_states_tfs.state_time, ...
                cond_based_tfs(i).tfs_avg_site(st, hs).time];
            % append nans to separate the states
            concat_states_tfs.powspctrm = cat(3, concat_states_tfs.powspctrm, ...
                nan(1, length(concat_states_tfs.freq), 100));
            concat_states_tfs.state_time = [concat_states_tfs.state_time, ...
                nan(1, 100)];

            state_info(st).onset_s = find(cond_based_tfs(i).tfs_avg_site(st, hs).time <= 0, 1, 'last');
            state_info(st).onset_t = 0;
            state_info(st).start_s = 1;
            state_info(st).start_t = cond_based_tfs(i).tfs_avg_site(st, hs).time(1);
            state_info(st).finish_s = length(cond_based_tfs(i).tfs_avg_site(st, hs).time);
            state_info(st).finish_t = cond_based_tfs(i).tfs_avg_site(st, hs).time(end);                    

            if st > 1
                state_info(st).start_s = length(cond_based_tfs(i).tfs_avg_site(st-1, hs).time) + ...
                    state_info(st).start_s + (st-1)*100;
                state_info(st).finish_s = length(cond_based_tfs(i).tfs_avg_site(st-1, hs).time) + ...
                    state_info(st).finish_s + (st-1)*100;
                state_info(st).onset_s = length(cond_based_tfs(i).tfs_avg_site(st-1, hs).time) + ...
                    state_info(st).onset_s + (st-1)*100;
            end

        end
        concat_states_tfs.time = 1:1:size(concat_states_tfs.powspctrm, 3);
        state_onsets = find(concat_states_tfs.state_time(1:end-1) .* ...
            concat_states_tfs.state_time(2:end) <= 0);
        state_samples = sort([state_info.start_s, state_info.onset_s, ...
            state_info.finish_s]);

        % now plot
        subplot(2,2,hs)
        cfg = [];
        cfg.baseline     = 'no'; %baseline_shift;                 % -400ms to -100ms before the onset of first state
        cfg.maskstyle    = 'saturation';
        cfg.interactive  = 'no';
        cfg.zlim         = [-1 1];
        cfg.channel  = concat_states_tfs.label;
        subplot(2,2,hs)
        ft_singleplotTFR(cfg, concat_states_tfs);
        % mark state onsets
        set(gca,'xtick',state_samples)
%                 xticklabels = [];
        for so = state_onsets
            line([so so], ylim); 
        end
        set(gca,'xticklabels', round(concat_states_tfs.state_time(state_samples), 2))
        xlabel('Time');
        ylabel('Frequency');
        title(sprintf('%s (ntrials = %g)', cfg.channel{1}, cond_based_tfs(i).ntrials(hs)));
        line([0 0], ylim, 'color', 'k');
    end
    plottitle = ['Session: ', states_lfp(i).session ', Target = ' states_lfp(i).target ', '  ...
    'Block ' num2str(cfg_condition.blocks) ', '];
    if cfg_condition.choice == 0
        plottitle = [plottitle 'Instructed trials'];
    else
        plottitle = [plottitle 'Choice trials'];
    end
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    saveas(gca, fullfile(cond_tfs_folder, [states_lfp(i).session cond_label '.png']));
    
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
