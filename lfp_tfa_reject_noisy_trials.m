function [filt_session_lfp, noisy_lfp_trials] = lfp_tfa_reject_noisy_trials( states_lfp, cfg_noise )
%lfp_tfa_reject_noisy_trials  - marks noisy lfp trials based on threshold 
%values of raw lfp amplitude and raw lfp derivative
%
% USAGE:
%	[sites, noisy_lfp_trials] = lfp_tfa_reject_noisy_trials( sites, cfg_nr )
%
% INPUTS:
%		states_lfp       - struct containing processed lfp data which is 
%                           created by lfp_tfa_read_LFP()
%		cfg_noise        - configuration for noise rejection
%           methods      : methods to be used for noise rejection ('raw',
%           'diff', 'std', 'pow')
%           amp_thr      : threshold for lfp raw amplitude (no of stdev)
%           amp_N        : n consecutive samples beyond the raw amplitude 
%           threshold to be considered as noisy
%           std_thr      : threshold for trial-wise stdev compared to stdev
%           of all trials
%           diff_thr     : threshold for lfp derivative (no of stdevs)
%           diff_N       : n consecutive samples beyond the lfp derivative 
%           threshold to be considered as noisy
%           pow_thr      : threshold for lfp spectral power (no of stdevs)
%           results_folder: root folder to save output plots and results
%
% OUTPUTS:
%		filt_session_lfp - struct containing processed lfp with
%                          additional field which marks noisy trials
%       noisy_lfp_trials - struct with detailed information about noisy 
%                          trials
%
% See also lfp_tfa_read_LFP
    
    sessionName = states_lfp(1).session;
    results_folder = fullfile(cfg_noise.results_folder, sessionName, date);
    results_folder_noise = fullfile(cfg_noise.results_folder, date, sessionName, 'LFP_Noise_Rejection');
    if ~exist(results_folder_noise, 'dir')
        mkdir(results_folder_noise);
    end
    
    fprintf('=============================================================\n');
    fprintf('Noise Rejection in LFP\n');
    fprintf('=============================================================\n');
    
    
    % struct to save info about noisz trials
    noisy_lfp_trials = struct();
    noisy_lfp_trials.cfg = cfg_noise;

    % first, calculate the mean and std of length-wise concatenated LFP from
    % each site
    % read in lfp data for each site
    for i = 1:length(states_lfp)
        fprintf('Processing site %s ...\n', states_lfp(i).site_ID);
        site = states_lfp(i);
        % cell arrays to store LFP data for all trials
        concat_site_lfp = cell(1,length(site.trials)); % store LFP raw data for each trial
        concat_site_time = cell(1, length(site.trials)); % store LFP timestamps for each trial
        concat_site_diff_lfp = cell(1,length(site.trials)); % store derivative of LFP for each trial
        concat_blocks = zeros(1, length(site.trials)); % store blocks for each trial
        concat_site_lfp_pow = cell(1,length(site.trials)); % store LFP spectral power for each trial
        concat_site_states = cell(1,length(site.trials)); % store state ids for each trial
        concat_site_states_onset = cell(1,length(site.trials)); % store state onset times for each trial
        
        for t = 1:length(site.trials)
            % read lfp data for each trial
            trial = site.trials(t);
            block = trial.block;
            lfp_data = trial.lfp_data;
            lfp_tfs = trial.tfs;
            lfp_time = trial.time;
            % create a new field to indicate if the trial is noisy
            states_lfp(i).trials(t).noisy = 0;
            
            trial_start_state = trial.states(1);
            trial_end_state = trial.states(end);
            trial_time = lfp_time(lfp_time >= trial.trialperiod(1) & ...
                lfp_time <= trial.trialperiod(2));
            trial_lfp_raw = lfp_data(lfp_time >= trial.trialperiod(1) & ...
                lfp_time <= trial.trialperiod(2));
            trial_lfp_pow = lfp_tfs.powspctrm(1,:,...
                lfp_tfs.time >= trial.trialperiod(1) & ...
                lfp_tfs.time <= trial.trialperiod(2));
            
            fs = trial.fsample;
            ts = 1/fs;

            % save lfp for this trial into cell array
            concat_site_lfp{t} = trial_lfp_raw;
            concat_site_diff_lfp{t} = [nan diff(trial_lfp_raw)];
            concat_site_time{t} = trial_time;
            concat_blocks(t) = block;            
            concat_site_lfp_pow{t} = trial_lfp_pow;
            concat_site_states{t} = [trial.states.id];
            concat_site_states_onset{t} = [trial.states.onset_t];

        end
        % concatenate all datapoints to a 1d array
        arr_concat_site_lfp = horzcat(concat_site_lfp{:});
        % now get the mean and std of LFP for each site
        site_lfp_mean = mean(arr_concat_site_lfp);
        site_lfp_std = std(arr_concat_site_lfp);
        % get lfp raw threshold
        lfp_raw_minbound = site_lfp_mean - cfg_noise.amp_thr * site_lfp_std;
        lfp_raw_maxbound = site_lfp_mean + cfg_noise.amp_thr * site_lfp_std;
        % save this data into the struct
        noisy_lfp_trials(i).lfp_mean = site_lfp_mean;
        noisy_lfp_trials(i).lfp_std = site_lfp_std;  
        
        % concatenate all derivatives to a 1d array
        arr_concat_diff_lfp = horzcat(concat_site_diff_lfp{:});
        % now get the mean and std of LFP for each site
        lfp_diff_mean = nanmean(arr_concat_diff_lfp);
        lfp_diff_std = nanstd(arr_concat_diff_lfp);
        % get lfp derivative threshold
        lfp_diff_minbound = lfp_diff_mean - cfg_noise.amp_thr * lfp_diff_std;
        lfp_diff_maxbound = lfp_diff_mean + cfg_noise.amp_thr * lfp_diff_std;
        
        % save this data into the struct
        noisy_lfp_trials(i).lfp_diff_mean = lfp_diff_mean;
        noisy_lfp_trials(i).lfp_diff_std = lfp_diff_std;  
        
        % find the spectral power average
        arr_concat_trials_tfr = cat(3,concat_site_lfp_pow{:});
        pow_mean_f = nanmean(arr_concat_trials_tfr, 3);
        pow_std_f = nanstd(arr_concat_trials_tfr, 0, 3);
        noisy_lfp_trials(i).lfp_pow_mean = pow_mean_f;
        noisy_lfp_trials(i).lfp_pow_std = pow_std_f;
        
               
        % arrays to store the trails marked as noisy by each method
        noisy_trials_lfp_mean   = zeros(1,length(concat_site_lfp));
        noisy_trials_lfp_std    = zeros(1,length(concat_site_lfp));
        noisy_trials_lfp_diff   = zeros(1,length(concat_site_lfp));
        noisy_trials_lfp_pow    = zeros(1,length(concat_site_lfp));
        
        % now find the noisy trials 
        
        % loop through each trial lfp again
        for t = 1:length(site.trials)
            trial_lfp_std = std(concat_site_lfp{t});
            % a trial is noisy if any LFP datapoint in the trial is beyond 
            % threshold for N consecutive samples
            n = 0; 
            for d = concat_site_lfp{t}
                if d < lfp_raw_minbound || d > lfp_raw_maxbound
                    n = n + 1;
                    if n >= cfg_noise.amp_N
                        states_lfp(i).trial(t).noisy = 1;
                        noisy_trials_lfp_mean(t) = 1;
                        break;
                    end
                else
                    n = 0;
                end
            end
            % a trial is noisy if the LFP std of trial is beyond 2*std for the
            % site
            if trial_lfp_std > 2*site_lfp_std
                states_lfp(i).trials(t).noisy = 1;
                noisy_trials_lfp_std(t) = 1;
            end            
            % a trial is noisy if for more than N consecutive samples,
            % derivative of LFP is greater or less than mean +/- 2*std of
            % LFP derivative of all trials

            for d = concat_site_diff_lfp{t}
                if d > lfp_diff_maxbound || d < lfp_diff_minbound
                    n = n + 1; % increment counter
                    if n >= cfg_noise.diff_N
                        states_lfp(i).trial(t).noisy = 1;
                        noisy_trials_lfp_diff(t) = 1; 
                        break;
                    end
                else % reset counter
                    n = 0;
                end
            end
                 
            % a trial is noisy if the spectral power for 50% of frequency bins at any time bin 
            % is greater than mean power + 2*std of all trials at that time bin
            for tbin = 1:size(concat_site_lfp_pow{t}, 3)
                pow_t = concat_site_lfp_pow{t}(:,:,tbin);
                noisyfbins = sum(pow_t > pow_mean_f + cfg_noise.pow_thr*pow_std_f);
                if noisyfbins / size(concat_site_lfp_pow{t}, 2) > 0.5
                    states_lfp(i).trial(t).noisy = 1;
                    noisy_trials_lfp_pow(t) = 1;
                    break;
                end
            end
            
            
            % plot trials
            if cfg_noise.plottrials 
                noisy_trial_folder = fullfile(results_folder_noise, '/noisy_trials/', site.site_ID);
                if ~exist(noisy_trial_folder, 'dir')
                    mkdir(noisy_trial_folder)
                end

                non_noisy_trial_folder = fullfile(results_folder_noise, '/non_noisy_trials/', site.site_ID);
                if ~exist(non_noisy_trial_folder, 'dir')
                    mkdir(non_noisy_trial_folder)
                end

                figure; colormap jet;
                % plot raw lfp amplitude
                subplot(311)
                plot(concat_site_time{t}, concat_site_lfp{t})
                xlabel('time (s)');
                title(sprintf('LFP raw amplitude, excluded = %g', noisy_trials_lfp_mean(t)));
                line(xlim, [lfp_raw_minbound lfp_raw_minbound], 'color', 'r');
                line(xlim, [lfp_raw_maxbound lfp_raw_maxbound], 'color', 'r');
                % mark all states
                %states_names = cell(1, length(concat_site_states{t}));
                for s = 1:length(concat_site_states{t})
                    state = concat_site_states{t}(s);
                    state_onset = concat_site_states_onset{t}(s);
                    if state == 6 || state == 62 || state == 21
                        line([state_onset state_onset], ylim, 'color', 'k', 'linestyle', '--');
                    end
                end
                % plot raw lfp derivative
                subplot(312)
                plot(concat_site_time{t}, concat_site_diff_lfp{t})
                xlabel('time (s)');
                title(sprintf('LFP derivative, excluded = %g', noisy_trials_lfp_diff(t)));
                line(xlim, [lfp_diff_minbound lfp_diff_minbound], 'color', 'r');
                line(xlim, [lfp_diff_maxbound lfp_diff_maxbound], 'color', 'r');
                % mark all states
                %states_names = cell(1, length(concat_site_states{t}));
                for s = 1:length(concat_site_states{t})
                    state = concat_site_states{t}(s);
                    state_onset = concat_site_states_onset{t}(s);
                    if state == 6 || state == 62 || state == 21
                        line([state_onset state_onset], ylim, 'color', 'k', 'linestyle', '--');
                    end
                end
                % plot lfp power spectrogram
                if sum(strcmp(cfg_noise.methods, 'pow'))
                    subplot(313)                       
                    % make a tfr struct compatible with FT
                    noisy_tfr = struct();
                    noisy_tfr.powspctrm = concat_site_lfp_pow{t};
                    noisy_tfr.freq = site.trials(t).tfs.freq;
                    noisy_tfr.time = linspace(concat_site_time{t} (1), ...
                        concat_site_time{t} (end), length(concat_site_lfp_pow{t}));
                    text = strrep([site.site_ID, '_Trial_' num2str(t)], '_', '\_');
                    noisy_tfr.label = {text};
                    % plot TFR
                    cfg          = [];
                    cfg.zlim     = [0 3e-8];
                    cfg.channel  = noisy_tfr.label;
                    ft_singleplotTFR(cfg, noisy_tfr);
                    
                    % mark all states
                    %states_names = cell(1, length(concat_site_states{t}));
                    for s = 1:length(concat_site_states{t})
                        state = concat_site_states{t}(s);
                        state_onset = concat_site_states_onset{t}(s);
                        if state == 6 || state == 62 || state == 21
                            line([state_onset state_onset], ylim, 'color', 'k', 'linestyle', '--');
                        end
                    end
                    title(sprintf('LFP time freq spectrogram, excluded = %g', noisy_trials_lfp_pow(t)));

                end
                % save figure
                if states_lfp(i).trials(t).noisy
                    saveas(gca, fullfile(results_folder_noise, '/noisy_trials/', site.site_ID, ...
                        ['Trial_', num2str(t), '_raw_', num2str(noisy_trials_lfp_mean(t)), ...
                        '_derivative_', num2str(noisy_trials_lfp_diff(t)), ...
                        '_tfs_', num2str(noisy_trials_lfp_pow(t)), '.png']));
                else
                    saveas(gca, fullfile(results_folder_noise, '/non_noisy_trials/', site.site_ID, ...
                        ['Trial_', num2str(t), '_raw_', num2str(noisy_trials_lfp_mean(t)), ...
                        '_derivative_', num2str(noisy_trials_lfp_diff(t)), ...
                        '_tfs_', num2str(noisy_trials_lfp_pow(t)), '.png']));
                end
                close;
            end
            
        end
        
        % print information
        fprintf('Total number trials: %g \n', length(site.trials));
        fprintf('No of rejected trials, raw LFP mean: %g \n', sum(noisy_trials_lfp_mean));
        fprintf('No of rejected trials, raw LFP std: %g \n', sum(noisy_trials_lfp_std));
        fprintf('No of rejected trials, LFP derivative: %g \n', sum(noisy_trials_lfp_diff));
        fprintf('No of rejected trials, LFP power: %g \n', sum(noisy_trials_lfp_pow));
        noisy_trials = noisy_trials_lfp_mean | noisy_trials_lfp_std | noisy_trials_lfp_diff | noisy_trials_lfp_pow;
        nnoisy = sum(noisy_trials_lfp_mean | noisy_trials_lfp_std | noisy_trials_lfp_diff | noisy_trials_lfp_pow);
        fprintf('Total no:of noisy trials detected = %g \n', nnoisy);
        states_lfp(i).ntrails = length(site.trials);
        states_lfp(i).noisytrails = sum(noisy_trials_lfp_mean);
        
        % save into struct
        noisy_lfp_trials(i).raw_amp = find(noisy_trials_lfp_mean);
        noisy_lfp_trials(i).raw_std = find(noisy_trials_lfp_std);
        noisy_lfp_trials(i).raw_diff = find(noisy_trials_lfp_diff);
        noisy_lfp_trials(i).lfp_pow = find(noisy_trials_lfp_pow);
        noisy_lfp_trials(i).all_noisy = find(noisy_trials);
        
        % get blockwise statistics of noise rejected trials for printing
        % and visualzation
%         unq_blocks = unique(concat_blocks);
%         % number samples in each block
%         block_nsamples = zeros(1, length(unq_blocks));
%         % no of trials in each block
%         block_ntrials = zeros(1, length(unq_blocks));
%         % no of rejected trials in each block
%         block_nrejtrials = zeros(1, length(unq_blocks));
%         
%         for b = 1:length(unq_blocks)
%             block_trials_lfp = horzcat(concat_site_lfp{concat_blocks == unq_blocks(b)}) ;
%             block_nsamples(b) = length(block_trials_lfp);
% %             for t = 1:length(site.trials)
% %                 block_nsamples(b) = block_nsamples(b)...
% %                     + length(concat_site_lfp{t});
% %             end
%             block_ntrials(b) = sum(concat_blocks == unq_blocks(b));
%             block_nrejtrials(b) = sum(concat_blocks == unq_blocks(b) & ...
%                 noisy_trials == 1);
%         end
        
        % plots
        % plot the concatenated lfp for each site with mean and std
%         figure
%         hold on
%         plot(arr_concat_site_lfp)
%         plot(arr_concat_diff_lfp)
%         line(xlim, [site_lfp_mean site_lfp_mean], 'color', 'k');
%         line(xlim, [site_lfp_mean + 6*site_lfp_std site_lfp_mean + 6*site_lfp_std], 'color', 'r');
%         line(xlim, [site_lfp_mean - 6*site_lfp_std site_lfp_mean - 6*site_lfp_std], 'color', 'r');
%         line(xlim, [lfp_diff_mean lfp_diff_mean], 'color', 'k');
%         line(xlim, [lfp_diff_mean + 3*lfp_diff_std lfp_diff_mean + 3*lfp_diff_std], 'color', 'g');
%         line(xlim, [site_lfp_mean - 3*lfp_diff_std lfp_diff_mean - 3*lfp_diff_std], 'color', 'g');
%         title(sprintf('Site = %s, (ntrials = %g, naccepted = %g)', strrep(site.site_ID, '_', '\_'), t, length(site.trials) - states_lfp(i).noisytrails));
%         % divide blocks by drawing vertical lines
%         block_end_sample = 0;
%         for b = 1:length(unq_blocks)
%             
%             block_end_sample = block_end_sample + block_nsamples(b);
%             line([block_end_sample block_end_sample], ylim, 'color', 'k', ...
%                 'linestyle', '--');
%             block_ann = ['Block ', unq_blocks(b), ' ntrials = ' num2str(block_ntrials(b)), ' nrejected = ' num2str(block_nrejtrials(b))];
%             %annotation('textbox', [1 0.6 0.1 0.3], 'String', block_ann);
%         end
%         saveas(gca, fullfile(fig_folder_noise, [states_lfp(i).site_ID '_concat_LFP.png']));
        
        % plot the concatenated lfp derivative for each site with mean and std
%         figure
%         hold on
%         plot(arr_concat_diff_lfp)
%         line(xlim, [lfp_diff_mean lfp_diff_mean], 'color', 'k');
%         line(xlim, [lfp_diff_mean + 3*lfp_diff_std lfp_diff_mean + 3*lfp_diff_std], 'color', 'r');
%         line(xlim, [site_lfp_mean - 3*lfp_diff_std lfp_diff_mean - 3*lfp_diff_std], 'color', 'r');
%         title(sprintf('Site = %s, (ntrials = %g, naccepted = %g)', strrep(site.site_ID, '_', '\_'), t, length(site.trials) - session_lfp(i).noisytrails));
%         % divide blocks by drawing vertical lines
%         block_end_sample = 0;
%         for b = 1:length(unq_blocks)
%             
%             block_end_sample = block_end_sample + block_nsamples(b);
%             line([block_end_sample block_end_sample], ylim, 'color', 'k', ...
%                 'linestyle', '--');
%             block_ann = ['Block ', unq_blocks(b), ' ntrials = ' num2str(block_ntrials(b)), ' nrejected = ' num2str(block_nrejtrials(b))];
%             %annotation('textbox', [1 0.6 0.1 0.3], 'String', block_ann);
%         end
%         saveas(gca, fullfile(fig_folder_noise, [session_lfp(i).site_ID '_concat_diff_LFP.png']));
        
        % plot the lfp TFR for each site with mean
%         figure
%         cfg = [];
%         %cfg.baseline     = [-3e-27 3e-27];                 % -400ms to -100ms before the state onset
%         cfg.baselinetype = 'no';%'absolute';  
%         cfg.maskstyle    = 'saturation';
%         cfg.zlim         = [-3e-10 3e-10];
%         cfg.masknans     = 'yes';
%         cfg.channel  = site.site_ID;
%         ft_singleplotTFR(cfg, TFR_avg);
%         xlabel('Time'), ylabel('Frequency');
%         title(sprintf('Mean power (Site = %s, ntrials = %g)', strrep(site.site_ID, '_', '\_'), t));
% 
%         saveas(gca, fullfile(fig_folder_noise, [states_lfp(i).site_ID '_mean_pow.png']));
        
        
    end
    
    filt_session_lfp = states_lfp;
    % save data
    save (fullfile(results_folder, 'states_lfp'), 'states_lfp');
    save (results_folder_noise, 'noisy_lfp_trials');
    
    fprintf('Noise rejection in LFP done. \n');
    fprintf('=============================================================\n');

end
