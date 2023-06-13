function [site_lfp] = lfp_tfa_reject_noisy_lfp_trials( site_lfp, cfg_noise )
%lfp_tfa_reject_noisy_trials  - marks noisy lfp trials based on threshold 
%values of raw lfp amplitude and raw lfp derivative
%
% USAGE:
%	[site_lfp] = lfp_tfa_reject_noisy_lfp_trials( site_lfp, cfg_noise )
%
% INPUTS:
%		site_lfp         - struct containing trial-wise LFP power
%		spectrograms, see lfp_tfa_process_lfp, lfp_tfa_compute_site_tfr
%		cfg_noise        - configuration for noise rejection, see
%		lfp_tfa_cfg.noise in settings/lfp_tfa_settings_example
%           methods      : methods to be used for noise rejection ('raw',
%                        'diff', 'std', 'pow')
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
%		site_lfp      - same as input struct site_lfp, but with an
%                     additional field which marks noisy trials
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_process_lfp,
% lfp_tfa_compute_site_tfr

    results_folder = cfg_noise.results_folder;
    results_folder_noise = fullfile(results_folder, 'LFP_Noise_Rejection');
    if ~exist(results_folder_noise, 'dir')
        mkdir(results_folder_noise);
    end
    
    fprintf('Noise Rejection in LFP trials\n');
    fprintf('-----------------------------\n');
    
    
    % struct to save info about noisz trials
    noisy_lfp_trials = struct();

    % first, calculate the mean and std of length-wise concatenated LFP from
    % each site
    % read in lfp data for each site
    %for i = 1:length(site_lfp)
        %fprintf('Processing site lfp %s ...\n', site_lfp.site_ID);
        % cell arrays to store LFP data for all trials
        concat_site_lfp = cell(1,length(site_lfp.trials)); % store LFP raw data for each trial
        concat_site_time = cell(1, length(site_lfp.trials)); % store LFP timestamps for each trial
        concat_site_diff_lfp = cell(1,length(site_lfp.trials)); % store derivative of LFP for each trial
        concat_blocks = zeros(1, length(site_lfp.trials)); % store blocks for each trial
        concat_site_lfp_pow = cell(1,length(site_lfp.trials)); % store LFP spectral power for each trial
        concat_site_states = cell(1,length(site_lfp.trials)); % store state ids for each trial
        concat_site_states_onset = cell(1,length(site_lfp.trials)); % store state onset times for each trial

        for t = 1:length(site_lfp.trials)
            % read lfp data for each trial
            trial = site_lfp.trials(t);
            block = trial.block;
            lfp_data = trial.lfp_data;
            lfp_tfs = trial.tfs;
            lfp_time = trial.time;
            % create a new field to indicate if the trial is noisy
            site_lfp.trials(t).noisy = 0;

            trial_start_state = trial.states(1);
            trial_end_state = trial.states(end);
            time_index=lfp_time >= trial.trialperiod(1) & lfp_time <= trial.trialperiod(2);
            trial_time = lfp_time(time_index);
            trial_lfp_raw = lfp_data(time_index);
            if isfield(lfp_tfs,'time') %% basically means, this has been resampled
                time_index_tfs=lfp_tfs.time >= trial.trialperiod(1) & lfp_tfs.time <= trial.trialperiod(2);
                trial_lfp_pow = lfp_tfs.powspctrm(1,:,time_index_tfs);
                trial_pow_tbins = lfp_tfs.time(time_index_tfs);
            else
                trial_lfp_pow = lfp_tfs.powspctrm(1,:,time_index);
                trial_pow_tbins = trial_time;
            end
            
            trial_pow_timestep = trial_pow_tbins(2) - trial_pow_tbins(1);
            trial_pow_fbins = lfp_tfs.freq;
            tsample = 1/trial.fsample;

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

        % concatenate all derivatives to a 1d array
        arr_concat_diff_lfp = horzcat(concat_site_diff_lfp{:});
        % now get the mean and std of LFP for each site
        lfp_diff_mean = nanmean(arr_concat_diff_lfp);
        lfp_diff_std = nanstd(arr_concat_diff_lfp);
        % get lfp derivative threshold
        lfp_diff_minbound = lfp_diff_mean - cfg_noise.diff_thr * lfp_diff_std;
        lfp_diff_maxbound = lfp_diff_mean + cfg_noise.diff_thr * lfp_diff_std;




        % arrays to store the trails marked as noisy by each method
        noisy_trials_lfp_mean   = zeros(1,length(concat_site_lfp));
        noisy_trials_lfp_std    = zeros(1,length(concat_site_lfp));
        noisy_trials_lfp_diff   = zeros(1,length(concat_site_lfp));
        noisy_trials_lfp_pow    = zeros(1,length(concat_site_lfp));

        % now find the noisy trials 

        % loop through each trial lfp again
        for t = 1:length(site_lfp.trials)
            trial_lfp_std = std(concat_site_lfp{t});
            % a trial is noisy if any LFP datapoint in the trial is beyond 
            % threshold for N consecutive samples
            n = 0; 
            for d = concat_site_lfp{t}
                if d < lfp_raw_minbound || d > lfp_raw_maxbound
                    n = n + 1;
                    if n >= cfg_noise.amp_N
                        site_lfp.trials(t).noisy = 1;
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
                site_lfp.trials(t).noisy = 1;
                noisy_trials_lfp_std(t) = 1;
            end            
            % a trial is noisy if for more than N consecutive samples,
            % derivative of LFP is greater or less than mean +/- 2*std of
            % LFP derivative of all trials

            for d = concat_site_diff_lfp{t}
                if d > lfp_diff_maxbound || d < lfp_diff_minbound
                    n = n + 1; % increment counter
                    if n >= cfg_noise.diff_N
                        site_lfp.trials(t).noisy = 1;
                        noisy_trials_lfp_diff(t) = 1; 
                        break;
                    end
                else % reset counter
                    n = 0;
                end
            end

            % a trial is noisy if the spectral power for 50% of frequency bins at any time bin 
            % is greater than mean power + 2*std of all trials at that time bin
            pow_mean_f = nanmean(concat_site_lfp_pow{t}, 3);
            pow_std_f = nanstd(concat_site_lfp_pow{t}, 0, 3);
            for tbin = 1:size(concat_site_lfp_pow{t}, 3)
                pow_t = concat_site_lfp_pow{t}(:,:,tbin);
                noisyfbins = sum(pow_t > pow_mean_f + cfg_noise.pow_thr*pow_std_f);
                if noisyfbins / size(concat_site_lfp_pow{t}, 2) > 0.5
                    site_lfp.trials(t).noisy = 1;
                    noisy_trials_lfp_pow(t) = 1;
                    break;
                end
            end


            % plot trials
            if cfg_noise.plottrials 
                noisy_trial_folder = fullfile(results_folder_noise, '/noisy_trials/', site_lfp.site_ID);
                if ~exist(noisy_trial_folder, 'dir')
                    mkdir(noisy_trial_folder)
                end

                non_noisy_trial_folder = fullfile(results_folder_noise, '/non_noisy_trials/', site_lfp.site_ID);
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

                    noisy_tfr.freq = site_lfp.trials(t).tfs.freq;
                    noisy_tfr.time = linspace(concat_site_time{t} (1), ...
                        concat_site_time{t} (end), length(concat_site_lfp_pow{t}));
                    % normalize
                    noisy_tfr.powspctrm_norm = (noisy_tfr.powspctrm - ...
                        repmat(pow_mean_f, 1, 1, size(noisy_tfr.powspctrm, 3))) ./ ...
                        repmat(pow_std_f, 1, 1, size(noisy_tfr.powspctrm, 3));

                    text = strrep([site_lfp.site_ID, '_Trial_' num2str(t)], '_', '\_');
                    noisy_tfr.label = {text};
                    % plot TFR
                    cfg          = [];
                    cfg.zlim     = [0 3e-8];
                    cfg.channel  = noisy_tfr.label;
                    imagesc(noisy_tfr.time, [1:numel(noisy_tfr.freq)], squeeze(noisy_tfr.powspctrm_norm));
                    axis xy; colorbar;
                    set(gca, 'ytick', ([1:8:numel(noisy_tfr.freq)]));
                    set(gca, 'YTickLabel', noisy_tfr.freq([1:8:numel(noisy_tfr.freq)]));
                    %ft_singleplotTFR(cfg, noisy_tfr);

                    % mark all states
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
                if site_lfp.trials(t).noisy
                    saveas(gca, fullfile(results_folder_noise, '/noisy_trials/', site_lfp.site_ID, ...
                        ['Trial_', num2str(t), '_raw_', num2str(noisy_trials_lfp_mean(t)), ...
                        '_derivative_', num2str(noisy_trials_lfp_diff(t)), ...
                        '_tfs_', num2str(noisy_trials_lfp_pow(t)), '.png']));
                else
                    saveas(gca, fullfile(results_folder_noise, '/non_noisy_trials/', site_lfp.site_ID, ...
                        ['Trial_', num2str(t), '_raw_', num2str(noisy_trials_lfp_mean(t)), ...
                        '_derivative_', num2str(noisy_trials_lfp_diff(t)), ...
                        '_tfs_', num2str(noisy_trials_lfp_pow(t)), '.png']));
                end
                close;
            end

        end

        % plot site-wise summary plots
        unq_blocks = unique([site_lfp.trials.block]);
        concat_blocks = [site_lfp.trials.block];
        block_nsamples = zeros(1, length(unq_blocks));
        for b = 1:length(unq_blocks)
            block_trials_lfp = horzcat(concat_site_lfp{concat_blocks == unq_blocks(b)}) ;
            block_nsamples(b) = length(block_trials_lfp);
        end
        figure;
        % concatenated raw LFP and LFP derivative
        subplot(321);
        hold on
        plot(arr_concat_site_lfp)
        line(xlim, [site_lfp_mean site_lfp_mean], 'color', 'k');
        line(xlim, [lfp_raw_maxbound lfp_raw_maxbound], 'color', 'r');
        line(xlim, [lfp_raw_minbound lfp_raw_minbound], 'color', 'r');
        title('LFP amplitude');
        %title(sprintf('Site = %s, (ntrials = %g, naccepted = %g)', strrep(site.site_ID, '_', '\_'), t, length(site.trials) - states_lfp(i).noisytrails));
        % divide blocks by drawing vertical lines
        block_end_sample = 0;
        for b = 1:length(unq_blocks)

            block_end_sample = block_end_sample + block_nsamples(b);
            line([block_end_sample block_end_sample], ylim, 'color', 'k', ...
                'linestyle', '--');
            % to use in future
            %block_ann = ['Block ', unq_blocks(b), ' ntrials = ' num2str(block_ntrials(b)), ' nrejected = ' num2str(block_nrejtrials(b))];
            %annotation('textbox', [1 0.6 0.1 0.3], 'String', block_ann);
        end
        %saveas(gca, fullfile(results_folder_noise, [states_lfp(i).site_ID '_concat_LFP.png']));
        subplot(323);
        hold on
        plot(arr_concat_diff_lfp)
        line(xlim, [lfp_diff_mean lfp_diff_mean], 'color', 'k');
        line(xlim, [lfp_diff_mean + cfg_noise.diff_thr*lfp_diff_std lfp_diff_mean + cfg_noise.diff_thr*lfp_diff_std], 'color', 'r');
        line(xlim, [lfp_diff_minbound lfp_diff_minbound], 'color', 'r');
        % divide blocks by drawing vertical lines
        block_end_sample = 0;
        for b = 1:length(unq_blocks)

            block_end_sample = block_end_sample + block_nsamples(b);
            line([block_end_sample block_end_sample], ylim, 'color', 'k', ...
                'linestyle', '--');
            % to use in future
            %block_ann = ['Block ', unq_blocks(b), ' ntrials = ' num2str(block_ntrials(b)), ' nrejected = ' num2str(block_nrejtrials(b))];
            %annotation('textbox', [1 0.6 0.1 0.3], 'String', block_ann);
        end
        subplot(322)
        hold on;
        hr = histogram(arr_concat_site_lfp);
        
        subplot(324)
        hold on;
        hd = histogram(arr_concat_diff_lfp);
        
        % plot concatenated spectrogram
        arr_concat_lfp_pow = cat(3, concat_site_lfp_pow{:});
        x = 1:trial_pow_timestep: trial_pow_timestep * (size(arr_concat_lfp_pow, 3) - 1);
        subplot(3,2,[5 6])
        imagesc(x, trial_pow_fbins(end:1), squeeze(arr_concat_lfp_pow));
        axis xy

        saveas(gca, fullfile(results_folder_noise, [site_lfp.site_ID '_concat_raw_and_deriv_LFP.png']));

        % print information
        fprintf('Total number trials: %g \n', length(site_lfp.trials));
        fprintf('No of rejected trials, raw LFP mean: %g \n', sum(noisy_trials_lfp_mean));
        fprintf('No of rejected trials, raw LFP std: %g \n', sum(noisy_trials_lfp_std));
        fprintf('No of rejected trials, LFP derivative: %g \n', sum(noisy_trials_lfp_diff));
        fprintf('No of rejected trials, LFP power: %g \n', sum(noisy_trials_lfp_pow));
        noisy_trials = noisy_trials_lfp_mean | noisy_trials_lfp_std | noisy_trials_lfp_diff | noisy_trials_lfp_pow;
        nnoisy = sum(noisy_trials_lfp_mean | noisy_trials_lfp_std | noisy_trials_lfp_diff | noisy_trials_lfp_pow);
        fprintf('Total no:of noisy trials detected = %g \n', nnoisy);
        site_lfp.ntrails = length(site_lfp.trials);
        site_lfp.noisytrails = sum(noisy_trials);

    %fprintf('Noise rejection in LFP done. \n');
    %fprintf('=============================================================\n');

end

