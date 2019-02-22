function [sites, noisy_lfp_trials] = rejectNoisyLFP( sites, cfg_nr )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % RejectNoisyTrials.m
    % Script to reject noisy trials
    % Method 1: A trial is rejected if 10 consecutive LFP samples are beyond 0.001 and
    % 99.999 percentile of all LFP samples across all trials from a site
    % Method 2: A trial is rejected LFP standard deviation for a trial is
    % greater than 2*std of LFP samples across all trials from a site
    % Method 3: A trial is rejected if 10 consecutive LFP samples are beyond 0.1 and
    % 99.9 percentile of all LFP samples across all trials from a site
    % Method 4: A trial is rejected if in the TFR at any time bin, more than
    % half of the frequency bins have a power larger than mean + 2*std of LFP
    % power over all trials 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %clear;
    close all;
    %fig_folder = 'C:\Data\MIP_timefreq_analysis\Figures\';
    
    %sessionName = sites(1).site_ID(1:end-8);
    %fig_folder_noise = fullfile(fig_folder, sessionName, date, 'LFP_Noise_Rejection');
    %mkdir(fig_folder_noise);
    
    % which method(s) to use for noise rejection
    %methods = [3];
    % if configuration is missing assign default configuration
    if isempty(cfg_nr) % empty configuration structure
        % create a configuration structurewith default values
        cfg_nr = struct();
        % methods to be used
        cfg_nr.methods = {'amp', 'std', 'diff'};
        % threshold for lfp raw amplitude
        cfg_nr.amp_thr = [0.001 99.999];
        % number of consecutive samples beyond threshold to be considered
        cfg_nr.amp_N = 10;
        % no of standard deviations of trial w.r.t complete LFP std
        cfg_nr.std_thr = 4;
        % threshold for lfp derivative in percentile
        cfg_nr.diff_thr = [0.1 99.9];
        % number of consecutive samples beyond threshold to be considered
        cfg_nr.diff_N = 10;
    end        
    
    % struct to store information about noisy lfps
    noisy_lfp_trials = struct();

    % read in lfp data for each site
    for i = 1:length(sites)
        site = sites(i);
        noisy_lfp_trials(i).site = site.site_ID;
        % cell arrays to store LFP raw data, derivative, blocks and LFP power for all trials
        concat_site_lfp = cell(1,length(site.trial));
        concat_site_diff_lfp = cell(1,length(site.trial)); % store derivative of LFP for each trial
        concat_blocks = zeros(1, length(site.trial));
        concat_site_lfp_pow = cell(1,length(site.trial));
        %concat_site_lfp_pow = cell(1,length(site.trial));
        % read lfp data for each trial
        for t = 1:length(site.trial)
            % create a new field to indicate if the trial is noisy
            sites(i).trial(t).noisy = 0;
            trial = site.trial(t);
            if trial.success
                block = trial.block;
                lfp_data = trial.LFP;
                fs = trial.TDT_LFPx_SR;
                ts = 1/fs;
                timestamps = linspace(0, ts*(length(lfp_data)), length(lfp_data));
                
    %             % structure to store state specific time and lfp
    %             state_time_lfp = [];
    %             for state = trial.states
    %                 state_time_lfp = [state_time_lfp; [timestamps(state.start_s:state.end_s); ...
    %                     lfp_data(state.start_s:state.end_s)]'];
    %                 %trial_lfp = [trial_lfp state_lfp(:,2)']; 
    %             end
    %             % get data with unique timestamps
    %             unq_state_time_lfp  = unique(state_time_lfp, 'rows');
    %             trial_lfp = unq_state_time_lfp(:,2)';
    %             timestamps = unq_state_time_lfp(:,1)';


                % save lfp for this trial into cell array
                concat_site_lfp{t} = lfp_data;
                concat_site_diff_lfp{t} = [nan diff(lfp_data)]; 
                concat_blocks(t) = block;
                    
                % calculate the TFR if method 'pow' is active
                if sum(strcmp(cfg_nr.methods, 'pow'))
                    
                    % per site time frequency analysis using FT for noise rejection
                    ft_lfp_data             = struct();
                    ft_lfp_data.trial       = {lfp_data};
                    ft_lfp_data.time        = {0:ts:ts*(length(lfp_data)-1)};
                    ft_lfp_data.fsample     = fs;
                    ft_lfp_data.sampleinfo  = [0 ts*(length(lfp_data)-1)];
                    ft_lfp_data.label       = {site.site_ID};

                    % configuration for TFR
                    cfg              = [];
                    cfg.output       = 'pow';
                    cfg.method       = 'mtmconvol';
                    cfg.taper        = 'hanning';
                    cfg.pad          = 'nextpow2';
                    cfg.foi          = 2:2:100;             % analysis 2 to 100 Hz in steps of 2 Hz
                    cfg.t_ftimwin    = ones(length(cfg.foi),1).*500*ts;    % length of time window = 0.5 sec
                    cfg.toi          = 0:50*ts:ts*(length(lfp_data)-1);% time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
                    cfg.channel      = ft_lfp_data.label;
                    TFR_hann_fixed   = ft_freqanalysis(cfg, ft_lfp_data);

                    % subtract the baseline from spectrogram
                    % get the baseline period
                    baseline_ref_t = trial.states_onset(trial.states == cfg_nr.pow_refstate); % cue
                    baseline_period = [baseline_ref_t + cfg_nr.pow_baseline(1)...
                        baseline_ref_t + cfg_nr.pow_baseline(2)];
                    TFR_hann_fixed.powspctrm = freqBaseLineNorm(TFR_hann_fixed, baseline_period, 'absolute');
                    concat_site_lfp_pow{t} = TFR_hann_fixed.powspctrm;
                    
                end
            end
        end
        % concatenate all datapoints to a 1d array
        arr_concat_site_lfp = horzcat(concat_site_lfp{:});
        % now get the 0.10 and 99.9 percentile of LFP for each site
        site_lfp_min = prctile(arr_concat_site_lfp, cfg_nr.amp_thr(1));
        site_lfp_max = prctile(arr_concat_site_lfp, cfg_nr.amp_thr(2));
%         site_lfp_mean = mean(arr_concat_site_lfp);
        % standard deviation of lfp values for this site
        site_lfp_std = std(arr_concat_site_lfp);
        % save this data into the struct
%         session_lfp(i).lfp_mean = site_lfp_mean;
%         session_lfp(i).lfp_std = site_lfp_std;  
        
        % concatenate all derivatives to a 1d array
        arr_concat_diff_lfp = horzcat(concat_site_diff_lfp{:});
        % get the 0.1 and 99.9 percentile of derivative for each site
        lfp_diff_minbound = prctile(arr_concat_diff_lfp, cfg_nr.diff_thr(1));
        lfp_diff_maxbound = prctile(arr_concat_diff_lfp, cfg_nr.diff_thr(2));
        % now get the mean and std of LFP for each site
        lfp_diff_mean = nanmean(arr_concat_diff_lfp);
        lfp_diff_std = nanstd(arr_concat_diff_lfp);
        % save this data into the struct
%         session_lfp(i).lfp_diff_mean = lfp_diff_mean;
%         session_lfp(i).lfp_diff_std = lfp_diff_std;  
        
%         % concatenate all power spectrums to a 3d array
%         arr_concat_site_lfp_pow = vertcat(concat_site_lfp_pow{:});
%         % now get the mean and std of LFP for each site
%         site_lfp_pow_mean = mean(arr_concat_site_lfp_pow, 1);
%         site_lfp_pow_std = std(arr_concat_site_lfp_pow, 1);
%         TFR_avg = TFR_hann_fixed;
%         TFR_avg.powspctrm = site_lfp_pow_mean;       
%         
        if sum(strcmp(cfg_nr.methods, 'pow'))
            
            % find the spectral power average
            arr_concat_trials_tfr = cat(3,concat_site_lfp_pow{:});
            pow_mean_f = nanmean(arr_concat_trials_tfr, 3);
            pow_std_f = nanstd(arr_concat_trials_tfr, 0, 3);
            
        end
        
        % arrays to store noisy trials from each method separately
        noisy_trials_lfp_mean = zeros(1,length(concat_site_lfp));
        noisy_trials_lfp_std = zeros(1,length(concat_site_lfp));
        noisy_trials_lfp_diff = zeros(1,length(concat_site_lfp));
        noisy_trials_lfp_pow = zeros(1, length(concat_site_lfp));     
        % now find the noisy trials 
        % a trial is noisy if any LFP datapoint in the trial is beyond mean +/-
        % 3std of the concatenated LFP for that site
        % loop through each trial lfp again       
    
        for t = 1:length(site.trial)
            if ~isempty(concat_site_lfp{t})
                trial_lfp_std = std(concat_site_lfp{t});
                n = 0; 
                for d = concat_site_lfp{t}
                    if d < site_lfp_min || d > site_lfp_max
                        n = n + 1;
                        if n >= cfg_nr.amp_N
                            sites(i).trial(t).noisy = 1;
                            noisy_trials_lfp_mean(t) = 1;
                            break;
                        end
                    else
                        n = 0;
                    end
                end
%                 if sum((concat_site_lfp{t} > ones(1, length(concat_site_lfp{t})) * ...
%                         (site_lfp_mean + 6*site_lfp_std)) + ...
%                         (concat_site_lfp{t} < ones(1, length(concat_site_lfp{t})) * ...
%                         (site_lfp_mean - 6*site_lfp_std)))
%                     sites(i).trial(t).noisy = 1;
%                     noisy_trials_lfp_mean(t) = 1;
%                 end
                % a trial is noisy if the LFP std of trial is beyond 2*std for the
                % site
                if trial_lfp_std > cfg_nr.std_thr*site_lfp_std
                    sites(i).trial(t).noisy = 1;
                    noisy_trials_lfp_std(t) = 1;
                end            
                % a trial is noisy if for more than maxN consecutive samples,
                % derivative of LFP is greater or less than mean +/- 2*std of
                % LFP derivative of all trials
                maxN = 10; % maximum no of consecutive samples
                N = 0; % counter
                for d = concat_site_diff_lfp{t}
                    if d > lfp_diff_maxbound || d < lfp_diff_minbound
                        N = N + 1; % increment counter
                        if N >= cfg_nr.diff_N
                            sites(i).trial(t).noisy = 1;
                            noisy_trials_lfp_diff(t) = 1; 
                            break;
                        end
                    else % reset counter
                        N = 0;
                    end
                end
                
                % a trial is noisy if at any time point the spectral power for
                % more than 50% of freq bins is larger than mean power + 2*std 
                % of all trials at that time bin
                if sum(strcmp(cfg_nr.methods, 'pow'))
                    for tbin = 1:size(concat_site_lfp_pow{t}, 3)
                        pow_t = concat_site_lfp_pow{t}(:,:,tbin);
                        noisyfbins = sum(pow_t > pow_mean_f + cfg_nr.pow_thr*pow_std_f);
                        if noisyfbins / size(concat_site_lfp_pow{t}, 2) > 0.5
                            sites(i).trial(t).noisy = 1;
                            noisy_trials_lfp_pow(t) = 1;
                            break;
                        end
                    end
                end

                % a trial is noisy if the spectral power for 50% of frequency bins at any time bin 
                % is greater than mean power + 2*std of all trials at that time bin
    %             for tbin = 1:size(concat_site_lfp_pow{t}, 3)
    %                 pow_t = concat_site_lfp_pow{t}(:,:,tbin);
    %                 noisyfbins = sum(pow_t > site_lfp_pow_mean(:,:,tbin) + ...
    %                     2*site_lfp_pow_std(:,:,tbin));
    %                 if noisyfbins / size(concat_site_lfp_pow{t}, 2) > 0.5
    %                     sites(i).trial(t).noisy = 1;
    %                     noisy_trials_lfp_mean(t) = 1;
    %                     break;
    %                 end
    %             end
                
            end
        end
        
        noisy_trials_all = noisy_trials_lfp_mean | noisy_trials_lfp_std | ...
            noisy_trials_lfp_diff | noisy_trials_lfp_pow;
        fprintf('Total number trials: %g \n', length(site.trial));
        fprintf('No of rejected trials, lfp mean: %g \n', sum(noisy_trials_lfp_mean));
        fprintf('No of rejected trials, std: %g \n', sum(noisy_trials_lfp_std));
        fprintf('No of rejected trials, lfp diff: %g \n', sum(noisy_trials_lfp_diff));
        fprintf('No of rejected trials, lfp power: %g \n', sum(noisy_trials_lfp_pow));
        
        % save into struct
        noisy_lfp_trials(i).raw_amp = find(noisy_trials_lfp_mean);
        noisy_lfp_trials(i).raw_std = find(noisy_trials_lfp_std);
        noisy_lfp_trials(i).raw_diff = find(noisy_trials_lfp_diff);
        noisy_lfp_trials(i).lfp_pow = find(noisy_trials_lfp_pow);
        noisy_lfp_trials(i).all_noisy = find(noisy_trials_all);
        %session_lfp(i).ntrails = length(site.trials);
        %session_lfp(i).noisytrails = sum(noisy_trials_lfp_mean);
        
        % get blockwise statistics of noise rejected trials for printing
        % and visualzation
        unq_blocks = unique(concat_blocks);
        % number samples in each block
        block_nsamples = zeros(1, length(unq_blocks));
        % no of trials in each block
        block_ntrials = zeros(1, length(unq_blocks));
        % no of rejected trials in each block
        block_nrejtrials = zeros(1, length(unq_blocks));
        
        for b = 1:length(unq_blocks)
            block_trials_lfp = horzcat(concat_site_lfp{concat_blocks == unq_blocks(b)}) ;
            block_nsamples(b) = length(block_trials_lfp);
%             for t = 1:length(site.trials)
%                 block_nsamples(b) = block_nsamples(b)...
%                     + length(concat_site_lfp{t});
%             end
            block_ntrials(b) = sum(concat_blocks == unq_blocks(b));
            block_nrejtrials(b) = sum(concat_blocks == unq_blocks(b) & ...
                noisy_trials_lfp_mean == 1);
        end
        
        % plots
        % plot the concatenated lfp for each site with mean and std
        figure
        hold on
        nsubplots = length(cfg_nr.methods) - 1; % std method not plotted
        subplot(nsubplots, 1, 1)
        plot(arr_concat_site_lfp)
        line(xlim, [site_lfp_min site_lfp_min], 'color', 'r');
        line(xlim, [site_lfp_max site_lfp_max], 'color', 'r');
        title(sprintf('Raw LFP, noisy trials = %g', sum(noisy_trials_lfp_mean)));
        % divide blocks by drawing vertical lines
        block_end_sample = 0;
        for b = 1:length(unq_blocks)
            
            block_end_sample = block_end_sample + block_nsamples(b);
            line([block_end_sample block_end_sample], ylim, 'color', 'k', ...
                'linestyle', '--');
            block_ann = ['Block ', unq_blocks(b), ' ntrials = ' num2str(block_ntrials(b)), ' nrejected = ' num2str(block_nrejtrials(b))];
            %annotation('textbox', [1 0.6 0.1 0.3], 'String', block_ann);
        end
        subplot(nsubplots, 1, 2)
        plot(arr_concat_diff_lfp)
        line(xlim, [lfp_diff_minbound lfp_diff_minbound], 'color', 'g');
        line(xlim, [lfp_diff_maxbound lfp_diff_maxbound], 'color', 'g');
        title(sprintf('LFP derivative, noisy trials = %g', sum(noisy_trials_lfp_diff)));
%         subplot(322)
%         histfit(arr_concat_site_lfp);
%         subplot(324)
%         histfit(arr_concat_diff_lfp);
        if sum(strcmp(cfg_nr.methods, 'pow'))
            subplot(nsubplots, 1, 3)
            % concatenate all noisy trials
            arr_noisy_lfp_pow = cat(3, concat_site_lfp_pow{find(noisy_trials_lfp_pow)});
            % make a tfr struct compatible with FT
            noisy_tfr = TFR_hann_fixed;
            noisy_tfr.powspctrm = arr_noisy_lfp_pow;
            noisy_tfr.time = linspace(0, ts*length(arr_noisy_lfp_pow), length(arr_noisy_lfp_pow));
            % plot TFR
            cfg = [];
            cfg.zlim         = [min(arr_noisy_lfp_pow) max(arr_noisy_lfp_pow)];
            cfg.channel  = noisy_tfr.label;
            ft_singleplotTFR(cfg, noisy_tfr);
            title(sprintf('LFP power, noisy trials = %g', sum(noisy_trials_lfp_pow)));
        end
        %spectrogram(arr_concat_site_lfp, 100, 50, 100, 1e3, 'yaxis');
        %line(xlim, [site_lfp_mean site_lfp_mean], 'color', 'k');
        
        %line(xlim, [lfp_diff_mean lfp_diff_mean], 'color', 'k');
%         line(xlim, [lfp_diff_minbound lfp_diff_minbound], 'color', 'g');
%         line(xlim, [lfp_diff_maxbound lfp_diff_maxbound], 'color', 'g');
%         line(xlim, [site_lfp_mean site_lfp_mean], 'color', 'k');
%         line(xlim, [site_lfp_mean + 6*site_lfp_std site_lfp_mean + 6*site_lfp_std], 'color', 'r');
%         line(xlim, [site_lfp_mean - 6*site_lfp_std site_lfp_mean - 6*site_lfp_std], 'color', 'r');
%         line(xlim, [lfp_diff_mean lfp_diff_mean], 'color', 'k');
%         line(xlim, [lfp_diff_mean + 6*lfp_diff_std lfp_diff_mean + 6*lfp_diff_std], 'color', 'g');
%         line(xlim, [site_lfp_mean - 6*lfp_diff_std lfp_diff_mean - 6*lfp_diff_std], 'color', 'g');
        plottitle = ['Site = ', strrep(site.site_ID, '_', '\_'), '(ntrials = '...
            num2str(length(site.trial)), ' ncompleted = ', num2str(sum([site.trial.success])) ...
            ' naccepted = ', num2str(sum([site.trial.success]) - sum(noisy_trials_all)), ')']
%             , length(site.trial), ...
%             sum([site.trial.success]), sum([site.trial.success]) - sum(noisy_trials_lfp)]
%         title(sprintf('Site = %s, (ntrials = %g, ncompleted = %g naccepted = %g)', ...
%             strrep(site.site_ID, '_', '\_'), length(site.trial), ...
%             sum([site.trial.success]), sum([site.trial.success]) - sum(noisy_trials_lfp)));
        annotation('textbox', [0 0.9 1 0.1], 'String', plottitle, ...
            'EdgeColor', 'none', 'HorizontalAlignment', 'center');
        
        %saveas(gca, fullfile(fig_folder_noise, [sites(i).site_ID '_concat_LFP.png']));
        
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
%         saveas(gca, fullfile(fig_folder_noise, [session_lfp(i).site_ID '_mean_pow.png']));
        
        
    end
    
    %save sites sites;         

end

