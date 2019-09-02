function [ shuffled_ecg_evoked ] = ecg_lfp_rand_shuffle_ECG_peaks( session_ecg, block_ECG, lfp_tfa_cfg )
%lfp_tfa_compute_site_tfr Summary of this function goes here
%   Detailed explanation goes here

nshuffles = 50;
dist_plot = true;

%loop through each run
for b = 1:length(block_ECG)
    fprintf('Extracting ECG for block %g\n-----------------------\n', b);
    % concatenate all trials for this run
    concat_time = []; % to concatenate sample time
    concat_ecg = [];
    trials_time = [];
    trials_idx = find([session_ecg.trials.block] == b);
    for t = 1:length(trials_idx) % ignore first trial
        
        ts = session_ecg.trials(trials_idx(t)).tsample;
        if ~isempty(concat_time) 
            trial_timestamps = ...
                linspace(ts, ts*length(session_ecg.trials(trials_idx(t)).time), ...
                length(site_lfp.trials(trials_idx(t)).time)) + concat_time(end);
        else
            trial_timestamps = ...
                linspace(0, ts*(length(session_ecg.trials(trials_idx(t)).time)-1), ...
                length(site_lfp.trials(trials_idx(t)).time));
        end
        
        
%         trial_timestamps = ...
%             linspace(0, (length(site_lfp.trials(trials_idx(t)).time) - 1)*ts, ...
%                 length(site_lfp.trials(trials_idx(t)).time));

        concat_time = [concat_time trial_timestamps];
        concat_ecg = [concat_ecg, session_ecg.trials(trials_idx(t)).ecg_data];
        
        % initialize ECG spikes
        %site_lfp.trials(trials_idx(t)).ECG_spikes = false(size(trial_timestamps));
    end
    
    % get ECG timestamps for this block
    ECG_timestamps = block_ECG(b).Rpeak_t;
    % ECG beat to beat interval
    ECG_b2btime = diff(ECG_timestamps);
    % plot distribution
    if dist_plot
        figure;
        histogram(ECG_b2btime);
    end
    
    % randomly shuffle ECG beat to beat intervals
    ecg_evoked_lfp = struct();
    for i = 1:nshuffles
        % ECG peaks
        shuffled_ECG_peaks = false(size(concat_time));
        % shuffled beat to beat interval
        shuffled_ECG_b2btime = ECG_b2btime(randperm(length(ECG_b2btime)));
        shuffled_ECG_timestamps = cumsum(shuffled_ECG_b2btime);
        shuffled_ECG_peak_samples = shuffled_ECG_timestamps / ts;
        shuffled_ECG_peaks(shuffled_ECG_peak_samples) = 1;
        
        ft_data_all = struct();
        ft_data_all.trial = [concat_ecg; shuffled_ECG_peaks];
        ft_data_all.time = concat_time;
        ft_data_all.labels = {'ECG_raw', 'Rpeak'};
        
        % now find the evoked LFP
        % evoked LFP average
        cfg                 = [];
        cfg.keeptrials      = 'yes';
        cfg.timwin          = [-0.4 0.4];
        cfg.channel         = ft_data_all.label(1); % LFP chan
        cfg.spikechannel    = ft_data_all.label(2); % ECG peak
        cfg.latency         = 'maxperiod';

        ecg_based_sta       = ft_spiketriggeredaverage(cfg, ft_data_all);

        ecg_evoked_lfp(i).lfp = ecg_based_sta.trial;
        ecg_evoked_lfp(i).lfp_time = ecg_based_sta.time;
        ecg_evoked_lfp(i).mean = ecg_based_sta.avg;
        ecg_evoked_lfp(i).std = permute(nanstd(ecg_based_sta.trial, 0, 1), [2 3 1]);
    end
    
    % mean of all shuffles
    shuffled_ecg_evoked.lfp = cat(1, ecg_evoked_lfp.lfp);
    shuffled_ecg_evoked.lfp_time = ecg_evoked_lfp(1).lfp_time;
    shuffled_ecg_evoked.mean = nanmean(cat(1, ecg_evoked_lfp.mean), 1);
    shuffled_ecg_evoked.std = nanstd(cat(1, ecg_evoked_lfp.mean), 0, 1);
    
    clear ecg_evoked_lfp;
    
    % plot
    plottitle = 'Shuffled ECG evoked response for ' + session_ecg.session;
    results_file = fullfile(lfp_tfa_cfg.analyse_lfp_folder, ...
        'ECG Analysis', 'Shuffled_ECG_evoked_LFP', session_ecg.session + '.png');
    lfp_tfa_plot_evoked_lfp( shuffled_ecg_evoked, lfp_tfa_cfg, plottitle, results_file )
        
end



end
