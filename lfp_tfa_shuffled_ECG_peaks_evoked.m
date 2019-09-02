function shuffled_ecg_evoked = lfp_tfa_shuffled_ECG_peaks_evoked( trials_signal, site_ID, cfg_ecg, signal_type )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

nshuffles = 50;
distplot = true;

shuffled_ecg_evoked.lfp = {};
shuffled_ecg_evoked.lfp_time = {};
shuffled_ecg_evoked.state = cfg_ecg{1};
shuffled_ecg_evoked.state_name = cfg_ecg{2};

signal = {};
timestamps = {};
ecg_peaks = {};
ecg_p2pt = [];

for t = 1:length(trials_signal)
    if trials_signal(t).completed == 0
        continue;
    end
    % check if ECG spike data exists for this trial
    if isempty(trials_signal(t).ECG_spikes)
        continue;
    end
    % get the LFP samples and timestamps for the trial
    if strcmp(signal_type, 'lfp')
        signal = [signal, trials_signal(t).lfp_data];
    elseif strcmp(signal_type, 'ecg')
        if isempty(trials_signal(t).ecg_data)
            continue;
        else
            signal = [signal, trials_signal(t).ecg_data];
        end
    end
    
    timestamps = [timestamps, trials_signal(t).time];
    ecg_peaks = [ecg_peaks, logical(trials_signal(t).ECG_spikes)];
    ecg_p2pt = [ecg_p2pt, diff(trials_signal(t).time(...
        logical(trials_signal(t).ECG_spikes)))];
    
end

ts = trials_signal(1).tsample;

if distplot
    h = figure;
    histogram(ecg_p2pt);

end


if length(signal) ~= length(ecg_peaks)
    if strcmp(signal_type, 'lfp')
        signal_ts = trials_signal(1).tsample;
    elseif strcmp(signal_type, 'ecg')
        signal_ts = trials_signal(1).tsample;
    end
    if strcmp(signal_type, 'ecg')
        resampled_ecg_peaks = false(size(signal));
    end
end

ecg_evoked_lfp = struct();
for i = 1:nshuffles
    shuffled_ecg_p2pt = ecg_p2pt(randperm(length(ecg_p2pt)));
    ft_data_all = struct();
    ft_data_all.trial = cell(1, length(signal));
    ft_data_all.time = cell(1, length(signal));
    ft_data_all.label = {site_ID, 'Rpeak'};
    peak_idx = 1;
    for s = 1:length(signal)
        shuffled_ecg_peaks = false(size(ecg_peaks{s}));
        npeaks = sum(ecg_peaks{s});
        shuffled_ecg_peak_times = cumsum(shuffled_ecg_p2pt...
            (peak_idx:peak_idx + npeaks - 2));
        shuffled_ecg_peak_times = shuffled_ecg_peak_times(...
            shuffled_ecg_peak_times < timestamps{s}(end));
        shuffled_ECG_peak_samples = round(shuffled_ecg_peak_times / ts);
        shuffled_ecg_peaks(shuffled_ECG_peak_samples) = 1;
        peak_idx = peak_idx + npeaks - 1;
        ft_data_all.trial{s} = [signal{s}; shuffled_ecg_peaks];
        ft_data_all.time{s} = timestamps{s};
    end
%     % ECG peaks
%     shuffled_ecg_peaks = false(size(ecg_peaks));
%     % shuffled beat to beat interval
%     shuffled_ecg_p2pt = ecg_p2pt(randperm(length(ecg_p2pt)));
%     shuffled_ecg_peak_times = cumsum(shuffled_ecg_p2pt);
%     shuffled_ECG_peak_samples = round(shuffled_ecg_peak_times / ts);
%     shuffled_ecg_peaks(shuffled_ECG_peak_samples) = 1;

    

    % now find the evoked LFP
    % evoked LFP average
    cfg                 = [];
    cfg.timwin          = [cfg_ecg{3} cfg_ecg{4}];
    cfg.channel         = ft_data_all.label(1); % LFP channel
    cfg.spikechannel    = ft_data_all.label(2); % ECG peak

    ecg_based_sta       = ft_spiketriggeredaverage(cfg, ft_data_all);

    ecg_evoked_lfp(i).lfp = [];
    ecg_evoked_lfp(i).lfp_time = ecg_based_sta.time;
    ecg_evoked_lfp(i).mean = ecg_based_sta.avg;
end

% mean of all shuffles
shuffled_ecg_evoked.lfp = cat(1, ecg_evoked_lfp.mean);
shuffled_ecg_evoked.lfp_time = ecg_evoked_lfp(1).lfp_time;
shuffled_ecg_evoked.mean = nanmean(cat(1, ecg_evoked_lfp.mean), 1);
shuffled_ecg_evoked.std = nanstd(cat(1, ecg_evoked_lfp.mean), 0, 1);

clear ecg_evoked_lfp;


end

