function shuffled_ecg_evoked = lfp_tfa_shuffled_Rpeak_evoked_ECG( trials_ecg, cfg_ecg, nshuffles )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%nshuffles = 100;
distplot = false;

shuffled_ecg_evoked.lfp = {};
shuffled_ecg_evoked.lfp_time = {};
shuffled_ecg_evoked.state = cfg_ecg{1};
shuffled_ecg_evoked.state_name = cfg_ecg{2};

ECG_raw = [];
%timestamps = [];
ecg_peaks = [];
ecg_p2pt = [];

for t = 1:length(trials_ecg)
    if trials_ecg(t).completed == 0
        continue;
    end
    % check if ECG spike data exists for this trial
    if isempty(trials_ecg(t).ECG_spikes)
        continue;
    end
    % get the ECG raw data for the trial
    if isempty(trials_ecg(t).ecg_data)
        continue;
    else
        trial_ECG_raw = trials_ecg(t).ecg_data;
        trial_ECG_raw(~trials_ecg(t).ECG_valid) = [];
        ECG_raw = [ECG_raw, trial_ECG_raw];
    end
    
    %timestamps = [timestamps, trials_ecg(t).time];
    ecg_peaks = [ecg_peaks, ...
        logical(trials_ecg(t).ECG_spikes(trials_ecg(t).ECG_valid))];
    ecg_p2pt = [ecg_p2pt, ...
        trials_ecg(t).ECG_b2btime(trials_ecg(t).ECG_spikes & ...
        trials_ecg(t).ECG_valid)];
    
end

mean_ecg_p2pt = mean(ecg_p2pt);

ts = trials_ecg(1).tsample;
timestamps = linspace(0, ts*(length(ECG_raw) - 1), length(ECG_raw));

if distplot
    h = figure;
    histogram(ecg_p2pt);

end
    
Rpeak_evoked_ECG = struct();
for i = 1:nshuffles
    fprintf('Shuffle %g\n', i);
    shuffled_ecg_p2pt = ecg_p2pt(randperm(length(ecg_p2pt)));
    peak_idx = 1;
    shuffled_ecg_peaks = false(size(ecg_peaks));
    npeaks = sum(ecg_peaks);
    shuffled_ecg_peak_times = cumsum(shuffled_ecg_p2pt...
        (peak_idx:peak_idx + npeaks - 2));
    shuffled_ecg_peak_times = shuffled_ecg_peak_times(...
        shuffled_ecg_peak_times < timestamps(end));
    shuffled_ECG_peak_samples = round(shuffled_ecg_peak_times / ts);
    % randomly shift shuffled ECG peaks
    shuffled_ECG_peak_samples = ...
        shuffled_ECG_peak_samples + ...
        round((2*rand(size(shuffled_ECG_peak_samples)) - 1) * (mean_ecg_p2pt/ts));  
    shuffled_ECG_peak_samples(shuffled_ECG_peak_samples > ...
        length(ecg_peaks) | shuffled_ECG_peak_samples < 1) = [];
    % shuffled ECG peaks
    shuffled_ecg_peaks(shuffled_ECG_peak_samples) = 1;
    
    ft_data_all = struct();
    ft_data_all.trial = {[ECG_raw; shuffled_ecg_peaks]};
    ft_data_all.time = {timestamps};
    ft_data_all.label = {'ECG_raw', 'Rpeak'};
    
    % now find the evoked LFP
    % evoked LFP average
    cfg                 = [];
    cfg.timwin          = [cfg_ecg{3} cfg_ecg{4}];
    cfg.channel         = ft_data_all.label(1); % LFP channel
    cfg.spikechannel    = ft_data_all.label(2); % ECG peak

    ecg_based_sta       = ft_spiketriggeredaverage(cfg, ft_data_all);
    
    Rpeak_evoked_ECG(i).ecg_time = ecg_based_sta.time;
    Rpeak_evoked_ECG(i).mean = ecg_based_sta.avg;    
    
%     for s = 1:length(ECG_raw)
%         shuffled_ecg_peaks = false(size(ecg_peaks{s}));
%         npeaks = sum(ecg_peaks{s});
%         shuffled_ecg_peak_times = cumsum(shuffled_ecg_p2pt...
%             (peak_idx:peak_idx + npeaks - 2));
%         shuffled_ecg_peak_times = shuffled_ecg_peak_times(...
%             shuffled_ecg_peak_times < timestamps{s}(end));
%         shuffled_ECG_peak_samples = round(shuffled_ecg_peak_times / ts);
%         % randomly shift shuffled ECG peaks
%         shuffled_ECG_peak_samples = ...
%             shuffled_ECG_peak_samples + ...
%             round((2*rand(size(shuffled_ECG_peak_samples)) - 1) * (mean_ecg_p2pt/ts));  
%         shuffled_ECG_peak_samples(shuffled_ECG_peak_samples > ...
%             length(ecg_peaks{s}) | shuffled_ECG_peak_samples < 1) = [];
%         shuffled_ecg_peaks(shuffled_ECG_peak_samples) = 1;
%         peak_idx = peak_idx + npeaks - 1;
%         ft_data_all.trial{s} = [ECG_raw{s}; shuffled_ecg_peaks];
%         ft_data_all.time{s} = timestamps{s};
%     end
%     % ECG peaks
%     shuffled_ecg_peaks = false(size(ecg_peaks));
%     % shuffled beat to beat interval
%     shuffled_ecg_p2pt = ecg_p2pt(randperm(length(ecg_p2pt)));
%     shuffled_ecg_peak_times = cumsum(shuffled_ecg_p2pt);
%     shuffled_ECG_peak_samples = round(shuffled_ecg_peak_times / ts);
%     shuffled_ecg_peaks(shuffled_ECG_peak_samples) = 1;

end



shuffled_ecg_evoked.lfp = cat(1, Rpeak_evoked_ECG.mean);
shuffled_ecg_evoked.dimord = 'nshuffles_time';
shuffled_ecg_evoked.lfp_time = Rpeak_evoked_ECG(1).ecg_time;
shuffled_ecg_evoked.mean = nanmean(cat(1, Rpeak_evoked_ECG.mean), 1);
shuffled_ecg_evoked.std = nanstd(cat(1, Rpeak_evoked_ECG.mean), 0, 1);

% mean of all shuffles
% shuffled_ecg_evoked.lfp = cat(1, shuffled_Rpeak_evoked_ECG.mean);
% shuffled_ecg_evoked.lfp_time = shuffled_Rpeak_evoked_ECG(1).lfp_time;
% shuffled_ecg_evoked.mean = nanmean(cat(1, shuffled_Rpeak_evoked_ECG.mean), 1);
% shuffled_ecg_evoked.std = nanstd(cat(1, shuffled_Rpeak_evoked_ECG.mean), 0, 1);

clear ecg_evoked_lfp;


end

