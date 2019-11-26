function shuffled_LFP_evoked = lfp_tfa_shuffled_Rpeak_evoked_LFP( trials_lfp, cfg_state, nshuffles )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%nshuffles = 100;
distplot = false;

shuffled_LFP_evoked.lfp = {};
shuffled_LFP_evoked.lfp_time = {};
shuffled_LFP_evoked.state = cfg_state{1};
shuffled_LFP_evoked.state_name = cfg_state{2};

LFP_raw = [];
%timestamps = [];
ecg_peaks = [];
ecg_p2pt = [];

for t = 1:length(trials_lfp)
    if trials_lfp(t).completed == 0
        continue;
    end
    % check if ECG spike data exists for this trial
    if isempty(trials_lfp(t).ECG_spikes)
        continue;
    end
    % get the ECG raw data for the trial
    if isempty(trials_lfp(t).lfp_data)
        continue;
    else
        trial_LFP_raw = trials_lfp(t).lfp_data;
        %trial_LFP_raw(~trials_lfp(t).ECG_valid) = [];
        LFP_raw = [LFP_raw, trial_LFP_raw];
    end
    
    %timestamps = [timestamps, trials_ecg(t).time];
    ecg_peaks = [ecg_peaks, ...
        logical(trials_lfp(t).ECG_spikes)];
    ecg_p2pt = [ecg_p2pt, ...
        trials_lfp(t).ECG_b2btime(trials_lfp(t).ECG_spikes & ...
        ~isnan(trials_lfp(t).ECG_b2btime))];
    
end

mean_ecg_p2pt = mean(ecg_p2pt);

ts = trials_lfp(1).tsample;
timestamps = linspace(0, ts*(length(LFP_raw) - 1), length(LFP_raw));

if distplot
    h = figure;
    histogram(ecg_p2pt);

end
    
Rpeak_evoked_LFP = struct();
for i = 1:nshuffles
    fprintf('Shuffle %g\n', i);
    shuffled_ecg_p2pt = ecg_p2pt(randperm(length(ecg_p2pt)));
    peak_idx = 1;
    shuffled_ecg_peaks = false(size(ecg_peaks));
    npeaks = numel(shuffled_ecg_p2pt);
    shuffled_ecg_peak_times = cumsum(shuffled_ecg_p2pt...
        (peak_idx:peak_idx + npeaks - 1));
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
    ft_data_all.trial = {[LFP_raw; shuffled_ecg_peaks]};
    ft_data_all.time = {timestamps};
    ft_data_all.label = {'LFP', 'Rpeak'};
    
    % now find the evoked LFP
    % evoked LFP average
    cfg                 = [];
    cfg.timwin          = [cfg_state{3} cfg_state{4}];
    cfg.channel         = ft_data_all.label(1); % LFP channel
    cfg.spikechannel    = ft_data_all.label(2); % ECG peak

    ecg_based_sta       = ft_spiketriggeredaverage(cfg, ft_data_all);
    
    Rpeak_evoked_LFP(i).lfp_time = ecg_based_sta.time;
    Rpeak_evoked_LFP(i).mean = ecg_based_sta.avg;    
    
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



shuffled_LFP_evoked.lfp = cat(1, Rpeak_evoked_LFP.mean);
shuffled_LFP_evoked.dimord = 'nshuffles_time';
shuffled_LFP_evoked.lfp_time = Rpeak_evoked_LFP(1).lfp_time;
shuffled_LFP_evoked.mean = nanmean(cat(1, Rpeak_evoked_LFP.mean), 1);
shuffled_LFP_evoked.std = nanstd(cat(1, Rpeak_evoked_LFP.mean), 0, 1);

% mean of all shuffles
% shuffled_ecg_evoked.lfp = cat(1, shuffled_Rpeak_evoked_ECG.mean);
% shuffled_ecg_evoked.lfp_time = shuffled_Rpeak_evoked_ECG(1).lfp_time;
% shuffled_ecg_evoked.mean = nanmean(cat(1, shuffled_Rpeak_evoked_ECG.mean), 1);
% shuffled_ecg_evoked.std = nanstd(cat(1, shuffled_Rpeak_evoked_ECG.mean), 0, 1);

%clear ecg_evoked_lfp;


end

