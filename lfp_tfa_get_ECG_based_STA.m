function ecg_triggered_evoked = lfp_tfa_get_ECG_based_STA( trials, site_ID, cfg_ecg, signal_type )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ecg_triggered_evoked.lfp = {};
ecg_triggered_evoked.lfp_time = {};
ecg_triggered_evoked.state = cfg_ecg{1};
ecg_triggered_evoked.state_name = cfg_ecg{2};

% create a FT datatype
ft_data_all = struct();
ft_data_all.label = {site_ID, 'ECG'};
ft_data_all.time = {}; % timestamps
ft_data_all.trial = {}; % evoked LFP response

for t = 1:length(trials)

    trialperiod           = trials(t).trialperiod;
    
    % check if ECG spike data exists for this trial
    if isempty(trials(t).ECG_spikes)
        continue;
    end
    
    % get the LFP samples and timestamps for the trial
    if strcmp(signal_type, 'lfp')
        signal = trials(t).lfp_data;
    elseif strcmp(signal_type, 'ecg')
        signal = trials(t).ecg_data;
    end
    
    ecg_peaks = trials(t).ECG_spikes;
    signal_time = trials(t).time;
    ft_data_all.trial = [ft_data_all.trial, ...
        [signal; ecg_peaks]];
    ft_data_all.time = [ft_data_all.time, signal_time];

end

if ~isempty(ft_data_all.trial)

    % evoked LFP average
    cfg                 = [];
    cfg.keeptrials      = 'yes';
    cfg.timwin          = [cfg_ecg{3} cfg_ecg{4}];
    cfg.channel         = ft_data_all.label(1); % LFP chan
    cfg.spikechannel    = ft_data_all.label(2); % ECG peak
    
    ecg_based_sta       = ft_spiketriggeredaverage(cfg, ft_data_all);
    
    ecg_triggered_evoked.lfp = ecg_based_sta.trial;
    ecg_triggered_evoked.dimord = 'npeaks_time';
    ecg_triggered_evoked.lfp_time = ecg_based_sta.time;
    ecg_triggered_evoked.mean = ecg_based_sta.avg;
    ecg_triggered_evoked.std = permute(nanstd(ecg_based_sta.trial, 0, 1), [2 3 1]);
%     ecg_based_sta.std = nanstd(arr_state_lfp, 0, 1);

    clear ft_data_all ecg_based_sta;
    
end

end

