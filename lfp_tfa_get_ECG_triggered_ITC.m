function ecg_triggered_itc = lfp_tfa_get_ECG_triggered_ITC( site_lfp, cond_trials, state, lfp_tfa_cfg )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

state_name = state{2};
epoch_ref_tstart = state{3};
epoch_ref_tend = state{4};
width = state{4} - state{3};

ecg_triggered_itc.fourierspctrm = {}; % fourier
ecg_triggered_itc.time = {}; % timebins fo spectrogram
ecg_triggered_itc.freq = {}; % freq bins
ecg_triggered_itc.state = [];
ecg_triggered_itc.state_name = state_name;

for t = find(cond_trials)

    trialperiod           = site_lfp.trials(t).trialperiod;
    
    % get the LFP samples and timestamps between start and end states
    lfp_tfs_fourier = site_lfp.trials(t).tfs.fourierspctrm;
    lfp_tfs_time = site_lfp.trials(t).tfs.time;
    lfp_tfs_freq = site_lfp.trials(t).tfs.freq;
    
    % ecg peak times
    lfp_time = site_lfp.trials(t).time; 
    ecg_peaks = site_lfp.trials(t).ECG_spikes; %...
%         ((site_lfp.trials(t).time >= trialperiod(1) & ...
%         site_lfp.trials(t).time <= trialperiod(2))); 
    ecg_peak_times = lfp_time(ecg_peaks);
    
    % lfp sample indices
    lfp_timebin_idx = 1:length(lfp_tfs_time);
    % width of each time bin
    lfp_ts = 1/site_lfp.trials(t).fsample;
    tbin_width = lfp_tfa_cfg.tfr.timestep * lfp_ts;
    
    % number fo time bins in the window
    ntbins_window = round(width/tbin_width);
    
    ecg_triggered_itc.cfg = site_lfp.trials(t).tfs.cfg;
    
    % now get the midpoint of windows around ECG peak
    window_mid_idx = [];
    for p = ecg_peak_times
        if min(abs(lfp_tfs_time - p)) > tbin_width
            continue;
        end
        window_mid_idx = [window_mid_idx, find(abs(lfp_tfs_time - p) == ...
            min(abs(lfp_tfs_time - p)))];
    end

    % loop through each window
    for w = 1:length(window_mid_idx)
        if window_mid_idx(w) + round(epoch_ref_tstart/tbin_width) < 1 || ...
                window_mid_idx(w) + round(epoch_ref_tend/tbin_width) > length(lfp_tfs_time)
            continue;
        end
        % power spectrum for this window
        ecg_triggered_itc.fourierspctrm = [ecg_triggered_itc.fourierspctrm, ...
            lfp_tfs_fourier(:, :, window_mid_idx(w) + round(epoch_ref_tstart/tbin_width):...
             window_mid_idx(w) + round(epoch_ref_tend/tbin_width))];
        % timestamps
        ecg_triggered_itc.time = ...
            lfp_tfs_time(window_mid_idx(w) + round(epoch_ref_tstart/tbin_width):...
             window_mid_idx(w) + round(epoch_ref_tend/tbin_width));
        % set mid-timestamp to zero
        ecg_triggered_itc.time = ecg_triggered_itc.time - ...
            lfp_tfs_time(window_mid_idx(w));

        ecg_triggered_itc.freq = lfp_tfs_freq;
    end
       
end


% find number of time bins in power
% spectrogram
ntimebins = min(cellfun('size', ecg_triggered_itc.fourierspctrm, 3));
nfreqbins = numel(ecg_triggered_itc.freq);
% crop each tfs to the ntimebins
for k = 1:length(ecg_triggered_itc.fourierspctrm)
    ecg_triggered_itc.fourierspctrm{k} = ecg_triggered_itc.fourierspctrm{k}(1,:,1:ntimebins);
    
end
ecg_triggered_itc.time = ecg_triggered_itc.time(1:ntimebins);

if ~isempty(ecg_triggered_itc.fourierspctrm)

    % find the average TFS for each state
    arr_state_fourier = cat(1, ecg_triggered_itc.fourierspctrm{:});
    nRpeak_epochs = size(arr_state_fourier, 1);
    
    % compute inter-trial phase coherence (itpc)
    ecg_triggered_itc.itpc      = arr_state_fourier./abs(arr_state_fourier); % divide by amplitude
    ecg_triggered_itc.itpc      = nansum(ecg_triggered_itc.itpc,1);   % sum angles
    ecg_triggered_itc.itpc      = abs(ecg_triggered_itc.itpc)/nRpeak_epochs;   % take the absolute value and normalize
    ecg_triggered_itc.itpc      = squeeze(ecg_triggered_itc.itpc); % remove the first singleton dimension
    
    % compute inter-trial linear coherence (itlc)
    ecg_triggered_itc.itlc      = nansum(arr_state_fourier) ./ ...
        (sqrt(nRpeak_epochs*nansum(abs(arr_state_fourier).^2)));
    ecg_triggered_itc.itlc      = abs(ecg_triggered_itc.itlc);     % take the absolute value, i.e. ignore phase
    ecg_triggered_itc.itlc      = squeeze(ecg_triggered_itc.itlc); % remove the first singleton dimension
    
    
end

end

