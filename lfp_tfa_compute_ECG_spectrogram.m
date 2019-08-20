function [ session_ecg ] = lfp_tfa_compute_ECG_spectrogram( session_ecg, lfp_tfa_cfg )
%lfp_tfa_compute_site_tfr Summary of this function goes here
%   Detailed explanation goes here

% create an ft_datatype_raw
ft_data_ecg = struct();

%loop through each run
for b = (unique([session_ecg.trials.block]))
    fprintf('Computing TFR for run %g\n-----------------------\n', b);
    % concatenate all trials for this run
    concat_ECG = []; % to concatenate LFP data
    concat_time = []; % to concatenate sample time
    trials_time = [];
    trials_idx = find([session_ecg.trials.block] == b);
    for t = (trials_idx)
        %concat_trial_start_samples = [concat_trial_start_samples, length(concat_LFP) + 1];
        
        ts = session_ecg.trials(t).tsample;
        if ~isempty(concat_time) 
            trial_timestamps = ...
                linspace(ts, ts*length(session_ecg.trials(t).ecg_data), length(session_ecg.trials(t).ecg_data)) + concat_time(end);
        else
            trial_timestamps = ...
                linspace(0, ts*(length(session_ecg.trials(t).ecg_data)-1), length(session_ecg.trials(t).ecg_data));
        end
            
        concat_ECG = [concat_ECG, session_ecg.trials(t).ecg_data];
        concat_time = [concat_time trial_timestamps];
        trials_time = [trials_time; [trial_timestamps(1), trial_timestamps(end)]];
    end
    % get trial start times
    %concat_trial_start_times = concat_time(concat_trial_start_samples);
    
    ft_data_ecg.trial = {concat_ECG};
    ft_data_ecg.time = {concat_time};
    ft_data_ecg.fsample = 1/ts;
    ft_data_ecg.sampleinfo = [concat_time(1) concat_time(end)];
    ft_data_ecg.label = {[session_ecg.session]};
        
    % calculate the TFR
    cfg              = [];
    cfg.method       = lfp_tfa_cfg.tfr.method; % method used to calculate spectra             
    cfg.width        = lfp_tfa_cfg.tfr.width; % width of wavelet in number of cycles
    cfg.taper        = lfp_tfa_cfg.tfr.taper; % taper used
    cfg.pad          = 'nextpow2'; % zero padding for FFT
    cfg.foi          = lfp_tfa_cfg.tfr.foi;   % freqs of interest
    cfg.t_ftimwin    = lfp_tfa_cfg.tfr.t_ftimwin;           % length of sliding time window for each foi
    cfg.toi          = concat_time(1):lfp_tfa_cfg.tfr.timestep*ts:...
        concat_time(end);% time window "slides" from start time to and time in steps of timestep number of LFP samples
    cfg.channel      = ft_data_ecg.label;
    TFR_wavelet      = ft_freqanalysis(cfg, ft_data_ecg);
    
    % now divide into trials
    for t = (trials_idx)
        trial_tfr = TFR_wavelet;
        trial_tfr.powspctrm = single(TFR_wavelet.powspctrm(1, :, ...
            TFR_wavelet.time>=trials_time(trials_idx == t, 1) & ...
            TFR_wavelet.time<=trials_time(trials_idx == t, 2)));
        trial_tfr.time = TFR_wavelet.time(...
            TFR_wavelet.time>=trials_time(trials_idx == t, 1) & ...
            TFR_wavelet.time<=trials_time(trials_idx == t, 2)) ...
            - trials_time(trials_idx == t, 1) + session_ecg.trials(t).time(1);
        session_ecg.trials(t).tfs = trial_tfr;
    end
        
end



end
