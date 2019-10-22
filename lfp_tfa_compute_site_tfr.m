function [ site_lfp ] = lfp_tfa_compute_site_tfr( site_lfp, lfp_tfa_cfg )
% lfp_tfa_compute_site_tfr - Computes the LFP time frequency spectrogram
% for all trials of a site. This function calls the ft_freqanalysis routine
% of Fieldtrip toolbox for the calculation of LFP power spectrogram
%
% USAGE:
%	site_lfp = lfp_tfa_compute_site_tfr( site_lfp, lfp_tfa_cfg )
%
% INPUTS:
%		site_lfp        - struct containing LFP signal for each trial for a
%		single site
%       lfp_tfa_cfg      - settings for TFR computation, see 
%       settings/lfp_tfa_settings_example
%           Required Fields: 
%               1. tfr.method   - method to be used for calculating the LFP
%               power spectra. Can be 'mtmconvol' or 'wavelet'
%               2. tfr.width    - width of the wavelets in number of
%               cycles (For method = 'wavelet', Ignored for 'mtmconvol')
%               3. tfr.taper    - taper (single or multiple) to be used.
%               Can be 'dpss', 'hanning' or many others (Used when
%               tfr.method = 'mtmconvol', ignored for 'wavelet'
%               4. tfr.foi      - frequencies of interest (in Hz), a vector
%               of freq values or a 1x2 vector with start and end freq
%               5. tfr.t_ftimwin - length of the sliding time-window in 
%               seconds, should be vector of length 1 x numfoi (Only used
%               when tfr.method = 'mtmconvol')
%               6. tfr.timestep - number of lfp samples to step for the 
%               sliding time window
% OUTPUTS:
%		site_lfp         - structure containing trial-wise LFP power
%		spectrograms
%
% See also lfp_tfa_process_lfp, settings/lfp_tfa_settings_example
% 
% Author(s):	S.Nair, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2019-02-15:	Created function (Sarath Nair)
% 2019-03-05:	First Revision
% ...
% $Revision: 1.0 $  $Date: 2019-03-05 17:18:00 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%

% create an ft_datatype_raw
ft_data_lfp = struct();

%loop through each run
for r = (unique([site_lfp.trials.run]))
    fprintf('Computing TFR for run %g\n-----------------------\n', r);
    % concatenate all trials for this run
    concat_LFP = []; % to concatenate LFP data
    concat_time = []; % to concatenate sample time
    trials_time = [];
    trials_idx = find([site_lfp.trials.run] == r);
    for t = (trials_idx)
        %concat_trial_start_samples = [concat_trial_start_samples, length(concat_LFP) + 1];
        
        ts = site_lfp.trials(t).tsample;
        if ~isempty(concat_time) 
            trial_timestamps = ...
                linspace(ts, ts*length(site_lfp.trials(t).lfp_data), length(site_lfp.trials(t).lfp_data)) + concat_time(end);
        else
            trial_timestamps = ...
                linspace(0, ts*(length(site_lfp.trials(t).lfp_data)-1), length(site_lfp.trials(t).lfp_data));
        end
            
        concat_LFP = [concat_LFP, site_lfp.trials(t).lfp_data];
        concat_time = [concat_time trial_timestamps];
        trials_time = [trials_time; [trial_timestamps(1), trial_timestamps(end)]];
    end
    % get trial start times
    %concat_trial_start_times = concat_time(concat_trial_start_samples);
    
    ft_data_lfp.trial = {concat_LFP};
    ft_data_lfp.time = {concat_time};
    ft_data_lfp.fsample = 1/ts;
    ft_data_lfp.sampleinfo = [concat_time(1) concat_time(end)];
    ft_data_lfp.label = {[site_lfp.site_ID]};
        
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
    cfg.channel      = ft_data_lfp.label;
    TFR_wavelet      = ft_freqanalysis(cfg, ft_data_lfp);
    
    % now divide into trials
    for t = (trials_idx)
        trial_tfr = TFR_wavelet;
        trial_tfr.powspctrm = single(TFR_wavelet.powspctrm(1, :, ...
            TFR_wavelet.time>=trials_time(trials_idx == t, 1) & ...
            TFR_wavelet.time<=trials_time(trials_idx == t, 2)));
        trial_tfr.time = TFR_wavelet.time(...
            TFR_wavelet.time>=trials_time(trials_idx == t, 1) & ...
            TFR_wavelet.time<=trials_time(trials_idx == t, 2)) ...
            - trials_time(trials_idx == t, 1) + site_lfp.trials(t).time(1);
        site_lfp.trials(t).tfs = trial_tfr;
    end
        
end



end
