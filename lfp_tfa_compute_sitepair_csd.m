function sitepair_crosspow = lfp_tfa_compute_sitepair_csd(site1_lfp, site2_lfp, lfp_tfa_cfg)
% lfp_tfa_compute_sitepair_csd - Computes the LFP-LFP cross power spectrogram
% between a pair of sites for all trials of a session. 
% This function calls the ft_freqanalysis routine
% of Fieldtrip toolbox for the calculation of LFP-LFP cross spectrum
%
% USAGE:
%	sitepair_crosspow = lfp_tfa_compute_sitepair_csd(site1_lfp, ...
%       site2_lfp, lfp_tfa_cfg)
%
% INPUTS:
%		site1_lfp        - struct containing LFP signal for each trial for
%		one site, see lfp_tfa_process_LFP
%       site2_lfp        - struct containing LFP signal for each trial for
%		another site, see lfp_tfa_process_LFP
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
%		sitepair_crosspow       - structure containing trial-wise LFP-LFP
%		cross power spectrograms
%
% See also lfp_tfa_process_lfp, settings/lfp_tfa_settings_example,
% lfp_tfa_sitepair_averaged_sync
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
    
% struct to store information about sitepair
sitepair_crosspow = struct();
sitepair_crosspow.session = site1_lfp.session;
sitepair_crosspow.sites = {site1_lfp.site_ID, site2_lfp.site_ID};
sitepair_crosspow.targets = {site1_lfp.target, site2_lfp.target};
% sitepair_crosspow.site_pos = {[site1_lfp.xpos, site1_lfp.ypos, site1_lfp.zpos], ...
%     [site2_lfp.xpos, site2_lfp.ypos, site2_lfp.zpos]};
% sitepair_crosspow.site_dist = norm(sitepair_crosspow.site_pos{2} - ...
%     sitepair_crosspow.site_pos{1});
sitepair_crosspow.trials = site1_lfp.trials;
sitepair_crosspow.trials = rmfield(sitepair_crosspow.trials, 'lfp_data');
sitepair_crosspow.trials = rmfield(sitepair_crosspow.trials, 'time');
%sitepair_crosspow.ntrials = site1_lfp.ntrails;

% concat_lfppair = []; % to concatenate LFP data
% concat_time = []; % to concatenate sample time
% trials_time = [];

%loop through each run
for r = (unique([site1_lfp.trials.run]))
    fprintf('Run %g\n -----------------\n', r);
    trials_idx = find([site1_lfp.trials.run] == r);
    concat_lfppair = []; % to concatenate LFP data
    concat_time = []; % to concatenate sample time
    trials_time = [];
    for t = (trials_idx)
        %concat_trial_start_samples = [concat_trial_start_samples, length(concat_LFP) + 1];
        
        ts = site1_lfp.trials(t).tsample;
        if ~isempty(concat_time) 
            trial_timestamps = ...
                linspace(ts, ts*length(site1_lfp.trials(t).lfp_data), length(site1_lfp.trials(t).lfp_data)) + concat_time(end);
        else
            trial_timestamps = ...
                linspace(0, ts*(length(site1_lfp.trials(t).lfp_data)-1), length(site1_lfp.trials(t).lfp_data));
        end
            
        concat_lfppair = [concat_lfppair, [site1_lfp.trials(t).lfp_data; site2_lfp.trials(t).lfp_data]];
        concat_time = [concat_time trial_timestamps];
        trials_time = [trials_time; [trial_timestamps(1), trial_timestamps(end)]];
    end
    % get trial start times
    %concat_trial_start_times = concat_time(concat_trial_start_samples);
    
    ft_data_lfp.trial       = {concat_lfppair};
    ft_data_lfp.time        = {concat_time};
    ft_data_lfp.fsample     = 1/ts;
    ft_data_lfp.sampleinfo  = [concat_time(1) concat_time(end)];
    ft_data_lfp.label       = {site1_lfp.site_ID, site2_lfp.site_ID};

    timestamps = ft_data_lfp.time{1};
    % Calculate the cross power spectrum
    cfg                     = [];
    cfg.method              = lfp_tfa_cfg.tfr.method;
    cfg.width               = lfp_tfa_cfg.tfr.width;
    cfg.output              = 'powandcsd';
    cfg.taper               = lfp_tfa_cfg.tfr.taper;
    cfg.pad                 = 'nextpow2';
    cfg.foi                 = lfp_tfa_cfg.tfr.foi;   % frequencies of interest
    %cfg.t_ftimwin          = ones(length(cfg.foi),1).*0.5; %lfp_tfa_cfg.tfr.t_ftimwin;           % number of cycles per time window
    cfg.toi                 = timestamps(1):lfp_tfa_cfg.tfr.timestep*(1/ft_data_lfp.fsample):...
        timestamps(end);% time period of interest
    %cfg.keeptrials         = 'yes';
    TFR_wavelet             = ft_freqanalysis(cfg, ft_data_lfp); 
        
    % now divide into trials
    for t = (trials_idx)
        % save cross spectrum for trial
        trial_tfr = TFR_wavelet;
        trial_tfr.powspctrm = single(TFR_wavelet.powspctrm(:, :, ...
            TFR_wavelet.time>=trials_time(trials_idx == t, 1) & ...
            TFR_wavelet.time<=trials_time(trials_idx == t, 2)));
        trial_tfr.crsspctrm = single(TFR_wavelet.crsspctrm(:, :, ...
            TFR_wavelet.time>=trials_time(trials_idx == t, 1) & ...
            TFR_wavelet.time<=trials_time(trials_idx == t, 2)));
        trial_tfr.time = single(TFR_wavelet.time(...
            TFR_wavelet.time>=trials_time(trials_idx == t, 1) & ...
            TFR_wavelet.time<=trials_time(trials_idx == t, 2)) ...
            - trials_time(trials_idx == t, 1) + site1_lfp.trials(t).time(1));
        sitepair_crosspow.trials(t).csd = trial_tfr; 
        sitepair_crosspow.trials(t).noisy = site1_lfp.trials(t).noisy | ...
            site2_lfp.trials(t).noisy;
    end
    
end

