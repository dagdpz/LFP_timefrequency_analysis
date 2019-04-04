function lfp_tfa_cfg = lfp_tfa_define_settings
%lfp_tfa_define_settings - Function to define LFP TFA settings 
%   Detailed explanation goes here
    % configuration structure
    lfp_tfa_cfg = [];    

    % first read in the information about states
    all_states = lfp_tfa_define_states();
    lfp_tfa_cfg.all_states = all_states;

    % load epochs
    lfp_tfa_cfg.epochs = lfp_tfa_define_epochs();

    % load LFP data for the selected session
    load(fullfile(lfp_tfa_cfg.data_folder, session_filename));

    lfp_tfa_cfg.trialinfo = struct();
    lfp_tfa_cfg.trialinfo.start_state = 'fxa';
    lfp_tfa_cfg.trialinfo.end_state = 'trh';

    lfp_tfa_cfg.lesional_hemispace = 'R';

    % maximum no:of sites to analyse
    maxsites = inf; % inf = analyse all sites
    lfp_tfa_cfg.maxsites = maxsites;

    %% Read the required fields from  the processed LFP data for the session
    % Configuration for calculating LFP time frequency spectrogram using
    % ft_freqanalysis function of the fieldtrip toolbox
    lfp_tfa_cfg.tfr.method          = 'wavelet'; % 
    lfp_tfa_cfg.tfr.width           = 6; % no:of cycles
    lfp_tfa_cfg.tfr.foi             = logspace(log10(2), log10(120), 60);
    lfp_tfa_cfg.tfr.taper           = [];
    lfp_tfa_cfg.tfr.t_ftimwin       = [];
    lfp_tfa_cfg.tfr.tapsmofrq       = [];
    lfp_tfa_cfg.tfr.timestep        = 25; % x sampling time

    %% Reject noisy trials
    % configuration for lfp noise rejection
    lfp_tfa_cfg.noise = [];
    % methods to be used
    lfp_tfa_cfg.noise.methods = {'amp', 'std', 'diff', 'pow'};
    % threshold for lfp raw amplitude
    lfp_tfa_cfg.noise.amp_thr = 5;
    % number of consecutive samples beyond threshold to be considered
    lfp_tfa_cfg.noise.amp_N = 5;
    % no of standard deviations of trial w.r.t complete LFP std
    lfp_tfa_cfg.noise.std_thr = 4;
    % threshold for lfp derivative in percentile
    lfp_tfa_cfg.noise.diff_thr = 4;
    % number of consecutive samples beyond threshold to be considered
    lfp_tfa_cfg.noise.diff_N = 5;
    % threshold for lfp power in standard deviations
    lfp_tfa_cfg.noise.pow_thr = 4;
    % folder to save results
    lfp_tfa_cfg.noise.results_folder = root_results_folder;
    %cfg_noise.results_folder = [pathname '\Figures'];
    % whether single trials should be plotted
    lfp_tfa_cfg.noise.plottrials = 0;

    %filt_session_lfp = rejectNoisyLFPTrials( session_lfp )

    %% Compute baseline
    lfp_tfa_cfg.baseline_ref_state = ''; % reference state around which baseline should be considered, leave empty to consider complete trial
    lfp_tfa_cfg.baseline_period = 'trial'; % period of interest for baseline calculation, trial = complete trial period
    lfp_tfa_cfg.baseline_block = 1; % consider only block 1 (control) for baseline calculation
    lfp_tfa_cfg.use_choice_trial = 0; % whether to consider choice or instructed trials
    lfp_tfa_cfg.results_folder = root_results_folder;

    %% Compute the TFR per site and average across sites
    lfp_tfa_cfg.add_conditions = [];
    lfp_tfa_cfg.add_conditions(1).blocks = 'inactivation'; % Analyse all inactivation blocks together
    % Leave empty for block-wise analysis of all blocks
    % define the peristates to analyse
    lfp_tfa_cfg.analyse_states = {6, 62};

    % baseline configuration
    cfg_baseline = [];
    lfp_tfa_cfg.baseline_method = 'zscore';   
    

end

