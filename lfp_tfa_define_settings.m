function lfp_tfa_cfg = lfp_tfa_define_settings(version)
%lfp_tfa_define_settings - Function to define LFP TFA settings 
%   Detailed explanation goes here
    % configuration structure
    lfp_tfa_cfg = [];
    
    %% Settings for data folders
    % versioning
    lfp_tfa_cfg.version = version;

    % Absolute path to the folder containing LFP data to be analyzed
    lfp_tfa_cfg.data_folder = 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data';

    % folder to save results, automatically created
    lfp_tfa_cfg.root_results_folder = fullfile(lfp_tfa_cfg.data_folder, '\LFP_TFA_Results', date, ['ver' num2str(version)]);
    if ~exist(lfp_tfa_cfg.root_results_folder, 'dir')
        mkdir(lfp_tfa_cfg.root_results_folder);
    end

    % first read in the information about states
    all_states = lfp_tfa_define_states();
    lfp_tfa_cfg.all_states = all_states;

    % load epochs
    lfp_tfa_cfg.epochs = lfp_tfa_define_epochs();

    % Specify events which mark trial start and end
    lfp_tfa_cfg.trialinfo = struct();
    lfp_tfa_cfg.trialinfo.start_state = 'fxa';
    lfp_tfa_cfg.trialinfo.end_state = 'trh';
    
    % Specify lesional hemisphere
    lfp_tfa_cfg.lesional_hemispace = 'R';

    % maximum no:of sites to analyse from each session
    maxsites = 2; % inf = analyse all sites
    lfp_tfa_cfg.maxsites = maxsites;

    %% Settings for ft_freqanalysis in FieldTrip
    % Configuration for calculating LFP time frequency spectrogram using
    % ft_freqanalysis function of the fieldtrip toolbox
    lfp_tfa_cfg.tfr.method          = 'wavelet'; % 
    lfp_tfa_cfg.tfr.width           = 6; % no:of cycles
    lfp_tfa_cfg.tfr.foi             = logspace(log10(2), log10(120), 60);
    lfp_tfa_cfg.tfr.taper           = [];
    lfp_tfa_cfg.tfr.t_ftimwin       = [];
    lfp_tfa_cfg.tfr.tapsmofrq       = [];
    lfp_tfa_cfg.tfr.timestep        = 25; % x sampling time

    %% Settings to detect noisy trials
    % configuration for lfp noise rejection
    lfp_tfa_cfg.noise = [];
    % methods to be used
    lfp_tfa_cfg.noise.methods = {'amp', 'std', 'diff', 'pow'};
    % threshold for lfp raw amplitude (x std from mean)
    lfp_tfa_cfg.noise.amp_thr = 5;
    % number of consecutive samples beyond threshold to be considered
    lfp_tfa_cfg.noise.amp_N = 5;
    % no of standard deviations of trial w.r.t complete LFP std
    lfp_tfa_cfg.noise.std_thr = 4;
    % threshold for lfp derivative (x std from mean)
    lfp_tfa_cfg.noise.diff_thr = 4;
    % number of consecutive samples beyond threshold to be considered
    lfp_tfa_cfg.noise.diff_N = 5;
    % threshold for lfp power in number of standard deviations from mean
    lfp_tfa_cfg.noise.pow_thr = 4;
    % folder to save results
    lfp_tfa_cfg.noise.results_folder = root_results_folder;
    %cfg_noise.results_folder = [pathname '\Figures'];
    % whether single trials should be plotted
    lfp_tfa_cfg.noise.plottrials = 0;

    %% Settings to compute baseline power
    % reference state around which baseline should be considered, leave empty to consider complete trial
    lfp_tfa_cfg.baseline_ref_state = ''; 
    % period of interest for baseline calculation, trial = complete trial period
    lfp_tfa_cfg.baseline_period = 'trial'; 
    % which blocks to consider for baseline calculation
    lfp_tfa_cfg.baseline_block = 1; 
    % whether to consider choice or instructed trials
    lfp_tfa_cfg.use_choice_trial = 0; 
    %lfp_tfa_cfg.results_folder = root_results_folder;

    %% Settings for averaging TFR and evoked LFP based on conditions
    lfp_tfa_cfg.add_conditions = [];
    lfp_tfa_cfg.add_conditions(1).blocks = 'control'; % Analyse all inactivation blocks together
    lfp_tfa_cfg.add_conditions(2).blocks = 'inactivation'; % Analyse all inactivation blocks together
    % Leave empty for block-wise analysis of all blocks
    % define the peristates to analyse
    lfp_tfa_cfg.analyse_states = {6, 62};

    % baseline configuration
    cfg_baseline = [];
    lfp_tfa_cfg.baseline_method = 'zscore';
    
    %% Settings for average across sessions
    lfp_tfa_cfg.compute_avg_across = {'sessions', 'sites'}; % 'sites'
    
    %% save struct
    save(fullfile(lfp_tfa_cfg.root_results_fldr, ['settings_ver' num2str(version) '.mat']), ...
        'lfp_tfa_cfg');   
    

end

