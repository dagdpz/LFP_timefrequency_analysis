%% Initialization

% call global states
%lfp_tfa_global_define_states;

% initialize configuration structure
lfp_tfa_cfg = [];
   
%% Settings for data folders
% versioning
lfp_tfa_cfg.version = 'SN_0.2';

% Absolute path to the folder containing LFP data to be analyzed
%lfp_tfa_cfg.data_folder = 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data';

% sorted neurons excel file
lfp_tfa_cfg.info_filepath = 'Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn\Lin_sorted_neurons.xls';

% dataset to be used
lfp_tfa_cfg.use_datasets = [31]; %???

% file list format ???
lfp_tfa_cfg.file_list = ...
    {'Magnus', '20190124', '';
    'Magnus', '20190314', ''};

lfp_tfa_cfg.results_folder = 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data\LFP_TFA_Results';

% Specify events which mark trial start and end
lfp_tfa_cfg.trialinfo = struct();
lfp_tfa_cfg.trialinfo.start_state = lfp_tfa_states.FIX_ACQ;
lfp_tfa_cfg.trialinfo.ref_tstart = -0.4; % s, start_state + ref_tstart will be taken
lfp_tfa_cfg.trialinfo.end_state = lfp_tfa_states.REWARD;
lfp_tfa_cfg.trialinfo.ref_tend = 0.6;

% reference hemisphere for hand-space labelling
lfp_tfa_cfg.ref_hemisphere = 'R'; % can be 'R' or 'L'

% maximum no:of sites to analyse from each session
maxsites = inf; % inf = analyse all sites
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
lfp_tfa_cfg.tfr.timestep        = 25; % number of lfp samples for one step 

%% Settings to detect noisy trials
% configuration for lfp noise rejection
lfp_tfa_cfg.noise = [];
% methods to be used
lfp_tfa_cfg.noise.methods = {'amp', 'std', 'diff', 'pow'};
% threshold for lfp raw amplitude (x std from mean)
lfp_tfa_cfg.noise.amp_thr = 6;
% number of consecutive samples beyond threshold to be considered
lfp_tfa_cfg.noise.amp_N = 10;
% no of standard deviations of trial w.r.t complete LFP std
lfp_tfa_cfg.noise.std_thr = 4;
% threshold for lfp derivative (x std from mean)
lfp_tfa_cfg.noise.diff_thr = 6;
% number of consecutive samples beyond threshold to be considered
lfp_tfa_cfg.noise.diff_N = 10;
% threshold for lfp power in number of standard deviations from mean
lfp_tfa_cfg.noise.pow_thr = 4;
%cfg_noise.results_folder = [pathname '\Figures'];
% whether single trials should be plotted
lfp_tfa_cfg.noise.plottrials = 0;

%% Settings to compute baseline power
% reference state around which baseline should be considered, leave empty to consider complete trial
lfp_tfa_cfg.baseline_ref_state = ''; 
% period of interest for baseline calculation, trial = complete trial period
lfp_tfa_cfg.baseline_ref_period = 'trial'; 
% which blocks to consider for baseline calculation
lfp_tfa_cfg.baseline_perturbation = 0; 
% whether to consider choice (1) or instructed trials (0)
lfp_tfa_cfg.baseline_use_choice_trial = 0; 

%% Settings for averaging TFR and evoked LFP based on conditions

% which type of trials (instructed / choice) trials to analyze
lfp_tfa_cfg.compare.targets = {'MIPa_R', 'MIPa_L'}; 
lfp_tfa_cfg.compare.types = [4];
lfp_tfa_cfg.compare.effectors = [6];
lfp_tfa_cfg.compare.choice_trials = 0; % 0 = only instructed, [0, 1] = both choice and instructed
lfp_tfa_cfg.compare.reach_hands = {'R', 'L'}; % for future use
lfp_tfa_cfg.compare.reach_spaces = {'R', 'L'}; % for future use
lfp_tfa_cfg.compare.perturbations = [0, 1]; % 0 = pre, 1 = post
lfp_tfa_cfg.compare.perturbation_groups = {0, 'all'}; % 'all', 'allbutone', 1xN int array

% define the states to analyse for LFP TFR and evoked LFP response
%lfp_tfa_cfg.analyse_states = {6, 62};
% {state_id, state_name, ref_tstart, ref_tend}
lfp_tfa_cfg.analyse_states = {6,    'Cue',      -1.0,   0.5;...
                             62,   'Reach',    -0.5,   0.5};

% define the states to analyse for LFP power spectrum
%lfp_tfa_cfg.analyse_states = {6, 62};
% {state_id, state_name, ref_tstart, ref_tend}
lfp_tfa_cfg.analyse_epochs = {6,    'FHol',    -0.3 ,    0  ;...
                              6,    'Cue' ,    0.05 ,    0.2 ; ...
                              8,    'EDel',    0.3 ,     0.6 ; ...
                              4,    'Del',     -0.3 ,    0  ; ...
                              62,   'PreR',    -0.3 ,    -0.05 ; ...
                              63,   'PeriR',   -0.2 ,    0.2 ; ...
                              20,   'THol',    -0.3 ,    0    };

lfp_tfa_cfg.mintrials_percondition = 5;

% baseline configuration
cfg_baseline = [];
lfp_tfa_cfg.baseline_method = 'zscore';
    
%% Settings for average across sessions
lfp_tfa_cfg.compute_avg_across = {'sessions', 'sites'}; % 'sites'