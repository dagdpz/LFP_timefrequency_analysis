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
lfp_tfa_cfg.use_datasets = [31];

% file list
lfp_tfa_cfg.file_list = ...
    {'Magnus', '20190130', 'Y:\Projects\PPC_pulv_body_signals\ephys\dPul_inactivation_20190130\sites_Magnus_20190130.mat';
    'Magnus', '20190213', 'Y:\Projects\PPC_pulv_body_signals\ephys\dPul_control_20190213\sites_Magnus_20190213.mat'; 
    'Magnus', '20190320', 'Y:\Projects\PPC_pulv_body_signals\ephys\MIP_control_20190320\sites_Magnus_20190320.mat'};
%     'Lin', '20170707', 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data\sites_Linus_20170707.mat';
%     'Lin', '20170713', 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data\sites_Linus_20170713.mat';
%     'Lin', '20170720', 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data\sites_Linus_20170720.mat';
%     'Lin', '20170802', 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data\sites_Linus_20170802.mat';
%     'Lin', '20170804', 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data\sites_Linus_20170804.mat';
%     'Lin', '20170818', 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data\sites_Linus_20170818.mat';
%     'Lin', '20171012', 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data\sites_Linus_20171012.mat'};

lfp_tfa_cfg.results_folder = 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data\LFP_TFA_Results';

% Specify events which mark trial start and end
lfp_tfa_cfg.trialinfo = struct();
lfp_tfa_cfg.trialinfo.start_state = lfp_tfa_states.FIX_ACQ;
lfp_tfa_cfg.trialinfo.ref_tstart = -0.4;
lfp_tfa_cfg.trialinfo.end_state = lfp_tfa_states.REWARD;
lfp_tfa_cfg.trialinfo.ref_tend = 0.6;

% reference hemisphere for hand-space labelling
lfp_tfa_cfg.ref_hemisphere = 'R';

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
lfp_tfa_cfg.tfr.timestep        = 25; % x sampling time

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
% whether to consider choice or instructed trials
lfp_tfa_cfg.baseline_use_choice_trial = 0; 

%% Settings for averaging TFR and evoked LFP based on conditions

% which type of trials (instructed / choice) trials to analyze
lfp_tfa_cfg.compare.types = [4];
lfp_tfa_cfg.compare.effectors = [6];
lfp_tfa_cfg.compare.targets = 'auto'; % 'auto' to automatically select
lfp_tfa_cfg.compare.choice_trials = 0; % 0 = only instructed, [0, 1] = both choice and instructed
lfp_tfa_cfg.compare.reach_hands = {'R', 'L'}; % for future use
lfp_tfa_cfg.compare.reach_spaces = {'R', 'L'}; % for future use
lfp_tfa_cfg.compare.perturbations = [0, 1]; % 0 = pre, 1 = post
lfp_tfa_cfg.compare.perturbation_groups = {0, 'all'}; % 'all', 'allbutone', 1xN int array

% define the peristates to analyse
% state id, state name, tbefore_onset(s), tafter_onset (s)
lfp_tfa_cfg.analyse_states = {6, 62};
%lfp_tfa_cfg.analyse_states = {6,    'Cue',      -1.0,   0.5;...
%                              62,   'Reach',    -0.5,   0.5};

% lfp_tfa_cfg.add_conditions = [];
% lfp_tfa_cfg.add_conditions(1).blocks = 'control'; % Analyse all inactivation blocks together
% lfp_tfa_cfg.add_conditions(2).blocks = 'inactivation'; % Analyse all inactivation blocks together
% Leave empty for block-wise analysis of all blocks

lfp_tfa_cfg.mintrials_percondition = 5;

% baseline configuration
cfg_baseline = [];
lfp_tfa_cfg.baseline_method = 'zscore';
    
%% Settings for average across sessions
lfp_tfa_cfg.compute_avg_across = {'sessions', 'sites'}; % 'sites'