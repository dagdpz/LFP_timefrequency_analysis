%% Initialization

% initialize configuration structure
lfp_tfa_cfg = [];
   
%% Settings for data folders
% versioning, a unique version for the settings file
% the results produced using this settings file would be saved under 
% the folder [lfp_tfa_cfg.results_folder, '\', date, '\ver_' lfp_tfa_cfg.version]
% eg: 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data\LFP_TFA_Results\20190506\ver_SN_0.2'
lfp_tfa_cfg.version = 'SN_0.2';

% sorted neurons excel file, from which information about sessions and
% individual sites can be obtained
lfp_tfa_cfg.info_filepath = 'Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn\Lin_sorted_neurons.xls';

% dataset to be used for analysis, see entry 'Set' in the sorted neurons excel file
% only those sessions belonging to 'Set' = lfp_tfa_cfg.use_datasets will be
% used for analysis
lfp_tfa_cfg.use_datasets = [31];

% list of sessions to be analysed
% { Monkey name, Date of session, Absolute path to the file containing LFP data for
% the session }
lfp_tfa_cfg.file_list = ...
    {'Magnus', '20190124', 'Y:\Projects\PPC_pulv_body_signals\ephys\MIP_inactivation_20190124\sites_Magnus_20190124.mat';
    'Magnus', '20190314', 'Y:\Projects\PPC_pulv_body_signals\ephys\MIP_inactivation_20190314\sites_Magnus_20190314.mat'};

% absolute path to the folder where the results of analysis should be stored
lfp_tfa_cfg.results_folder = 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data\LFP_TFA_Results';

% Specify events which mark trial start and end
lfp_tfa_cfg.trialinfo = struct();

% ID of the reference state which indicates start of a trial
% Example:
% lfp_tfa_cfg.trialinfo.start_state = lfp_tfa_states.FIX_ACQ; reference for 
% trial start is the onset of fixation acquisition
lfp_tfa_cfg.trialinfo.start_state = lfp_tfa_states.FIX_ACQ;

% offset to be considered from the onset of
% trial start reference state for calculating the trial start time
% i.e., trial start time = onset of trial start state + start offset
% Example:
% 1. lfp_tfa_cfg.trialinfo.ref_tstart = -0.5;
% trial start time = onset time of lfp_tfa_cfg.trialinfo.start_state - 0.5;
% 1. lfp_tfa_cfg.trialinfo.ref_tstart = 0.5;
% trial start time = onset time of lfp_tfa_cfg.trialinfo.start_state + 0.5;
lfp_tfa_cfg.trialinfo.ref_tstart = -0;

% ID of the reference state which indicates start of a trial
% Example:
% lfp_tfa_cfg.trialinfo.end_state = lfp_tfa_states.TAR_HOL; reference for 
% trial start is the onset of target hold
lfp_tfa_cfg.trialinfo.end_state = lfp_tfa_states.TAR_HOL;

% offset to be considered from the onset of
% trial end reference state for calculating the trial end time
% i.e., trial end time = onset of trial end state + end offset
% Example:
% 1. lfp_tfa_cfg.trialinfo.ref_tend = 0.5;
% trial start time = onset time of lfp_tfa_cfg.trialinfo.end_state + 0.5;
% 1. lfp_tfa_cfg.trialinfo.ref_tend = -0.5;
% trial start time = onset time of lfp_tfa_cfg.trialinfo.end_state - 0.5;
lfp_tfa_cfg.trialinfo.ref_tend = 0;

% reference hemisphere for hand-space labelling
% can be 'R' (for right hemisphere) or 'L' (for left hemisphere)
% ref_hemisphere is used for labelling contra and ipsi hand and space
% set ref_hemisphere to lesioned hemishere for ipsi lesional and contra
% lesional labeling
% set ref_hemisphere to recorded hemishere for ipsi lateral and contra
% lateral labeling
lfp_tfa_cfg.ref_hemisphere = 'R'; 

% maximum no:of sites to analyse from each session
% If maxsites < total number of sites in a session, only maxsite number of
% sites will be analysed
% Examples:
% 1. lfp_tfa_cfg.maxsites = 2; only first two sites will be analysed from 
% each session
% 1. lfp_tfa_cfg.maxsites = inf; all the sites will be analysed from 
% each session
lfp_tfa_cfg.maxsites = inf; % inf = analyse all sites

%% Settings for ft_freqanalysis in FieldTrip
% Configuration for calculating LFP time frequency spectrogram using
% ft_freqanalysis function of the fieldtrip toolbox
lfp_tfa_cfg.tfr.method          = 'wavelet'; % 
lfp_tfa_cfg.tfr.width           = 6; % no:of cycles
lfp_tfa_cfg.tfr.foi             = logspace(log10(2), log10(120), 60);
lfp_tfa_cfg.tfr.taper           = [];
lfp_tfa_cfg.tfr.t_ftimwin       = [];
lfp_tfa_cfg.tfr.tapsmofrq       = [];
lfp_tfa_cfg.tfr.timestep        = 25; % number of lfp samples to step

%% Settings to detect noisy trials
% configuration for lfp noise rejection
lfp_tfa_cfg.noise = [];
% whether or not to apply noise rejection - future use
lfp_tfa_cfg.noise.detect = 0;
% combination of methods to be used - future use
% currently all methods are used together 
lfp_tfa_cfg.noise.methods = {'amp', 'std', 'diff', 'pow'};
% threshold for lfp raw amplitude (number of std deviations from mean)
lfp_tfa_cfg.noise.amp_thr = 6;
% number of consecutive samples beyond threshold to be considered for marking 
% a noisy trial
lfp_tfa_cfg.noise.amp_N = 10;
% no of standard deviations of trial LFP w.r.t LFP std of all trials
lfp_tfa_cfg.noise.std_thr = 4;
% threshold for lfp derivative (number of std deviations from mean)
lfp_tfa_cfg.noise.diff_thr = 6;
% number of consecutive samples beyond threshold to be considered for marking 
% a noisy trial
lfp_tfa_cfg.noise.diff_N = 10;
% threshold for lfp power in number of standard deviations from mean
lfp_tfa_cfg.noise.pow_thr = 4;
% whether single trials should be plotted
lfp_tfa_cfg.noise.plottrials = 0;

%% Settings to compute baseline power

% ID of the reference state around which baseline should be considered, see
% lfp_tfa_global_define_states
% Examples
% 1. lfp_tfa_cfg.baseline_ref_state = lfp_tfa_states.CUE_ON; considers cue
% onset as the reference state for baseline period
% 2. lfp_tfa_cfg.baseline_ref_state = ''; consider the whole trial period
% for baseline
lfp_tfa_cfg.baseline_ref_state = ''; 

% period of interest relative to onset of baseline_ref_state for baseline power calculation, 
% Examples: 
% 1. lfp_tfa_cfg.baseline_ref_period = [-0.5 -0.1]; considers the time 
% period -0.5 s to -0.1 s from the onset of baseline_ref_state (i.e., if 
% lfp_tfa_cfg.baseline_ref_state = lfp_tfa_states.CUE_ON, time period from 
% -0.5 s to -0.1s from the cue onset is considered as baseline period)
% 2. lfp_tfa_cfg.baseline_ref_period = 'trial'; to consider the complete trial period
% for baseline power calculation
lfp_tfa_cfg.baseline_ref_period = 'trial'; 

% which perturbation blocks to be considered for baseline power calculation
% set to 0 for considering only pre-injection blocks
% Examples: 
% 1. lfp_tfa_cfg.baseline_perturbation = 0; only perturbation block 0
% (pre-injection) is used for baseline power calculation
% 2. lfp_tfa_cfg.baseline_perturbation = [0, 2]; combines perturbation blocks
% 0 (pre-injection) and 2 (post-injection), but is not recommended
% 3. lfp_tfa_cfg.baseline_perturbation = [2, 3]; combines perturbation blocks
% 2 and 3
lfp_tfa_cfg.baseline_perturbation = 0; 

% whether to consider choice (1) or instructed trials (0) in baseline power
% calculation 
% Examples:
% 1. lfp_tfa_cfg.baseline_use_choice_trial = 0; % consider only instructed trials
% 2. lfp_tfa_cfg.baseline_use_choice_trial = 1; % consider only choice trials
lfp_tfa_cfg.baseline_use_choice_trial = 0; 

%% Settings for averaging TFR and evoked LFP based on conditions

% targets to be included in the analysis
% should be a cell array of strings which indicate the target names
% the target names should be same as the target field in the LFP data
% structure
% Those targets which are not in the analysed sessions will be ignored
% Example:
% 1. lfp_tfa_cfg.compare.targets = {'MIPa_R', 'MIPa_L', 'dPul_R', 'dPul_L'}; 
lfp_tfa_cfg.compare.targets = {'MIPa_R', 'MIPa_L'}; 

% trial types to be included in the analysis
% should be a vector of integers specifying the types
% entries should be same as one of the type specified in the 'type' field
% in the input LFP data
% Example: 
% 1. lfp_tfa_cfg.compare.types = [4, 2]; % analyse trials with type = 4
% and type = 2 separately
% 2. lfp_tfa_cfg.compare.types = nan; Ignore trial type (trials with any
% type value are combined)
lfp_tfa_cfg.compare.types = [4];

% effectors to be included in the analysis
% should be a vector of integers specifying the effectors
% entries should be same as one of the type specified in the 'effector' field
% in the input LFP data
% Example: 
% 1. lfp_tfa_cfg.compare.types = [4, 6]; % analyse trials with effector = 4
% and effector = 6 separately
% 2. lfp_tfa_cfg.compare.types = nan; Ignore effector (trials with any
% effector value are combined)
lfp_tfa_cfg.compare.effectors = [6];

% which type of choice trials are to be included in the analysis
% Examples:
% 1. lfp_tfa_cfg.compare.choice_trials = 0; % analyse only instructed trials
% 2. lfp_tfa_cfg.compare.choice_trials = 1; % analyse only choice trials
% 3. lfp_tfa_cfg.compare.choice_trials = [0, 1]; % analyse choice and 
% instructed trials separately
% 3. lfp_tfa_cfg.compare.choice_trials = nan; % ignore choice (both choice
% and instructed trials are combined)
lfp_tfa_cfg.compare.choice_trials = 0; 

% reach hands to be included for analysis
% should be nan or a cell array that contain only values 'R', 'L'
% Examples:
% 1. lfp_tfa_cfg.compare.reach_hands = {'L'}; include only those trials in
% which reach hand is left
% 2. lfp_tfa_cfg.compare.reach_hands = {'R'}; include only those trials in
% which reach hand is right
% 3. lfp_tfa_cfg.compare.reach_hands = {'L', 'R'}; analyse the trials in
% which reach hand is left and right separately
% 4. lfp_tfa_cfg.compare.reach_hands = nan; ignore hand label (trial with
% any hand label is combined)
lfp_tfa_cfg.compare.reach_hands = {'R', 'L'}; % for future use

% reach space to be included for analysis
% should be nan or a cell array that contain only values 'R', 'L'
% Examples:
% 1. lfp_tfa_cfg.compare.reach_hands = {'L'}; include only those trials in
% which reach hand is left
% 2. lfp_tfa_cfg.compare.reach_hands = {'R'}; include only those trials in
% which reach hand is right
% 3. lfp_tfa_cfg.compare.reach_hands = {'L', 'R'}; analyse the trials in
% which reach hand is left and right separately
% 4. lfp_tfa_cfg.compare.reach_hands = nan; ignore hand label (trial with
% any hand label is combined)
lfp_tfa_cfg.compare.reach_spaces = {'R', 'L'}; 

% perturbations to be included in the analysis
% should be nan, 0, 1 or [0, 1]
% Examples:
% lfp_tfa_cfg.compare.perturbations = 0; consider only those trials with
% perturbation value = lfp_tfa_cfg.compare.perturbation_groups(0)
% lfp_tfa_cfg.compare.perturbations = 1; consider only those trials with
% perturbation value = lfp_tfa_cfg.compare.perturbation_groups(1)
% lfp_tfa_cfg.compare.perturbations = [0, 1]; consider the trials with
% perturbation value = lfp_tfa_cfg.compare.perturbation_groups(0) or 
% lfp_tfa_cfg.compare.perturbation_groups(1) separately
% lfp_tfa_cfg.compare.perturbations = nan; combine the trials with
% any perturbation value 
lfp_tfa_cfg.compare.perturbations = [0, 1]; 

% perturbation groups to be considered for pre- and post- injection
% analysis
% Should be a cell array of same size as that of lfp_tfa_cfg.compare.perturbations
% First element can be a vector of integers and corresponds to the 
% perturbation blocks to be considered for pre-injection analysis, typically 0
% Second element can be a vector of integers or 'all' or 'allbutone' and 
% corresponds to the perturbation blocks to be considered for
% post-injection analysis
% Examples:
% 1. lfp_tfa_cfg.compare.perturbation_groups = {0, 'all'}; 
% consider trials with perturbation value = 0 for pre-injection and any
% perturbation value not equal to zero for post injection
% 2. lfp_tfa_cfg.compare.perturbation_groups = {0, 'allbutone'}; 
% similar to example 1, but the trials with first post-injection 
% perturbation value is excluded (considers only from the second
% perturbation block)
% 3. lfp_tfa_cfg.compare.perturbation_groups = {0, [2, 3, 4]}; 
% consider trials with perturbation value = 0 for pre-injection and
% perturbation value = 2, 3 or 4 for post injection
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