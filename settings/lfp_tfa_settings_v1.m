%% Initialization

% initialize configuration structure
lfp_tfa_cfg = [];
   
%% Settings for data folders
% versioning, a unique version for the settings file and analysis results
% the results produced using this settings file would be saved under 
% the folder [lfp_tfa_cfg.results_folder, '\', date, '\ver_' lfp_tfa_cfg.version]
% eg: 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data\LFP_TFA_Results\20190506\ver_SN_0.2'
lfp_tfa_cfg.version = 'Magnus_dPul_control';%'Linus_inactivation';

% sorted neurons excel file, from which information about sessions and
% individual sites can be obtained
lfp_tfa_cfg.info_filepath = 'Y:\Projects\PPC_pulv_body_signals\ephys\MIP_inactivation_20190314\Mag_sorted_neurons.xls';

% dataset to be used for analysis, see entry 'Set' in the sorted neurons excel file
% only those sessions belonging to 'Set' = lfp_tfa_cfg.use_datasets will be
% used for analysis
lfp_tfa_cfg.use_datasets = [31];

% absolute path to the folder where the results of analysis should be stored
lfp_tfa_cfg.results_folder = 'Y:\Personal\Sarath\Results\LFP_TFA_Results';

% info about sessions to be analysed
% should be a 1 x N struct, N = number of sessions to analyse
% the struct should contain the following fields:
%       Monkey:         name of monkey (string)
%       Date:           recording date (string of format YYYYMMDD)
%       Input_file:     Absolute path to the file containing LFP data for the session
%       Preinj_blocks:  Blocks to be considered for pre-injection,
%       typically 0
%       Postinj_blocks: Blocks to be considered for post-injection, can be
%       integer, array of integers, 'all' or 'allbutone' (if 'all' is
%       specified, all post-injection blocks will be combined; if
%       'allbutfirst', all blocks from the second post-injection block will
%       be combined)

lfp_tfa_cfg.session_info(1) = ...
    struct('Monkey',        'Magnus', ...
           'Date',          '20190131', ...
           'Input',         {{'Y:\Projects\PPC_pulv_body_signals\ephys\dPul_control_20190131_rest\sites_Magnus_20190131.mat', ...
                            'Y:\Projects\PPC_pulv_body_signals\ephys\dPul_control_20190131\sites_Magnus_20190131.mat'}}, ...
           'Input_ECG',     'Y:\Projects\PhysiologicalRecording\Data\Magnus\20190131\20190131_ecg.mat', ...
           'Preinj_blocks',  [0 2], ...
           'Postinj_blocks', []);
% lfp_tfa_cfg.session_info(1) = ...
%     struct('Monkey',        'Lin', ...
%            'Date',          '20170622', ...
%            'Input',         'Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn\sites_Linus_20170622.mat', ...
%            'Preinj_blocks',  0, ...
%            'Postinj_blocks', 'allbutfirst');
% lfp_tfa_cfg.session_info(2) = ...
%     struct('Monkey',        'Lin', ...
%            'Date',          '20170629', ...
%            'Input',         'Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn\sites_Linus_20170629.mat', ...
%            'Preinj_blocks',  0, ...
%            'Postinj_blocks', 'allbutfirst');
% lfp_tfa_cfg.session_info(3) = ...
%     struct('Monkey',        'Lin', ...
%            'Date',          '20170707', ...
%            'Input',         'Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn\sites_Linus_20170707.mat', ...
%            'Preinj_blocks',  0, ...
%            'Postinj_blocks', 'allbutfirst');
% lfp_tfa_cfg.session_info(4) = ...
%     struct('Monkey',        'Lin', ...
%            'Date',          '20170713', ...
%            'Input',         'Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn\sites_Linus_20170713.mat', ...
%            'Preinj_blocks',  0, ...
%            'Postinj_blocks', 'allbutfirst');
% lfp_tfa_cfg.session_info(5) = ...
%     struct('Monkey',        'Lin', ...
%            'Date',          '20170720', ...
%            'Input',         'Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn\sites_Linus_20170720.mat', ...
%            'Preinj_blocks',  0, ...
%            'Postinj_blocks', 'allbutfirst');
% lfp_tfa_cfg.session_info(6) = ...
%     struct('Monkey',        'Lin', ...
%            'Date',          '20170802', ...
%            'Input',         'Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn\sites_Linus_20170802.mat', ...
%            'Preinj_blocks',  0, ...
%            'Postinj_blocks', 'allbutfirst');
%        
% lfp_tfa_cfg.session_info(7) = ...
%     struct('Monkey',        'Lin', ...
%            'Date',          '20170804', ...
%            'Input',         'Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn\sites_Linus_20170804.mat', ...
%            'Preinj_blocks',  0, ...
%            'Postinj_blocks', 'allbutfirst');
%        
% lfp_tfa_cfg.session_info(8) = ...
%     struct('Monkey',        'Lin', ...
%            'Date',          '20170818', ...
%            'Input',         'Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn\sites_Linus_20170818.mat', ...
%            'Preinj_blocks',  0, ...
%            'Postinj_blocks', 'allbutfirst');
%        
% lfp_tfa_cfg.session_info(9) = ...
%     struct('Monkey',        'Lin', ...
%            'Date',          '20171012', ...
%            'Input',         'Y:\Projects\PPC_pulv_eye_hand\ephys\MIP_dPul_inj_working_post_sfn\sites_Linus_20171012.mat', ...
%            'Preinj_blocks',  0, ...
%            'Postinj_blocks', 'allbutfirst');
       
       
% To add a new session to analyse, increment the counter by 1 and add a new
% value into the lfp_tfa_cfg.session_info struct
% Example: 
% lfp_tfa_cfg.session_info(3) = ...
%     struct('Monkey',        'Magnus', ...
%            'Date',          '20190208', ...
%            'Input',         'Y:\Projects\PPC_pulv_body_signals\ephys\MIP_inactivation_20190208\sites_Magnus_20190208.mat', ...
%            'Preinj_blocks',  0, ...
%            'Postinj_blocks', 'allbutfirst');

% what kind of analyses should be done on LFP
% should be a cell array of strings which indicate which kind of analyses
% should be performed on LFP
% Currently supported analyses are 'tfs', 'evoked', 'pow', and 'sync'
%       'tfs'       - LFP time frequency spectrogram average for given conditions and time windows
%       'evoked'    - LFP evoked response average for given conditions and time windows
%       'pow'       - LFP power spectrum average for given conditions and epochs
%       'sync'      - LFP-LFP phase synchronization measure for given conditions and
%           time windows
lfp_tfa_cfg.analyses = {'tfs', 'evoked', 'pow'}; % , 'sync', 'syncspctrm'

% targets to be included in the analysis
% should be a cell array of strings which indicate the target names
% the target names should be same as the target field in the LFP data
% structure
% Those targets which are not in the analysed sessions will be ignored
% Example:
% 1. lfp_tfa_cfg.compare.targets = {'MIPa_R', 'MIPa_L', 'dPul_R', 'dPul_L'}; 
lfp_tfa_cfg.compare.targets = {'dPul_R'};%{'MIP_R', 'MIP_L'}; 

% target pairs to be included for LFP-LFP sychronization
% should be a 1xN cell array of 1x2 cell array of strings which indicate
% the target pairs between which the LFP-LFP phase synchronization should
% be calculated - valid only if LFP-LFP phase sync should be calculated
if any(strcmp(lfp_tfa_cfg.analyses, 'sync') | strcmp(lfp_tfa_cfg.analyses, 'syncspctrm'))
    lfp_tfa_cfg.compare.target_pairs = {{'MIP_R', 'MIP_R'}, {'MIP_R', 'MIP_L'}, ...
        {'MIP_L', 'MIP_L'}}; 
end

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

% random seed for random number generator for reproducibility
lfp_tfa_cfg.random_seed = rng;

%% Settings for averaging TFR and evoked LFP based on conditions

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
% 4. lfp_tfa_cfg.compare.reach_hands = {'any'}; ignore hand label (trial with
% any hand label is combined)
lfp_tfa_cfg.compare.reach_hands = {'L', 'R'};

% reach space to be included for analysis
% should be a cell array that contain only values 'R', 'L', or 'any'
% Examples:
% 1. lfp_tfa_cfg.compare.reach_spaces = {'L'}; include only those trials in
% which acquired target is on left
% 2. lfp_tfa_cfg.compare.reach_hands = {'R'}; include only those trials in
% which acquired target is on right
% 3. lfp_tfa_cfg.compare.reach_hands = {'L', 'R'}; analyse the trials in
% which acquired target is on left and on right separately
% 4. lfp_tfa_cfg.compare.reach_hands = {'any'}; ignore space label (trial with
% any acquired target position is combined)
lfp_tfa_cfg.compare.reach_spaces = {'L', 'R'}; 

% hand space combinations to be excluded from analysis
% should be a cell array with each element containing the hand and space
% label to be excluded
% if no hand-space conditions are to be excluded, leave empty
% Example:
% 1. lfp_tfa_cfg.compare.exclude_handspace = {'LR', 'RL'};
% exclude left hand right space and right hand left space trials
% 2. lfp_tfa_cfg.compare.exclude_handspace = {'LL', 'RR'};
% exclude left hand left space and right hand right space trials
% 3. lfp_tfa_cfg.compare.exclude_handspace = {'RL', 'RR'};
% exclude right hand left space and right hand right space trials
lfp_tfa_cfg.compare.exclude_handspace = {};

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

% differences in conditions to be analysed
% add new entries for further difference calculations
% each entry should be a cell array of type 1xN cell where N should be even. 
% The entries should be made as a string representing the condition 
% (eg: 'perturbation', 'choice') followed by a cell array containing 
% the two values of the given condition to be compared. 
% Examples:
% lfp_tfa_cfg.diff_condition(1) = {{'perturbation', {0, 1}}};
% This computes the difference between trials with perturbation = 1 (post-injection) and perturbation = 0 (pre-injection) when other conditions (eg. task type, choice) remains the same. 
% lfp_tfa_cfg.diff_condition(2) = {{'choice', {0, 1}}};
% This computes the difference between choice = 1 and choice = 0 (instructed) trials when other conditions (eg. task type, perturbation) remains the same. 
% lfp_tfa_cfg.diff_condition(3) = {{'type_eff', {[4 4], [4 6]}}};
% This computes the difference between given task types trials when other conditions (eg. task type, perturbation) remains the same. 
% lfp_tfa_cfg.diff_condition(4) = {{'perturbation', {0, 1}, ...
%    'choice', {0, 1}}};
% Compute difference between difference between post and pre-injection trials of choice trials and that of instructed trials     

lfp_tfa_cfg.diff_condition = {};
lfp_tfa_cfg.diff_condition(1) = {{'perturbation', {0, 1}}};
% lfp_tfa_cfg.diff_condition(2) = {{'choice', {0, 1}}};
% lfp_tfa_cfg.diff_condition(3) = {{'type_eff', {[4 4], [4 4]}}};
% lfp_tfa_cfg.diff_condition(3) = {{'perturbation', {0, 1}, ...
%     'choice', {0, 1}}};

%% Time information

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
lfp_tfa_cfg.trialinfo.end_state = lfp_tfa_states.SUCCESS;

% offset to be considered from the onset of
% trial end reference state for calculating the trial end time
% i.e., trial end time = onset of trial end state + end offset
% Example:
% 1. lfp_tfa_cfg.trialinfo.ref_tend = 0.5;
% trial start time = onset time of lfp_tfa_cfg.trialinfo.end_state + 0.5;
% 1. lfp_tfa_cfg.trialinfo.ref_tend = -0.5;
% trial start time = onset time of lfp_tfa_cfg.trialinfo.end_state - 0.5;
lfp_tfa_cfg.trialinfo.ref_tend = 0;


%% Settings for ft_freqanalysis in FieldTrip
% Configuration for calculating LFP time frequency spectrogram using
% ft_freqanalysis function of the fieldtrip toolbox

% method for calculating the LFP power spectra
% can be 'mtmfft', 'mtmconvol', 'wavelet'
% Example:
% 1. lfp_tfa_cfg.tfr.method = 'wavelet'; % implements wavelet time frequency
%                        transformation (using Morlet wavelets) based on
%                        multiplication in the frequency domain.
% 2. lfp_tfa_cfg.tfr.method = 'mtmfft', analyses an entire spectrum for the entire data
%                        length, implements multitaper frequency transformation.
% see http://www.fieldtriptoolbox.org/reference/ft_freqanalysis/ for more
% details
lfp_tfa_cfg.tfr.method          = 'wavelet';  

% frequencies of interest (in Hz)
% Example: 
% 1. lfp_tfa_cfg.tfr.foi = logspace(log10(2), log10(120), 60); 60 logspaced
% frequencies from 2Hz to 120 Hz
lfp_tfa_cfg.tfr.foi             = logspace(log10(2), log10(120), 60);

% number of lfp samples to step for the sliding time window
% Example:
% lfp_tfa_cfg.tfr.timestep  = 25; 
% the sliding time window steps by an amount equal to 25 lfp samples. 
lfp_tfa_cfg.tfr.timestep        = 25; 

% depending on the method chosen, other configurations vary

% For method = 'wavelet', Ignored for other methods
% width of the wavelets in number of cycles
% Making the value smaller will increase the temporal resolution at the expense of frequency resolution and vice versa
% Wavelet duration = width / F / pi (F = frequency), wavelet duration
% decreases with frequency
% for width = 6, frequency F = 30 Hz, wavelet duration = 6/30/pi = 0.0637 s
% Example: 
% 1. lfp_tfa_cfg.tfr.width = 6; % wavelet of width 6 cycles 
lfp_tfa_cfg.tfr.width           = 6; 

% For method = 'mtmfft' or 'mtmconvol', Ignored for method = 'wavelet'

% taper (single or multiple) to be used 
% can be 'dpss', 'hanning' or many others
% 'hanning' - conventional single taper
% 'dpss' - multiple tapers based on discrete prolate spheroidal sequences 
% (DPSS), also known as the Slepian sequence
lfp_tfa_cfg.tfr.taper           = [];

% the width of frequency smoothing in Hz (fw)
% Note that 4 Hz smoothing means plus-minus 4 Hz, i.e. a 8 Hz smoothing box.
% should be a vector of size 1 x numfoi, Leave empty for method = 'wavelet'
% Example: 
% 1. lfp_tfa_cfg.tfr.tapsmofrq  = 0.4 *cfg.foi; 
%the smoothing will increase with frequency.
lfp_tfa_cfg.tfr.tapsmofrq       = [];

% length of the sliding time-window in seconds (= tw)
% should be vector of length 1 x numfoi
% Following relation must hold: 
% K = 2twfw-1, where K is required to be larger than 0
% K is the number of tapers applied
% Example:
% lfp_tfa_cfg.tfr.t_ftimwin  = 5./cfg.foi; % 5 cycles per window
% window length decreases with frequency
lfp_tfa_cfg.tfr.t_ftimwin       = [];

%% settings for LFP-LFP sync analysis
% measure of LFP-LFP phase synchronization
% can be only 'ppc' currently
% 'ppc' calculates pairwise phase consistency
% entry will be used as cfg.method for performing ft_connectivityanalysis
lfp_tfa_cfg.sync.measure = 'ppc';

%% Settings to detect noisy trials
% configuration for lfp noise rejection
lfp_tfa_cfg.noise = [];
% whether or not to apply noise rejection 
% Set to 0 to accept all trials
% Set to 1 to run the noise trial detection methods
lfp_tfa_cfg.noise.detect = 1;
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

if isempty(lfp_tfa_cfg.baseline_ref_state)
	lfp_tfa_cfg.baseline_ref_period = 'trial';
else
	lfp_tfa_cfg.baseline_ref_period = []; % SET LIMITS OF baseline_ref_period here
end

% which perturbation blocks to be considered for baseline power calculation
% set to 0 for considering only pre-injection blocks
% in case trials from a single perturbation block is analysed, the same
% block will be used for baseline calculation
% Examples: 
% 1. lfp_tfa_cfg.baseline_perturbation = 0; only perturbation block 0
% (pre-injection) is used for baseline power calculation
% 2. lfp_tfa_cfg.baseline_perturbation = [0, 2]; combines perturbation blocks
% 0 (pre-injection) and 2 (post-injection), but is not recommended
% 3. lfp_tfa_cfg.baseline_perturbation = [2, 3]; combines perturbation blocks
% 2 and 3
if length(lfp_tfa_cfg.compare.perturbations) == 1
    lfp_tfa_cfg.baseline_perturbation = lfp_tfa_cfg.compare.perturbations;
else
    lfp_tfa_cfg.baseline_perturbation = 0; % set the perturbation block(s) to be used for computing baseline
end

% whether to consider choice (1) or instructed trials (0) in baseline power
% calculation 
% in case, only either choice/instructed is being analysed, the trials from
% the analysed choice is used for baseline
% Examples:
% 1. lfp_tfa_cfg.baseline_use_choice_trial = 0; % consider only instructed trials
% 2. lfp_tfa_cfg.baseline_use_choice_trial = 1; % consider only choice trials
% 3. lfp_tfa_cfg.baseline_use_choice_trial = [0, 1]; % consider both instructed and choice trials
if length(lfp_tfa_cfg.compare.choice_trials) == 1
    lfp_tfa_cfg.baseline_use_choice_trial = lfp_tfa_cfg.compare.choice_trials;
else
    % set choice(1) and/or instructed(0) to be used for computing baseline
    lfp_tfa_cfg.baseline_use_choice_trial = 0; 
end

% which task type to be used in baseline power calculation 
% in case, only one task type is being analysed, the trials from
% the analysed task type is used for baseline
% Examples:
% 1. lfp_tfa_cfg.baseline_use_type = 4; % consider only trials of type=4
% 2. lfp_tfa_cfg.baseline_use_type = [4, 1]; % consider trials of
% both type=4 and type=1(not implemented - for future use)
% 3. lfp_tfa_cfg.baseline_use_type = inf; % consider trials of
% any type 
if length(lfp_tfa_cfg.compare.types) == 1
    lfp_tfa_cfg.baseline_use_type = lfp_tfa_cfg.compare.types;
else
    % set task type(s) to be used for computing baseline
    lfp_tfa_cfg.baseline_use_type = 1; 
end

% which effector to be used in baseline power calculation 
% in case, only one efffector is being analysed, the trials from
% the analysed effector is used for baseline
% Examples:
% 1. lfp_tfa_cfg.baseline_use_effector = 0; % consider only trials with
% effector = 0
% 2. lfp_tfa_cfg.baseline_use_effector = inf; % consider trials with any
% effector value
% 3. lfp_tfa_cfg.baseline_use_effector = [0, 6]; % consider only trials with
% effector = 0 and effector = 6 (not implemented - for future use)
if length(lfp_tfa_cfg.compare.effectors) == 1
    lfp_tfa_cfg.baseline_use_effector = lfp_tfa_cfg.compare.effectors;
else
    % set effector(s) to be used for computing baseline
    lfp_tfa_cfg.baseline_use_effector = 0; 
end

% define the time windows to analyse for LFP TFR and evoked LFP response
% Must be a Nx4 cell array, N = number of windows to analyse
% Each row corresponds to one state and contain following elements
% 1. Identifier of state around which the window is referenced, 
% see lfp_tfa_global_states, Example:  lfp_tfa_states.CUE_ON
% 2. Name of the reference state (window) - string (used for labeling 
% purposes in plots) eg: 'Cue'
% 3. Start time offset - offset in seconds from reference state onset for 
% the start of time window
% start time = Reference state onset time + Start time offset
% 4. End time offset - offset in seconds from ref. state onset for end of
% time window
% end time = Ref. state onset time + end time offset
% 
% Example row: 
%   lfp_tfa_states.CUE_ON,     'Cue',    -1.0 ,    0.5
lfp_tfa_cfg.analyse_states = {lfp_tfa_states.CUE_ON,    'Cue',      -0.5,   0.9;...
                             lfp_tfa_states.REA_INI,    'Reach',    -0.3,   0.5};                    

% define the epochs to analyse for LFP power spectrum
% Must be a Nx4 cell array, N = number of epochs to analyse
% Each row corresponds to one epoch and contain following elements
% 1. Identifier of state to which the epoch is referred, see lfp_tfa_global_states, Example:  lfp_tfa_states.CUE_ON
% 2. Name of the epoch - string (used for labeling purposes in plots) eg: 'FHol'
% 3. Start time offset - offset in seconds from reference state onset for 
% the epoch start
% Epoch start time = Reference state onset time + Start time offset
% 4. End time offset - offset in seconds from ref. state onset for epoch
% end
% Epoch end time = Ref. state onset time + end time offset
% Example row: 
%   lfp_tfa_states.CUE_ON,     'FHol',    -0.3 ,    0
lfp_tfa_cfg.analyse_epochs = {lfp_tfa_states.CUE_ON,     'FHol',    -0.3 ,    0  ;...
                              lfp_tfa_states.CUE_ON,     'Cue' ,    0.05 ,    0.2 ; ...
                              lfp_tfa_states.DEL_PER,    'EDel',    0.3 ,     0.6 ; ...
                              lfp_tfa_states.TAR_ACQ,    'Del',     -0.3 ,    0  ; ...
                              lfp_tfa_states.REA_INI,    'PreR',    -0.3 ,    -0.05 ; ...
                              lfp_tfa_states.REA_END,    'PeriR',   -0.2 ,    0.2 ; ...
                              lfp_tfa_states.SUCCESS,    'THol',    -0.3 ,    0    };

% minimum number of trials per condition to be satisfied to consider a site
% for averaging, if for a site, for any condition, the  number of valid 
% (non-noisy) trials is less than mintrials_percondition, the site is not considered for averaging
% Set lfp_tfa_cfg.mintrials_percondition = 1 to consider a site if atleast
% one valid trial is obtained (keep minimum value of 1)
% Example:
% consider those sites with atleast 5 trials for each condition
% lfp_tfa_cfg.mintrials_percondition = 5; 
% By condition, we mean a combination of choice/instr, pre/post-injection, type and effector, hand-space
 lfp_tfa_cfg.mintrials_percondition = 0;

% method to be used for baseline normalization
% can be 'zscore', 'relchange', 'subtraction', 'division'
% 'zscore' - computes Z-score for each time freq bin
% Z(t,f) = (P(t, f) - mu_P(f)) / (sigma_P(f))
% 'relchange' - relative power change w.r.t. the baseline power
% P_norm(t,f) = (P(t, f) - mu_P(f)) / (mu_P(f))
% 'subtraction' - absolute increase in power w.r.t. the baseline
% P_norm(t,f) = (P(t, f) - mu_P(f))
% 'division' - relative increase in power w.r.t. the baseline
% P_norm(t,f) = (P(t, f)) / (mu_P(f))
% Example:
% lfp_tfa_cfg.baseline_method = 'relchange';
lfp_tfa_cfg.baseline_method = 'zscore';

% flag to indicate if LFP TFR average should be computed - for future use
% Set to 0 if LFP TFR average should not be computed, else set to 1
% lfp_tfa_cfg.compute_tfr = 1;

% flag to indicate if LFP evoked response average should be computed - for future use
% Set to 0 if LFP evoked response average should not be computed, else set to 1
% lfp_tfa_cfg.compute_evoked = 1;

% flag to indicate if LFP power spectrum average should be computed - for future use
% Set to 0 if LFP power spectrum average should not be computed, else set to 1
% lfp_tfa_cfg.compute_pow = 1;

    
%% Settings for averaging across sessions or sites

% how to average data across multiple sessions/sites
% 'sessions' - average the session averages (a session average is the
% average of site averages within a session)
% 'sites' - average across sites, regardless of which session they come from
% Example: lfp_tfa_cfg.compute_avg_across = 'sites'
% Example: lfp_tfa_cfg.compute_avg_across = {'sessions', 'sites'};  compute
% both averages across session averages and across site averages
lfp_tfa_cfg.compute_avg_across = {'sessions', 'sites'}; 