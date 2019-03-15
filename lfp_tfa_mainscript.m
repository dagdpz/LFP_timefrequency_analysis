%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading processed LFP data, rejection of noise trials
% and task specific analysis using TFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION

clear;

close all;

% select the session folder
[session_filename, pathname, ~] = uigetfile('*.mat', 'Select the mat file containing processed LFP data for the session to analyse', ...
    'MultiSelect', 'off');

% folder to save figures
root_fig_folder = [pathname '\Figures'];


% folder to save results
root_results_folder = [pathname '\Results'];

% first read in the information about states
lfp_tfa_define_states;
load('..\all_states.mat');

% load LFP data for the selected session
load(fullfile(pathname, session_filename));

trialinfo = struct();
trialinfo.start_state = 'fxa';
trialinfo.end_state = 'trh';

% maximum no:of sites to analyse
maxsites = 1; % inf = analyse all sites

%% Read the required fields from  the processed LFP data for the session
states_lfp = lfp_tfa_read_LFP(fullfile(pathname, session_filename), all_states, maxsites, root_results_folder);

%% Reject noisy trials
% configuration for lfp noise rejection
% methods to be used
cfg_noise.methods = {'amp', 'std', 'diff', 'pow'};
% threshold for lfp raw amplitude
cfg_noise.amp_thr = 5;
% number of consecutive samples beyond threshold to be considered
cfg_noise.amp_N = 5;
% no of standard deviations of trial w.r.t complete LFP std
cfg_noise.std_thr = 4;
% threshold for lfp derivative in percentile
cfg_noise.diff_thr = 4;
% number of consecutive samples beyond threshold to be considered
cfg_noise.diff_N = 5;
% threshold for lfp power in standard deviations
cfg_noise.pow_thr = 2;
% folder to save results
cfg_noise.results_folder = root_results_folder;
%cfg_noise.results_folder = [pathname '\Figures'];
% whether single trials should be plotted
cfg_noise.plottrials = 0;

[states_lfp, noisy_trials] = lfp_tfa_reject_noisy_trials(states_lfp, cfg_noise);
%filt_session_lfp = rejectNoisyLFPTrials( session_lfp )

%% Compute baseline
cfg_tfs = [];
cfg_tfs.baseline_ref_state = ''; % reference state around which baseline should be considered, leave empty to consider complete trial
cfg_tfs.baseline_period = 'trial'; % period of interest for baseline calculation, trial = complete trial period
cfg_tfs.baseline_block = 1; % consider only block 1 (control) for baseline calculation
cfg_tfs.choice_trial = 0; % whether to consider choice or instructed trials
cfg_tfs.results_folder = root_results_folder;

[states_lfp, baseline] = lfp_tfa_compute_baseline(states_lfp, cfg_tfs);

%% Prepare FT datatype fpr TFR analysis
%[ft_data_sites, session_lfp] = prepareFTdatatype(sites(5:9), analyse_states, all_states, maxsites, choice, inactivation, blocks, baseline);

%% Compute the TFR per site and average across sites
cfg_condition = [];
% define the peristates to analyse
analyse_states = {6, 62};
% analyse choice or instructed trials
cfg_condition.choice = 0; % 0 = instructed, 1 = choice, nan = both
% analyse control or inactivation trials
cfg_condition.perturbation = 3; % 1 = inactivation, 0 = control, nan = both
% blocks to be analysed, 
cfg_condition.blocks = 3;
% recorded hemispace
cfg_condition.recorded_hemispace = 'L';

% baseline configuration
cfg_baseline = [];
cfg_baseline.method = 'zscore';

[cond_based_tfs] = lfp_tfa_compute_plot_tfr(states_lfp, analyse_states, cfg_condition, cfg_baseline, root_results_folder );

%% Evoked LFP analysis
%evoked_lfp = evokedLFPAnalysis(ft_data_sites, root_fig_folder);

%% LFP power spectrum
%lfp_spectrum = plotLFPPowerSpectrum(ft_data_sites, pathname);


