%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading processed LFP data, rejection of noise trials
% and task specific analysis using TFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% INITIALIZATION

clear;

close all; 

version = 1;

lfp_tfa_cfg = lfp_tfa_define_settings(version);

%% Function calls

% Read and process LFP data
lfp_datafiles = dir(fullfile(lfp_tfa_cfg.data_folder, 'sites_*.mat'));

%lfp_tfa_cfg.datafile_path = {};
for file = lfp_datafiles
    fprintf('Processing LFP for session %s', file.name(7:end-4));
    data_filepath = fullfile(lfp_tfa_cfg.data_folder, file.name);    
    lfp_tfa_process_LFP(data_filepath, lfp_tfa_cfg);
    %lfp_tfa_cfg.data_filepath = [lfp_tfa_cfg.data_filepath; data_filepath];
end

% Detect noisy trials
%[noisy_trials] = lfp_tfa_reject_noisy_lfp(sites_lfp_folder,
%lfp_tfa_cfg.noise);

% Compute baseline power
[sites_lfp_folder, baseline] = lfp_tfa_compute_baseline_power(sites_lfp_folder, lfp_tfa_cfg);

% Average TFR
[cond_based_tfs] = lfp_tfa_plot_average_tfr(sites_lfp_folder, lfp_tfa_cfg.analyse_states, lfp_tfa_cfg );

% Evoked LFP analysis
[ cond_based_evoked ] = lfp_tfa_plot_average_evoked_LFP( sites_lfp_folder, lfp_tfa_cfg.analyse_states, lfp_tfa_cfg ) ;

% LFP power spectrum
[ cond_based_psd ] = lfp_tfa_plot_average_powspctrum( sites_lfp_folder, lfp_tfa_cfg);


