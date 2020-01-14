%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading LFP data, rejection of noise trials
% and condition specific analysis using TFR, evoked, spectra, sync
% spectrograms and sync spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; 

% folder containing settings file 
cfg_folderpath = 'C:\Users\snair\Documents\GitHub\LFP_timefrequency_analysis\settings';
% files containing settings for LFP analysis
settings_filepaths = {'lfp_decoding_settings_example.m'};

settings_files = {};
if strcmp(settings_files, '*.m')
    settings_files = dir(fullfile(cfg_folderpath, '*.m'));
    for f = 1:length(settings_files)
        file = settings_files(f).name;
        settings_files = [settings_files, file];
    end
else
    settings_files = settings_filepaths;
end

% loop through each file in the cfg folder
for f = 1:length(settings_files)
    settings_filepath = fullfile(cfg_folderpath, settings_files{f});
    fprintf('%s\n', settings_filepath);
    % INITIALIZATION
    close all;

    % read settings file
    lfp_decode_cfg = lfp_tfa_decode_define_settings(settings_filepath);

    % Read processed LFP for decoding
    if lfp_decode_cfg.read_decode_LFP
        lfp_decode = lfp_tfa_decode_get_conditions_lfp( lfp_decode_cfg );
    else
        fprintf('Loading LFP data to decode ...\n'); 
        load(lfp_decode_cfg.decode_lfp_file, 'lfp_decode');
    end

    % Decode LFP
    fprintf('Decoding LFP ...\n');
    lfp_decode_accuracy = lfp_tfa_decode_predict_classes( lfp_decode, lfp_decode_cfg, lfp_decode_cfg.analyses);
    fprintf('done.\n');
end

