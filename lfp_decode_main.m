%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script for LFP time frequency analysis
% Runs functions for reading LFP data, rejection of noise trials
% and condition specific analysis using TFR, evoked, spectra, sync
% spectrograms and sync spectra
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; 

% folder containing settings file 
cfg_folderpath = 'C:\Users\snair\Documents\GitHub\LFP_timefrequency_analysis\settings\LFP_Decoding_Linus_8sessions';
% file containing settings for LFP analysis
% should have the same format as settings/lfp_tfa_settings_example.m
% settings_files = {'lfp_decoding_settings_MIPR_Instr_control_IH_ISvsCS.m', ...
%     'lfp_decoding_settings_MIPR_Instr_control_CH_ISvsCS.m', ...
%     'lfp_decoding_settings_MIPR_Instr_control_IHvsCH_IS.m', ...
%     'lfp_decoding_settings_MIPR_Instr_control_IHvsCH_CS.m', ...
%     'lfp_decoding_settings_MIPR_Instr_control_LHvsRH.m', ...
%     'lfp_decoding_settings_MIPR_Instr_control_LSvsRS.m', ...
%     'lfp_decoding_settings_MIPR_Instr_inactivation_IH_ISvsCS.m', ...
%     'lfp_decoding_settings_MIPR_Instr_inactivation_CH_ISvsCS.m', ...
%     'lfp_decoding_settings_MIPR_Instr_inactivation_IHvsCH_IS.m', ...
%     'lfp_decoding_settings_MIPR_Instr_inactivation_IHvsCH_CS.m', ...
%     'lfp_decoding_settings_MIPR_Instr_inactivation_LHvsRH.m', ...
%     'lfp_decoding_settings_MIPR_Instr_inactivation_LSvsRS.m'};

% loop through each file in the cfg folder
settings_files = dir(fullfile(cfg_folderpath, '*.m'));
for f = 1:length(settings_files)
    settings_filepath = fullfile(cfg_folderpath, settings_files(f).name);
    fprintf('%s\n', settings_filepath);
    % INITIALIZATION
    close all;

    % read settings file
    lfp_tfa_cfg = lfp_tfa_decode_define_settings(settings_filepath);

    % Read processed LFP for decoding
    if lfp_tfa_cfg.read_decode_LFP
        lfp_decode = lfp_tfa_decode_get_conditions_lfp( lfp_tfa_cfg );
    else
        fprintf('Loading LFP data to decode ...\n'); 
        load(lfp_tfa_cfg.decode_lfp_file, 'lfp_decode');
    end

    % Decode LFP
    fprintf('Decoding LFP ...\n');
    lfp_decode = lfp_tfa_decode_predict_classes( lfp_decode, lfp_tfa_cfg, {'lfp_tfs'});
    fprintf('done.\n');
end

