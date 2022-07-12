function lfp_tfa_cfg = lfp_tfa_define_settings(settings_filepath)
%lfp_tfa_define_settings - Function to define LFP time frequency analysis settings 
%
% USAGE:
%	lfp_tfa_cfg = lfp_tfa_define_settings(settings_filepath)
%
% INPUTS:
%       settings_filepath         - absolute path to the matlab script file
%       where LFP TFA settings are defined, see
%       settings/lfp_tfa_settings_example
%
% OUTPUTS:
%		lfp_tfa_cfg               - structure containing all settings
%
% REQUIRES:	lfp_tfa_read_info_file, lfp_tfa_compare_conditions,
%
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_read_info_file, 
% lfp_tfa_compare_conditions, lfp_tfa_mainscript
%
% Author(s):	S.Nair, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2019-02-15:	Created function (Sarath Nair)
% 2019-03-05:	First Revision
% ...
% $Revision: 1.0 $  $Date: 2019-03-05 17:18:00 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%

    % add external functions to path
    % hopefully this is not needed any more!!
    temp_path=which('lfp_tfa_define_settings');
    end_idx=strfind(temp_path,filesep);
    temp_path=temp_path(1:end_idx(end));
    addpath(genpath(fullfile(temp_path, 'external')));

    % define state IDs
    lfp_tfa_global_define_states;    

    % load the specified settings file
    run(settings_filepath);
    
    % set random seed for reproducibility
    rng(lfp_tfa_cfg.random_seed);
    
    % read info excel file (Sorted neurons file)
    % lfp_tfa_cfg.sites_info = lfp_tfa_read_info_file(lfp_tfa_cfg);
     
    % create a root folder to save results of the analysis
    % root_results_folder = [lfp_tfa_cfg.results_folder, '\', date, '\ver_' lfp_tfa_cfg.version]
    % eg: 'C:\Data\MIP_timefreq_analysis\LFP_timefrequency_analysis\Data\LFP_TFA_Results\20190506\ver_SN_0.2'
    lfp_tfa_cfg.root_results_fldr = fullfile(lfp_tfa_cfg.results_folder, ...
        lfp_tfa_cfg.version);
    if ~exist(lfp_tfa_cfg.root_results_fldr, 'dir')
        mkdir(lfp_tfa_cfg.root_results_fldr);
    end
    
    % get conditions to be included in the analysis
    lfp_tfa_cfg.conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg);
    
    % folder to save noise rejection results
    lfp_tfa_cfg.noise.results_folder = lfp_tfa_cfg.root_results_fldr;
    % folder to save baseline results
    lfp_tfa_cfg.results_folder = lfp_tfa_cfg.root_results_fldr;
    
    % folder to save LFP processing results
    if ~lfp_tfa_cfg.process_LFP && ~exist(lfp_tfa_cfg.proc_lfp_folder, 'dir')
        lfp_tfa_cfg.process_LFP = true;
    end
    if lfp_tfa_cfg.process_LFP
        lfp_tfa_cfg.proc_lfp_folder = [lfp_tfa_cfg.root_results_fldr filesep 'Processed LFP'];
    end
    % folder to save LFP analysis results
    if ~lfp_tfa_cfg.compute_site_average && ~exist(lfp_tfa_cfg.analyse_lfp_folder, 'dir')
        lfp_tfa_cfg.compute_site_average = true;
    end
    %lfp_tfa_cfg.analyse_lfp_folder = [lfp_tfa_cfg.root_results_fldr filesep 'LFP Analysis'];
    if lfp_tfa_cfg.compute_site_average
        lfp_tfa_cfg.analyse_lfp_folder = [lfp_tfa_cfg.root_results_fldr filesep 'LFP Analysis'];
    end
    % folder to store LFP-LFP sync results
    lfp_tfa_cfg.lfp_sync_folder = [lfp_tfa_cfg.analyse_lfp_folder filesep 'Condition_based_Sync'];
    % folder to store session-wise analysis results
    for i = 1:length(lfp_tfa_cfg.session_info)
        lfp_tfa_cfg.session_info(i).session = ...
            [lfp_tfa_cfg.session_info(i).Monkey, '_', lfp_tfa_cfg.session_info(i).Date];
        lfp_tfa_cfg.session_info(i).proc_results_fldr = ...
                fullfile(lfp_tfa_cfg.proc_lfp_folder, lfp_tfa_cfg.session_info(i).session);
        lfp_tfa_cfg.session_info(i).lfp_sync_results_fldr = [lfp_tfa_cfg.analyse_lfp_folder ...
            filesep lfp_tfa_cfg.session_info(i).session filesep 'Condition_based_Sync' ];
        lfp_tfa_cfg.session_info(i).lfp_syncspctrm_results_fldr = [lfp_tfa_cfg.analyse_lfp_folder ...
            filesep lfp_tfa_cfg.session_info(i).session filesep 'Condition_based_Syncspectrum' ];
        lfp_tfa_cfg.session_info(i).lfp_tfs_results_fldr = [lfp_tfa_cfg.analyse_lfp_folder ...
            filesep lfp_tfa_cfg.session_info(i).session filesep 'Condition_based_TFS' ];
        lfp_tfa_cfg.session_info(i).lfp_evoked_results_fldr = [lfp_tfa_cfg.analyse_lfp_folder ...
            filesep lfp_tfa_cfg.session_info(i).session filesep 'Condition_based_Evoked_LFP' ];
        lfp_tfa_cfg.session_info(i).lfp_pow_results_fldr = [lfp_tfa_cfg.analyse_lfp_folder ...
            filesep lfp_tfa_cfg.session_info(i).session filesep 'Condition_based_LFP_Power' ];
    end

    % save settings struct
    save(fullfile(lfp_tfa_cfg.root_results_fldr, ['lfp_tfa_settings_' lfp_tfa_cfg.version '.mat']), ...
        'lfp_tfa_cfg');

end

