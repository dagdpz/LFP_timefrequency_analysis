function lfp_tfa_cfg = lfp_tfa_define_settings(settings_filepath, maxsites)
%lfp_tfa_define_settings - Function to define LFP TFA settings 
%   Detailed explanation goes here

    lfp_tfa_global_define_states;    

    % load the specified settings file
    run(settings_filepath);
    
    % read info excel file
    lfp_tfa_cfg.sites_info = lfp_tfa_read_info_file(lfp_tfa_cfg);
     
%     lfp_tfa_cfg.data_folder = data_folder;
    % create a folder to save results
    lfp_tfa_cfg.root_results_fldr = fullfile(lfp_tfa_cfg.results_folder, ...
        date, ['ver_' num2str(lfp_tfa_cfg.version)]);
    if ~exist(lfp_tfa_cfg.root_results_fldr, 'dir')
        mkdir(lfp_tfa_cfg.root_results_fldr);
    end
    
    % load states
    %lfp_tfa_cfg.all_states = lfp_tfa_define_states();
    % first read in the information about states
    all_states = lfp_tfa_define_states();
    lfp_tfa_cfg.all_states = all_states;
    
    % load epochs
    lfp_tfa_cfg.epochs = lfp_tfa_define_epochs();
    
    % maxsites
    if nargin > 2
        lfp_tfa_cfg.maxsites = maxsites;
    end
    
    % folder to save noise rejection results
    lfp_tfa_cfg.noise.results_folder = lfp_tfa_cfg.root_results_fldr;
    % folder to save baseline results
    lfp_tfa_cfg.results_folder = lfp_tfa_cfg.root_results_fldr;

    % save settings struct
    save(fullfile(lfp_tfa_cfg.root_results_fldr, ['lfp_tfa_settings_ver_' num2str(lfp_tfa_cfg.version) '.mat']), ...
        'lfp_tfa_cfg');

end

