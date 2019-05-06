function [ usable_sites_table ] = lfp_tfa_read_info_file( cfg )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    info_table = readtable(cfg.info_filepath);
    % filter datasets
    filt_dataset_table = info_table(ismember(info_table.Set, cfg.use_datasets), :);
    
    usable_sites_table = unique(filt_dataset_table(:,{'Monkey', 'Session', 'Date', ...
        'Site_ID', 'Block', 'Perturbation', 'Set', 'Usable', 'task', 'Target', ...
        'Hemisphere', 'Comments'}), 'rows');
    
    usable_sites_table = usable_sites_table(usable_sites_table.Usable == 1, :);

end

