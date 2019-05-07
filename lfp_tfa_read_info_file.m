function [ usable_sites_table ] = lfp_tfa_read_info_file( cfg )
%lfp_tfa_define_settings - Function to define LFP time frequency analysis settings 
%
% USAGE:
%	usable_sites_table  = lfp_tfa_read_info_file( cfg );
%
% INPUTS:
%       cfg                     - settings for reading the info file
%           Required Fields: see lfp_tfa_settings_v1
%               1. info_filepath    - absolute path to the info excel file
%               2. use_datasets     - which datasets to use
%
% OUTPUTS:
%		usable_sites_table      - table containing all entries from the
%		info excel file for which entry 'Set' is equal to one of the
%		specified values in cfg.use_datasets and entry 'Usable' is 1. 
%
% REQUIRES:	
%
% See also settings/lfp_tfa_settings_v1
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

    info_table = readtable(cfg.info_filepath);
    % filter rows based on 'Set'
    filt_dataset_table = info_table(ismember(info_table.Set, cfg.use_datasets), :);
    
    usable_sites_table = unique(filt_dataset_table(:,{'Monkey', 'Session', 'Date', ...
        'Site_ID', 'Block', 'Perturbation', 'Set', 'Usable', 'task', 'Target', ...
        'Hemisphere', 'Comments'}), 'rows');
    
    % filter rows based on 'Usable'
    usable_sites_table = usable_sites_table(usable_sites_table.Usable == 1, :);

end

