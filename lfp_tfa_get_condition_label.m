function condition_label = lfp_tfa_get_condition_label(cfg_condition, label_length)
%lfp_tfa_get_condition_label - Get a label for a given trial condition (A
%trial condition is a permutation of perturbation, choice, type and effector)
%
% USAGE:
%   condition_label = lfp_tfa_get_condition_label(cfg_condition)
%	condition_label = lfp_tfa_get_condition_label(cfg_condition, label_length)
%
% INPUTS:
%       cfg_condition   - structure containing the conditions
%           Required fields:
%           type            - trial type
%           effector        - trial effector
%           choice          - choice/instructed trial
%           perturbation    - pre or post injection trial
%       label_length    - length of the condition label ('long' or
%       'short'), 'short' for a shorter label
%
% OUTPUTS:
%		condition_label  - label for the given trial condition
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

switch cfg_condition.type
    case 1
        typ='Fixation';
        typs='F';
    case 2
        typ='Visually_guided';
        typs='V';
    case 3
        typ='Memory';
        typs='M';
    case 4
        typ='Delay';
        typs='D';
    case 5
        typ='M2S';
        typs='M2';
    case 6
        typ='S2S';
        typs='S2';
end

switch cfg_condition.effector
    case 0
        eff='saccades';
        effs='sac';
    case 1
        eff='free_gaze_reaches';
        effs='fgr';
    case 2
        eff='joint';
        effs='joi';
    case 3
        eff='diss_saccades';
        effs='dsa';
    case 4
        eff='diss_reaches';
        effs='dre';
    case 5
        eff='sen_reaches';
        effs='sen';
    case 6
        eff='central_fix_reaches';
        effs='cfr';
end

if isnan(cfg_condition.choice)
    ch = 'All choice';
    chs = 'Allch';
elseif cfg_condition.choice == 1
    ch='Choice';
    chs='Choi';
elseif cfg_condition.choice == 0
    ch='Instructed';
    chs='Instr';
elseif strcmp(cfg_condition.choice, 'diff1')
    ch='(Instructed - Choice)';
    chs = '(Instr - Choi)';      
end 

if isnan( cfg_condition.perturbation)
    pert = 'All perturbation';
    perts = 'Allperturb';
elseif cfg_condition.perturbation == 0
    pert='Pre-injection';
    perts='Pre';
elseif cfg_condition.perturbation ==  1
    pert='Post-injection';
    perts='Post';
elseif strcmp(cfg_condition.perturbation,  'diff1')
    pert='(Post-injection - Pre-injection)';
    perts = '(Post - Pre)';    
        
end

condition_label=[typ, ' ', eff, ' ', ch, ' ', pert];
condition_label_short=[typs, effs, ' ', chs, ' ', perts];

if strcmp(label_length, 'short')
    condition_label = condition_label_short;
end    

end