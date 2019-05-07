function all_states = lfp_tfa_define_states()
% TO be removed after testing
% lfp_tfa_define_states - function which defines configuration for states 
%
% USAGE: 
%   all_states = lfp_tfa_define_states(results_fldr)
%
% INPUTS:
%   results_fldr - folder to save the .mat file containing state info
% OUTPUTS:
%   all_states - structure containing information about states
%       state_ID: state identifier (integer) - required
%       state_name: name of the state - required
%       state_desc: description of the state (optional)
%       tminus_onset: time in seconds before the state onset where LFP 
%       analysis should start
%       tplus_onset: time in seconds up to which the LFP analysis should be 
%       done from the state onset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    all_states = struct();

    % Fixation acquisition
    all_states(1).state_ID = 2;
    all_states(1).state_name = 'fxa';
    all_states(1).state_desc = 'Fixation acquisition';
    all_states(1).tminus_onset = 0.4; % time to be considered for analysis before state onset
    all_states(1).tplus_onset = 0.6; % time to be considered for analysis after state onset

    % Fixation hold
    all_states(2).state_ID = 3;
    all_states(2).state_name = 'fxh';
    all_states(2).state_desc = 'Fixation hold';
    all_states(2).tminus_onset = 0.4; % time to be considered for analysis before state onset
    all_states(2).tplus_onset = 0.6; % time to be considered for analysis after state onset

    % Cue onset
    all_states(3).state_ID = 6;
    all_states(3).state_name = 'cue';
    all_states(3).state_desc = 'Cue onset';
    all_states(3).tminus_onset = 1.0; % time to be considered for analysis before state onset
    all_states(3).tplus_onset = 0.5; % time to be considered for analysis after state onset

    % Reach start
    all_states(4).state_ID = 62;
    all_states(4).state_name = 'reach';
    all_states(4).state_desc = 'Reach start';
    all_states(4).tminus_onset = 0.5; % time to be considered for analysis before state onset
    all_states(4).tplus_onset = 0.5; % time to be considered for analysis after state onset

    % Target hold
    all_states(5).state_ID = 5;
    all_states(5).state_name = 'trh';
    all_states(5).state_desc = 'Target hold';
    all_states(5).tminus_onset = 0.4; % time to be considered for analysis before state onset
    all_states(5).tplus_onset = 0.6; % time to be considered for analysis after state onset

    % Reward
    all_states(5).state_ID = 21;
    all_states(5).state_name = 'rwd';
    all_states(5).state_desc = 'Reward';
    all_states(5).tminus_onset = 0.4; % time to be considered for analysis before state onset
    all_states(5).tplus_onset = 0.6; % time to be considered for analysis after state onset
    
    % Success
    all_states(6).state_ID = 20;
    all_states(6).state_name = 'scs';
    all_states(6).state_desc = 'Success';
    all_states(6).tminus_onset = 0.4; % time to be considered for analysis before state onset
    all_states(6).tplus_onset = 0.6; % time to be considered for analysis after state onset
    
    % Reward
    all_states(7).state_ID = 4;
    all_states(7).state_name = 'traq';
    all_states(7).state_desc = 'Target Acquisition';
    all_states(7).tminus_onset = 0.4; % time to be considered for analysis before state onset
    all_states(7).tplus_onset = 0.6; % time to be considered for analysis after state onset
    
    % Reward
    all_states(8).state_ID = 8;
    all_states(8).state_name = 'del';
    all_states(8).state_desc = 'Delay Period';
    all_states(8).tminus_onset = 0.4; % time to be considered for analysis before state onset
    all_states(8).tplus_onset = 0.6; % time to be considered for analysis after state onset
    
    % Reward
    all_states(9).state_ID = 63;
    all_states(9).state_name = 'reach_end';
    all_states(9).state_desc = 'Reach End';
    all_states(9).tminus_onset = 0.4; % time to be considered for analysis before state onset
    all_states(9).tplus_onset = 0.6; % time to be considered for analysis after state onset

    %save(fullfile(results_fldr, 'all_states.mat'), 'all_states');
    
end
    
