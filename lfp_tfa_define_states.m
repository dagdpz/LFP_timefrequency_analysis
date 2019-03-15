%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DefineStates.m 
% Script which defines information about states in a struct
% Define all states to be analysed
% state_ID: state identifier (integer) - required
% state_name: name of the state - required
% state_desc: description of the state
% tminus_onset: time in seconds before the state onset where LFP analysis
% should start (for future use)
% tplus_onset: time in seconds up to which the LFP analysis should be done
% from the state onset (for future use)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear;

all_states = struct();

% Trial start
% all_states(1).state_ID = 1;
% all_states(1).state_name = 'trial_ini';
% all_states(1).state_desc = 'Initialization of a trial';
% all_states(1).tminus_onset = 0; % time to be considered for analysis before state onset
% all_states(1).tplus_onset = 0.5; % time to be considered for analysis after state onset
% 
% % Fixation acq
% all_states(2).state_ID = 2;
% all_states(2).state_name = 'fix_acq';
% all_states(2).state_desc = 'Fixation acquisition';

% % Fixation hold
% all_states(3).state_ID = 3;
% all_states(3).state_name = 'fix_hold';
% all_states(3).state_desc = 'Fixation hold';
% 

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
all_states(3).tminus_onset = 0.6; % time to be considered for analysis before state onset
all_states(3).tplus_onset = 0.6; % time to be considered for analysis after state onset

% Reach start
all_states(4).state_ID = 62;
all_states(4).state_name = 'reach';
all_states(4).state_desc = 'Reach start';
all_states(4).tminus_onset = 0.3; % time to be considered for analysis before state onset
all_states(4).tplus_onset = 0.6; % time to be considered for analysis after state onset

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

save all_states;
    
