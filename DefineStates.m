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

% Cue onset
all_states(1).state_ID = 6;
all_states(1).state_name = 'cue_on';
all_states(1).state_desc = 'Cue onset';
all_states(1).tminus_onset = 0.4; % time to be considered for analysis before state onset
all_states(1).tplus_onset = 0.6; % time to be considered for analysis after state onset

% Reach start
all_states(2).state_ID = 62;
all_states(2).state_name = 'reach_start';
all_states(2).state_desc = 'Reach start';
all_states(2).tminus_onset = 0.4; % time to be considered for analysis before state onset
all_states(2).tplus_onset = 0.6; % time to be considered for analysis after state onset

save all_states;