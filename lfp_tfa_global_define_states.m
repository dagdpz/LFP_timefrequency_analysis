% Script which defines the global states

global lfp_tfa_states;

lfp_tfa_states.INI_TRI       = 1; % initialize trial
lfp_tfa_states.FIX_ACQ       = 2; % fixation acquisition
lfp_tfa_states.FIX_HOL       = 3; % fixation hold
lfp_tfa_states.TAR_ACQ       = 4; % target acquisition
lfp_tfa_states.TAR_HOL       = 5; % target hold
lfp_tfa_states.CUE_ON        = 6; % cue on
lfp_tfa_states.MEM_PER       = 7; % memory period
lfp_tfa_states.DEL_PER       = 8; % delay period
lfp_tfa_states.TAR_ACQ_INV   = 9; % target acquisition invisible
lfp_tfa_states.TAR_HOL_INV   = 10; % target hold invisible
lfp_tfa_states.MAT_ACQ       = 11; % target acquisition in sample to match
lfp_tfa_states.MAT_HOL       = 12; % target acquisition in sample to match
lfp_tfa_states.MAT_ACQ_MSK   = 13; % target acquisition in sample to match
lfp_tfa_states.MAT_HOL_MSK   = 14; % target acquisition in sample to match
lfp_tfa_states.SEN_RET       = 15; % return to sensors for poffenberger
lfp_tfa_states.ABORT         = 19;
lfp_tfa_states.SUCCESS       = 20;
lfp_tfa_states.REWARD        = 21;
lfp_tfa_states.ITI           = 50;
lfp_tfa_states.SAC_INI       = 60; % saccade start
lfp_tfa_states.SAC_END       = 61;
lfp_tfa_states.REA_INI       = 62; %reach start
lfp_tfa_states.REA_END       = 63;
lfp_tfa_states.TRI_END       = 90;
lfp_tfa_states.ITI_END       = 98;
lfp_tfa_states.CLOSE         = 99;



