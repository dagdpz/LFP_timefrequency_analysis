function epochs = lfp_tfa_define_epochs()
% TO be removed after testing
% lfp_tfa_define_epochs - function which defines configuration for epochs 
% of interest 
%
% USAGE: 
%   epochs = lfp_tfa_define_epochs
%
% OUTPUTS:
%   epochs - structure containing information about epochs of interest
%       name: (string) Name of the epoch of interest
%       ref_state_ID: (int) ID of the state to refer to
%       state_desc: description of the state (optional)
%       ref_period: (1x2 double) start and end time of epoch period with respect the the reference state onset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    epochs = struct();

    % Fixation hold
    epochs(1).name = 'Fhol';
    epochs(1).ref_state_ID = 6;
    epochs(1).ref_period = [-0.3 0];

    % Cue
    epochs(2).name = 'Cue';
    epochs(2).ref_state_ID = 6;
    epochs(2).ref_period = [0.05 0.2]; % time to be considered for analysis w.r.t state onset

    % Extended Delay
    epochs(3).name = 'EDel';
    epochs(3).ref_state_ID = 8;
    epochs(3).ref_period = [0.3 0.6]; % time to be considered for analysis w.r.t state onset

    % Delay
    epochs(4).name = 'Del';
    epochs(4).ref_state_ID = 4;
    epochs(4).ref_period = [-0.3 0]; % time to be considered for analysis w.r.t state onset

    % Pre Reach
    epochs(5).name = 'PreR';
    epochs(5).ref_state_ID = 62;
    epochs(5).ref_period = [-0.3 -0.05]; % time to be considered for analysis w.r.t state onset

    % Peri Reach
    epochs(6).name = 'PeriR';
    epochs(6).ref_state_ID = 63;
    epochs(6).ref_period = [-0.2 0.2]; % time to be considered for analysis w.r.t state onset
    
    % Target Hold
    epochs(7).name = 'THol';
    epochs(7).ref_state_ID = 20;
    epochs(7).ref_period = [-0.3 0]; % time to be considered for analysis w.r.t state onset

    %save(fullfile(results_fldr, 'all_states.mat'), 'all_states');
    
end
    
