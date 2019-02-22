function [ ft_data_sites, session_lfp ] = prepareFTdatatype( sites, analyse_states, ...
    all_states, maxsites, choice, inactivation, baseline )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PrepareLFPDataPerSite.m
    % Script to read in the sitewise LFP data and store the required
    % information into a new struct sites_lfp. 
    % Struct members:
    % session: name of session
    % site_ID: site identifier
    % ipsilateral: ipsilateral hand for this site
    % ipsilesional: ipsilesional hand for this site (for future use)
    % trial: (1xM) struct to store information about trials recorded from the site
    % 
    %       block: block/run number
    %       choice_trial: 0 (instructed), 2 (choice trial with 2 targets)
    %       reach_hand: hand used to reach the target
    %       reach_space: space in which target which is reached appeared
    %       lfp_data: lfp data recorded for this trial (1xN)
    %       time: timestamps of lfp data (1xN)
    %       fsample: sampling freq of LFP data
    %       tsample: sampling time of LFP data
    %       hndspc_lbl: hand-space tuning for the trial
    %       states: (1xK) struct which stores the information about states to
    %       be analysed
    %               name: name of the state (as defined in DefineStates.m_
    %               onset_t: state onset time
    %               onset_s: sample number of state onset
    %               start_s: sample at which LFP analysis should start
    %               end_s: sample at which LFP analysis should end
    %               
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    close all; 
    
    % prepare a root data structure for all sites
    session_lfp = struct();

    % Read in LFP data for the session
    fprintf('Reading processed LFP data');
    %load(session_filename);

    % save data inside struct 
    % first loop through each site
    for i = 1:min(length(sites), maxsites)
        session_lfp(i).session = sites(i).site_ID(1:end-8);
        session_lfp(i).site_ID = sites(i).site_ID;
        session_lfp(i).target = sites(i).target;
        session_lfp(i).ipsilateral = session_lfp(i).target(end);
        % create a ft_datatype_raw struct for storing LFP data for all trials
        % see ft_rawdatatype
        % sites_lfp(i).lfp_data = struct();
        % now loop through each trial for this site
        comp_trial = 0; % iterator for completed trials
        for t = 1:length(sites(i).trial)
            success = sites(i).trial(t).success;
            noisy = sites(i).trial(t).noisy;
            if success && ~noisy
                comp_trial = comp_trial + 1;
                block = sites(i).trial(t).block;
                choice_trial = sites(i).trial(t).choice;
                reach_hand = sites(i).trial(t).reach_hand; % 1 = left, 2 = right
                perturbation = sites(i).trial(t).perturbation; % 0 = control
                
                % falg indicating whether this trial should be considered
                % or not based on whether instructed or choice, control or
                % inactivation
                consider_trial = true;
                % select choice or instructed trials based on user
                % selection
                if ~isnan(choice)
                    if choice % user wants to select choice trials
                        consider_trial = consider_trial & choice_trial ~= 0;
                    else
                        consider_trial = consider_trial & choice_trial == 0;
                    end
                end
                % select choice or instructed trials based on user
                % selection
%                 if inactivation % user wants to select inactivat trials
%                     consider_trial = consider_trial & perturbation ~= 0;
%                 else
%                     consider_trial = consider_trial & perturbation == 0;
%                 end
                if ~consider_trial
                    comp_trial = comp_trial - 1;
                    continue;
                end
                
                % check if the trials have all the states to be analysed
                trial_states = sites(i).trial(t).states;
                
                if ~all(ismember([analyse_states{:}], trial_states)) 
                    comp_trial = comp_trial - 1;
                    continue
                end                    
                
                if reach_hand == 1, reach_hand = 'L'; else reach_hand = 'R'; end
                tar_pos = sites(i).trial(t).tar_pos;
                fix_pos = sites(i).trial(t).fix_pos;
                if sign(real(tar_pos) - real(fix_pos)) == -1
                    reach_space = 'L'; 
                else
                    reach_space = 'R';
                end
                % assign hand-space
                if session_lfp(i).ipsilateral == reach_hand  
                    if session_lfp(i).ipsilateral == reach_space
                        hs_label = 'IH IS';
                    else
                        hs_label = 'IH CS';
                    end
                else 
                    if session_lfp(i).ipsilateral == reach_space
                        hs_label = 'CH IS';
                    else
                        hs_label = 'CH CS';
                    end
                end

                %states = trial(i).states;
                start_time = (sites(i).trial(t).TDT_LFPx_tStart); % trial start time
                fs = sites(i).trial(t).TDT_LFPx_SR; % sample rate
                LFP = sites(i).trial(t).LFP; % LFP data
                ts = (1/fs); % sample time
                nsamples = numel(LFP);
                end_time = start_time + (ts*(nsamples-1));
                timestamps = linspace(start_time, end_time, nsamples);

                % save retrieved data into struct
                session_lfp(i).trials(comp_trial).block = block;
                session_lfp(i).trials(comp_trial).choice_trial = choice_trial;
                session_lfp(i).trials(comp_trial).lfp_data = LFP;
                session_lfp(i).trials(comp_trial).time = timestamps;
                session_lfp(i).trials(comp_trial).fsample  = fs;
                session_lfp(i).trials(comp_trial).tsample = ts;
                session_lfp(i).trials(comp_trial).reach_hand  = reach_hand;
                session_lfp(i).trials(comp_trial).reach_space  = reach_space;
                session_lfp(i).trials(comp_trial).hndspc_lbl  = hs_label;
                session_lfp(i).trials(comp_trial).perturbation  = perturbation;


    %             % get movement onset time
    %             reach_start = sites(i).trial(t).states_onset(sites(i).trial(t).states == 62);
    %             % give cue_onset zero timestamp 
    %             reach_start_sample = find(abs(timestamps - reach_start) == min(abs(timestamps - reach_start)));

                
                % get state onset times and onset samples
                session_lfp(i).trials(comp_trial).states = struct();
                for s = 1:length(analyse_states)
                    % get state onset time
                    state_onset = sites(i).trial(t).states_onset(sites(i).trial(t).states == all_states(s).state_ID);
                    % get sample number of state onset time
                    state_onset_sample = find(abs(timestamps - state_onset) == min(abs(timestamps - state_onset)));
                    % save into struct
                    session_lfp(i).trials(comp_trial).states(s).name  = all_states(s).state_name;
                    session_lfp(i).trials(comp_trial).states(s).onset_t  = state_onset;
                    session_lfp(i).trials(comp_trial).states(s).onset_s  = state_onset_sample;
                    session_lfp(i).trials(comp_trial).states(s).start_s  = max(1, state_onset_sample - 400);
                    session_lfp(i).trials(comp_trial).states(s).end_s = min(length(timestamps), state_onset_sample + 600);
                end
                
                % define baseline
                session_lfp(i).trials(comp_trial).baseline = struct();
                baseline_ref_s = session_lfp(i).trials(comp_trial).states(...
                    strcmp({session_lfp(i).trials(comp_trial).states.name}, ...
                    baseline.ref_state)).onset_s;
                session_lfp(i).trials(comp_trial).baseline.ref_s = ...
                    baseline_ref_s;
                session_lfp(i).trials(comp_trial).baseline.start_s = ...
                    baseline_ref_s + sample(baseline.period(1), 0, ts);
                session_lfp(i).trials(comp_trial).baseline.end_s = ...
                    baseline_ref_s + sample(baseline.period(2), 0, ts);

                
            end
        end
    end
        
    %sites_lfp(i).lfp_data.fsample = fs;
    ft_data_sites = struct();
    for i = 1:length(session_lfp)
        site = session_lfp(i);
        ft_data_sites(i).session = site.session;
        ft_data_sites(i).target = site.target;
        ft_data_sites(i).siteID = site.site_ID;
        ft_data_sites(i).ipsilateral = site.ipsilateral;
        ft_data_sites(i).states = struct();
        hs_labels = unique({session_lfp(i).trials.hndspc_lbl});
        ft_data_sites(i).TFR_hs = cell(length(analyse_states),length(hs_labels));
        for st = 1:length(analyse_states) % add 1 for baseline
            % store name of the state
            ft_data_sites(i).states(st).name = all_states([all_states.state_ID] == analyse_states{st}).state_name;
            % struct to save hand-space tuned FT type LFP data for this state
            ft_data_sites(i).states(st).hndspc = struct();
            % initialize hand-space tuning structures
            for hs = 1:length( hs_labels)
                ft_data_sites(i).states(st).hndspc(hs).label = hs_labels(hs);
                ft_data_sites(i).states(st).hndspc(hs).trial = {};
                ft_data_sites(i).states(st).hndspc(hs).time  = {};
            end
        end
        % initialize struct to store baseline
        ft_data_sites(i).baseline = struct();
        % struct to save hand-space tuned FT type LFP data for this state
        ft_data_sites(i).baseline.hndspc = struct();
        % initialize hand-space tuning structures
        for hs = 1:length( hs_labels)
            ft_data_sites(i).baseline.hndspc(hs).label = hs_labels(hs);
            ft_data_sites(i).baseline.hndspc(hs).trial = {};
            ft_data_sites(i).baseline.hndspc(hs).time  = {};
            ft_data_sites(i).baseline.hndspc(hs).sampleinfo = [];
        end
        % loop through each trial for this site
        blocks = [];
        ntrials = 0;
        for t = 1:length(site.trials)
            trial = site.trials(t);
            
            ntrials = ntrials + 1;
            blocks = [blocks, trial.block];
            hs_tuning = trial.hndspc_lbl;
            hs_idx = find(strcmp(hs_labels, hs_tuning));

            % loop through states
            for st = 1:length(trial.states)
                state = trial.states(st);
                st_idx = find(strcmp({ft_data_sites(i).states.name}, state.name));
                %ft_data_sites(i).states(st_idx).hndspc(hs_idx).label = hs_tuning;
                ft_data_sites(i).states(st_idx).hndspc(hs_idx).trial = [ft_data_sites(i).states(st_idx).hndspc(hs_idx).trial; ...
                    trial.lfp_data(state.start_s:state.end_s)];
                ft_data_sites(i).states(st_idx).hndspc(hs_idx).time = [ft_data_sites(i).states(st_idx).hndspc(hs_idx).time; ...
                    trial.time(state.start_s:state.end_s) - trial.time(state.onset_s)];
                ft_data_sites(i).states(st_idx).hndspc(hs_idx).fsample = trial.fsample;
            end   
            
            % baseline
            baseline = trial.baseline;
            
            %st_idx = length(ft_data_sites(i).states)+1;
            %ft_data_sites(i).states(st_idx).name = 'baseline';
            ft_data_sites(i).baseline.hndspc(hs_idx).label = {hs_tuning};
            ft_data_sites(i).baseline.hndspc(hs_idx).trial = [ft_data_sites(i).baseline.hndspc(hs_idx).trial; ...
                trial.lfp_data(baseline.start_s:baseline.end_s)];
            ft_data_sites(i).baseline.hndspc(hs_idx).time = [ft_data_sites(i).baseline.hndspc(hs_idx).time; ...
                trial.time(baseline.start_s:baseline.end_s) - trial.time(baseline.ref_s)];
            ft_data_sites(i).baseline.hndspc(hs_idx).fsample = trial.fsample;
            
        end
        ft_data_sites(i).nblocks = length(unique(blocks));
        ft_data_sites(i).ntrials = ntrials;
        
    end

    % save data
    save session_lfp session_lfp;
    save ft_data_sites ft_data_sites;

end

