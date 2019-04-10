function [ cond_based_psd ] = lfp_tfa_plot_average_powspctrum( sites_lfp_folder, lfp_tfa_cfg ) 

% lfp_tfa_plot_average_powspctrum  - plots lfp power spectrum for
% different hand-space tuning conditions for each site and across sites
%
% USAGE:
%	[ cond_based_psd ] = lfp_tfa_plot_average_powspctrum( sites_lfp_folder, lfp_tfa_cfg )
%
% INPUTS:
%		sites_lfp_folder  	- folder containing lfp data
%       lfp_tfa_cfg     - struct containing configuration for TFR 
%           Required fields:
%           trial_condition.blocks              - blocks to be analysed, 
%           leave empty to analyse each block separately
%           trial_condition.baseline_method     - method to be used for 
%           baseline normalization ('zscore', 'subtraction', 'division',
%           'relchange')
%           results_folder                      - folder to save results
%
% OUTPUTS:
%		cond_based_psd	- output structure which saves the average power spectrum for  
%                         trials of a given condition for different handspace 
%                         tunings and periods around the states analysed
%                         same datastructure as input ft_data_sites, but
%                         with additional field to store condition-wise lfp
%                         time freq response		
%
%
%
% See also lfp_tfa_compute_baseline, lfp_tfa_define_states
    
    % suppress warning for xticklabel
    warning ('off', 'MATLAB:hg:willberemoved');

    % make a folder to save figures
    %sessionName = states_lfp(1).session;
    results_folder_tfr = fullfile(lfp_tfa_cfg.results_folder, 'Condition_based_LFP_Power');
    if ~exist(results_folder_tfr, 'dir')
        mkdir(results_folder_tfr);
    end
    
    states_lfp = struct();
    % load site lfps
    site_lfp_files = dir([sites_lfp_folder, '\*.mat']);
    for i = 1:length( site_lfp_files )
        load(fullfile(sites_lfp_folder, site_lfp_files(i).name));
        if i == 1
            states_lfp = site_lfp;
        else
            states_lfp = [states_lfp, site_lfp];
        end
    end
    
    recorded_hemispace = unique([states_lfp.recorded_hemispace]);
    choice = unique([states_lfp(1).trials.choice_trial]);
    perturbation = unique([states_lfp(1).trials.perturbation]);
    blocks = unique([states_lfp(1).trials.block]);
    
    % create conditions
    cfg_conditions = struct();
       
    i = 0;
    for rec_hem = recorded_hemispace        
        for c = choice
            for b = blocks
                i = i + 1;
                cfg_conditions(i).recorded_hemispace = rec_hem;
                cfg_conditions(i).choice = c;
                cfg_conditions(i).block = b;
                cfg_conditions(i).perturbation = ...
                    unique([states_lfp(1).trials([states_lfp(1).trials.block] == b).perturbation]);
                cond_label = [];
                if cfg_conditions(i).recorded_hemispace == 'L'
                    cond_label = [cond_label 'Left_hemisphere_'];
                else
                    cond_label = [cond_label 'Right_hemisphere_'];
                end
                if cfg_conditions(i).choice == 0
                    cond_label = [cond_label 'Instructed_'];
                else
                    cond_label = [cond_label 'Choice_'];
                end
                if cfg_conditions(i).perturbation == 0
                    cond_label = [cond_label 'Control_'];
                else
                    cond_label = [cond_label 'Inactivation_'];
                end
                cond_label = [cond_label, 'Block_', num2str(cfg_conditions(i).block)];
                cfg_conditions(i).label = cond_label;
                
                % create a folder for storing results for this condition
                cfg_conditions(i).results_folder = fullfile(results_folder_tfr); %, cfg_conditions(i).label
                if ~exist(cfg_conditions(i).results_folder, 'dir')
                    mkdir(cfg_conditions(i).results_folder)
                end
                                
            end
        end
    end
    
    
    if isfield(lfp_tfa_cfg, 'add_conditions')
        for c = 1:length(lfp_tfa_cfg.add_conditions)
            if ~isempty(lfp_tfa_cfg.add_conditions(c))
                if ~isempty(lfp_tfa_cfg.add_conditions(c).blocks)
                    for rec_hem = recorded_hemispace        
                        for ch = choice
                            i = i + 1;
                            if strcmp(lfp_tfa_cfg.add_conditions(c).blocks, 'inactivation')
                                cfg_conditions(i).block = blocks(perturbation ~= 0);
                                cfg_conditions(i).perturbation = 1;
                            else
                                cfg_conditions(i).block = lfp_tfa_cfg.add_conditions(c).blocks;
                                cfg_conditions(i).perturbation = perturbation(blocks == lfp_tfa_cfg.add_conditions(c).blocks(1));
                                
                            end                    
                            cfg_conditions(i).choice = ch;
                            if isfield(lfp_tfa_cfg.add_conditions(c), 'perturbation')
                                cfg_conditions(i).perturbation = lfp_tfa_cfg.add_conditions(c).perturbation;
                            end
                            cfg_conditions(i).recorded_hemispace = rec_hem;
                            cond_label = [];
                            if cfg_conditions(i).recorded_hemispace == 'L'
                                cond_label = [cond_label 'Left_hemispace_'];
                            else
                                cond_label = [cond_label 'Right_hemispace_'];
                            end
                            if cfg_conditions(i).choice == 0
                                cond_label = [cond_label 'Instructed_'];
                            else
                                cond_label = [cond_label 'Choice_'];
                            end
                            if cfg_conditions(i).perturbation == 0
                                cond_label = [cond_label 'Control_'];
                            else
                                cond_label = [cond_label 'Inactivation_'];
                            end
                            cond_label = [cond_label, 'Block_', num2str(cfg_conditions(i).block)];
                            cfg_conditions(i).label = cond_label;
                            
                            % create a folder for storing results for this condition
                            cfg_conditions(i).results_folder = fullfile(results_folder_tfr); %, cfg_conditions(i).label
                            if ~exist(cfg_conditions(i).results_folder, 'dir')
                                mkdir(cfg_conditions(i).results_folder)
                            end
                        end
                    end
                end
            end
        end
    end
    
       
    % condition based TFS
    cond_based_psd = struct();
    for cn = 1:length(cfg_conditions)
        
        cond_based_psd(cn).label = cfg_conditions(cn).label;
        cond_based_psd(cn).cfg_cond = cfg_conditions(cn);
        
        % states to be analysed
    %     analyse_states = {states_lfp(1).states.name};

        % hand-space tuning of LFP
        hs_labels = unique({states_lfp(1).trials.hndspc_lbl});
        
        % current trial condition analysed
        cfg_condition = cfg_conditions(cn);
        
        
        % cell array to store time frequency average across sites
        %TFR_avg = cell(length(analyse_states),length(hs_labels));

        % loop through each site
        nsites = length(states_lfp);
        
        cond_based_psd(cn).sites = struct();

        % loop through each site
        for i = 1:length(states_lfp)            
            
            %cond_based_tfs(cn).sites(i) = struct();
            % consider site based on recorded hemispace
            if strcmp(states_lfp(i).recorded_hemispace, cfg_conditions(cn).recorded_hemispace) 
                cond_based_psd(cn).sites(i).site_ID = states_lfp(i).site_ID;
                cond_based_psd(cn).sites(i).session = states_lfp(i).session;
                cond_based_psd(cn).sites(i).target = states_lfp(i).target;
                % make a struct for concatenating TFR for all states
                cond_based_psd(cn).sites(i).tfs_avg_site = struct(); 
                % struct to store evoked LFP
                cond_based_psd(cn).sites(i).site_evoked_lfp = struct();
                % struct to store LFP power spectrum
                cond_based_psd(cn).sites(i).site_lfp_psd = struct();
                %cond_based_tfs(i).tfs_avg_site.powspctrm = cell(length(analyse_states),length(hs_labels));
                cond_based_psd(cn).sites(i).ntrials = zeros(1,length(hs_labels));

                for hs = 1:length(hs_labels)
                    cond_trials = zeros(1, length(states_lfp(i).trials));
                    % get the trials for given condition and this hs label
%                     if ~isnan(cfg_conditions(cn).perturbation)
%                         cond_trials = cond_trials & ...
%                             ([states_lfp(i).trials.perturbation]) == ...
%                             (cfg_conditions(cn).perturbation);
%                     end
                    if ~isnan(cfg_conditions(cn).block)
                        for b = cfg_conditions(cn).block
                            cond_trials = cond_trials | ...
                                ([states_lfp(i).trials.block] == b);
                        end
                    end
                    if ~isnan(cfg_conditions(cn).choice)
                        cond_trials = cond_trials & ...
                            ([states_lfp(i).trials.choice_trial] == ...
                            cfg_conditions(cn).choice);
                    end
                    cond_trials = cond_trials & ...
                        strcmp({states_lfp(i).trials.hndspc_lbl}, hs_labels(hs));
                    cond_based_psd(cn).sites(i).ntrials(hs) = sum(cond_trials);
                                        
                    fprintf('Condition %s - %s\n', cfg_conditions(cn).label, hs_labels{hs});
                    fprintf('Total number of trials %g\n', sum(cond_trials));
                    
                    cond_based_psd(cn).sites(i).noisytrials(hs) = ...
                        sum(cond_trials & [states_lfp(i).trials.noisy]); 
                                        
                    % consider only non noisy trials
                    fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                        & [states_lfp(i).trials.noisy]));
                    cond_trials = cond_trials & ~[states_lfp(i).trials.noisy];
                    
                                 
                    % loop through trials 

                    for ep = 1:length(lfp_tfa_cfg.epochs)
                        epoch_tfs.psd = {}; % power spectrum
                        epoch_tfs.psd_f = {}; % power spectrum freq

                        for t = find(cond_trials)

                            epochs          = states_lfp(i).trials(t).epochs;
                            %states          = states_lfp(i).trials(t).states;
                            epoch_onset_t   = epochs(ep).onset_t;
                            epoch_start_t   = epochs(ep).start_t;
                            epoch_end_t     = epochs(ep).end_t;
                            % sampling frequency
                            %fs = states_lfp(i).trials(t).fsample;

                            % LFP power spectrum
                            epoch_powspctrm = states_lfp(i).trials(t).tfs.powspctrm(1, :, ...
                                (states_lfp(i).trials(t).tfs.time >= epoch_start_t & ...
                                states_lfp(i).trials(t).tfs.time <= epoch_end_t));  
                            epoch_tfs.psd = [epoch_tfs.psd sum(epoch_powspctrm, 3)];
                            epoch_tfs.psd_f = states_lfp(i).trials(t).tfs.freq;

                        end
                                               
                        if ~isempty(epoch_tfs.psd)
                            % power spectrum average
                            arr_state_psd = vertcat(epoch_tfs.psd{:});
                            epoch_psd_mean = nanmean(arr_state_psd, 1);
                            
                            % save LFP power spectrum
                            cond_based_psd(cn).sites(i).site_lfp_psd(ep, hs).mean = epoch_psd_mean;
                            cond_based_psd(cn).sites(i).site_lfp_psd(ep, hs).freq = epoch_tfs.psd_f;
                            cond_based_psd(cn).sites(i).site_lfp_psd(ep, hs).trials = find(cond_trials);
                            cond_based_psd(cn).sites(i).site_lfp_psd(ep, hs).hs_label = hs_labels(hs);
                            cond_based_psd(cn).sites(i).site_lfp_psd(ep, hs).epoch_name = lfp_tfa_cfg.epochs(ep).name;
                            
                        end

                    end

                end
            else
                continue;

            end

        end
        
    
    
        
        
        % calculate the average LFP power spectrum across sites for this condition 
        if isfield(cond_based_psd(cn).sites, 'site_lfp_psd')
            % struct to store average evoked LFP across sites
            cond_based_psd(cn).session_lfp_psd = struct();%cell(length(analyse_states),length(hs_labels));
            
            
            for i = 1:length(cond_based_psd(cn).sites)

                for hs = 1:length(hs_labels)
                    for ep = 1:length(lfp_tfa_cfg.epochs)
                        if ~isempty(cond_based_psd(cn).sites(i).site_lfp_psd(ep, hs))
                            
                            if i == 1
                                cond_based_psd(cn).session_lfp_psd(ep, hs).mean = ...
                                    (1/length(cond_based_psd(cn).sites))*cond_based_psd(cn).sites(i).site_lfp_psd(ep, hs).mean ;
                                cond_based_psd(cn).session_lfp_psd(ep, hs).freq = ...
                                    cond_based_psd(cn).sites(i).site_lfp_psd(ep, hs).freq;
                            else
                                if ~isempty(cond_based_psd(cn).session_lfp_psd(ep, hs).mean)
                                    nfreqs = length(cond_based_psd(cn).session_lfp_psd(ep, hs).freq);
                                    if nfreqs > length(cond_based_psd(cn).sites(i).site_lfp_psd(ep, hs).freq)
                                        nfreqs = length(cond_based_psd(cn).sites(i).site_lfp_psd(ep, hs).freq);
                                    end                               
                                    cond_based_psd(cn).session_lfp_psd(ep, hs).mean = ...
                                        cond_based_psd(cn).session_lfp_psd(ep, hs).mean(1:nfreqs) + ...
                                        (1/length(cond_based_psd(cn).sites)) * cond_based_psd(cn).sites(i).site_lfp_psd(ep, hs).mean(1:nfreqs) ;
                                else
                                    cond_based_psd(cn).session_lfp_psd(ep, hs).mean = ...
                                        (1/length(cond_based_psd(cn).sites))*cond_based_psd(cn).sites(i).site_lfp_psd(ep, hs).mean ;
                                end
                            end
                            cond_based_psd(cn).session_lfp_psd(ep, hs).hs_label = ...
                                cond_based_psd(cn).sites(i).site_lfp_psd(ep, hs).hs_label;
                            cond_based_psd(cn).session_lfp_psd(ep, hs).epoch_name = ...
                                cond_based_psd(cn).sites(i).site_lfp_psd(ep, hs).epoch_name;
                            cond_based_psd(cn).session_lfp_psd(ep, hs).nsites = ...
                                length(cond_based_psd(cn).sites);
                        end
                        
                    end
                end
            end
        end
    
        % plots

        % site averages
        for i = 1:length(cond_based_psd(cn).sites) 
            % TFR
            if isfield(cond_based_psd(cn).sites, 'tfs_avg_site')

                plottitle = ['Site ID: ', cond_based_psd(cn).sites(i).site_ID ...
                    ', Target = ' cond_based_psd(cn).sites(i).target ', '  ...
                '(block ' num2str(cfg_conditions(cn).block) '), '];
                if cfg_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                else
                    plottitle = [plottitle 'Choice trials'];
                end
                result_folder = fullfile(cfg_conditions(cn).results_folder, cond_based_psd(cn).sites(i).site_ID);
                if ~exist(result_folder, 'dir')
                    mkdir(result_folder);
                end
                result_file = fullfile(cfg_conditions(cn).results_folder, cond_based_psd(cn).sites(i).site_ID, ...
                    ['Avg_Pow_Spctrum_' cond_based_psd(cn).sites(i).site_ID '_' cfg_conditions(cn).label '.png']);
                lfp_tfa_plot_hs_tuned_psd(cond_based_psd(cn).sites(i).site_lfp_psd, ...
                    lfp_tfa_cfg, plottitle, result_file);
            end
        end

        % plot average power spectrum across sites for this session
        plottitle = ['Session: ', cond_based_psd(cn).sites(1).session ...
            ', Target = ' cond_based_psd(cn).sites(1).target ', '  ...
            'Block ' num2str(cfg_conditions(cn).block) ', '];
        if cfg_conditions(cn).choice == 0
            plottitle = [plottitle 'Instructed trials'];
        else
            plottitle = [plottitle 'Choice trials'];
        end
        results_file = fullfile(cfg_conditions(cn).results_folder, ...
            ['Avg_Pow_Spctrum_' cond_based_psd(cn).sites(1).session '_' cfg_conditions(cn).label '.png']);
        lfp_tfa_plot_hs_tuned_psd(cond_based_psd(cn).session_lfp_psd, ...
                    lfp_tfa_cfg, plottitle, results_file);
        
    end
    
    close all;
    % save mat files
    save(fullfile(results_folder_tfr, 'cfg_conditions.mat'), 'cfg_conditions');
    save(fullfile(results_folder_tfr, 'cond_based_psd.mat'), 'cond_based_psd');
    save(fullfile(results_folder_tfr, 'lfp_tfa_cfg.mat'), 'lfp_tfa_cfg');
end
        
