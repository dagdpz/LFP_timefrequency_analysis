function [ session_pow ] = lfp_tfa_plot_site_powspctrum( states_lfp, lfp_tfa_cfg ) 

% lfp_tfa_plot_average_powspctrum  - plots lfp power spectrum for
% different hand-space tuning conditions for each site and across sites
%
% USAGE:
%	[ session_pow ] = lfp_tfa_plot_average_powspctrum( sites_lfp_folder, lfp_tfa_cfg )
%
% INPUTS:
%		states_lfp  	- structure containing processed lfp data for all
%		sites of one session, output of lfp_tfa_process_lfp or
%		lfp_tfa_reject_noisy_lfp or lfp_tfa_compute_baseline_power
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
%		session_pow     - output structure which saves the average power spectrum for  
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
    results_folder_psd = fullfile(lfp_tfa_cfg.session_results_fldr, 'Condition_based_LFP_Power');
    if ~exist(results_folder_psd, 'dir')
        mkdir(results_folder_psd);
    end
    
%     states_lfp = struct();
%     load site lfps
%     site_lfp_files = dir([sites_lfp_folder, '\*.mat']);
%     for i = 1:length( site_lfp_files )
%         load(fullfile(sites_lfp_folder, site_lfp_files(i).name));
%         if i == 1
%             states_lfp = site_lfp;
%         else
%             states_lfp = [states_lfp, site_lfp];
%         end
%     end
    
%     recorded_hemisphere = unique([states_lfp.recorded_hemisphere]);
%     choice = unique([states_lfp(1).trials.choice_trial]);
%     perturbation = unique([states_lfp(1).trials.perturbation]);
%     blocks = unique([states_lfp(1).trials.block]);
%     
%     % create conditions
%     cfg_conditions = struct();
%        
%     i = 0;
%     for rec_hem = recorded_hemisphere        
%         for c = choice
%             for b = blocks
%                 i = i + 1;
%                 cfg_conditions(i).recorded_hemisphere = rec_hem;
%                 cfg_conditions(i).choice = c;
%                 cfg_conditions(i).block = b;
%                 cfg_conditions(i).perturbation = ...
%                     unique([states_lfp(1).trials([states_lfp(1).trials.block] == b).perturbation]);
%                 cond_label = [];
%                 if cfg_conditions(i).recorded_hemisphere == 'L'
%                     cond_label = [cond_label 'Left_hemisphere_'];
%                 else
%                     cond_label = [cond_label 'Right_hemisphere_'];
%                 end
%                 if cfg_conditions(i).choice == 0
%                     cond_label = [cond_label 'Instructed_'];
%                 else
%                     cond_label = [cond_label 'Choice_'];
%                 end
%                 if cfg_conditions(i).perturbation == 0
%                     cond_label = [cond_label 'Control_'];
%                 else
%                     cond_label = [cond_label 'Inactivation_'];
%                 end
%                 cond_label = [cond_label, 'Block_', num2str(cfg_conditions(i).block)];
%                 cfg_conditions(i).label = cond_label;
%                 
%                 % create a folder for storing results for this condition
% %                 cfg_conditions(i).results_folder = fullfile(results_folder_tfr); %, cfg_conditions(i).label
% %                 if ~exist(cfg_conditions(i).results_folder, 'dir')
% %                     mkdir(cfg_conditions(i).results_folder)
% %                 end
%                                 
%             end
%         end
%     end
%     
%     
%     if isfield(lfp_tfa_cfg, 'add_conditions')
%         for c = 1:length(lfp_tfa_cfg.add_conditions)
%             if ~isempty(lfp_tfa_cfg.add_conditions(c))
%                 if ~isempty(lfp_tfa_cfg.add_conditions(c).blocks)
%                     for rec_hem = recorded_hemisphere        
%                         for ch = choice
%                             i = i + 1;
%                             if strcmp(lfp_tfa_cfg.add_conditions(c).blocks, 'inactivation')
%                                 cfg_conditions(i).block = blocks(perturbation ~= 0);
%                                 cfg_conditions(i).perturbation = 1;
%                             else
%                                 cfg_conditions(i).block = lfp_tfa_cfg.add_conditions(c).blocks;
%                                 cfg_conditions(i).perturbation = perturbation(blocks == lfp_tfa_cfg.add_conditions(c).blocks(1));
%                                 
%                             end                    
%                             cfg_conditions(i).choice = ch;
%                             if isfield(lfp_tfa_cfg.add_conditions(c), 'perturbation')
%                                 cfg_conditions(i).perturbation = lfp_tfa_cfg.add_conditions(c).perturbation;
%                             end
%                             cfg_conditions(i).recorded_hemispace = rec_hem;
%                             cond_label = [];
%                             if cfg_conditions(i).recorded_hemispace == 'L'
%                                 cond_label = [cond_label 'Left_hemispace_'];
%                             else
%                                 cond_label = [cond_label 'Right_hemispace_'];
%                             end
%                             if cfg_conditions(i).choice == 0
%                                 cond_label = [cond_label 'Instructed_'];
%                             else
%                                 cond_label = [cond_label 'Choice_'];
%                             end
%                             if cfg_conditions(i).perturbation == 0
%                                 cond_label = [cond_label 'Control_'];
%                             else
%                                 cond_label = [cond_label 'Inactivation_'];
%                             end
%                             cond_label = [cond_label, 'Block_', num2str(cfg_conditions(i).block)];
%                             cfg_conditions(i).label = cond_label;
%                             
%                             % create a folder for storing results for this condition
% %                             cfg_conditions(i).results_folder = fullfile(results_folder_tfr); %, cfg_conditions(i).label
% %                             if ~exist(cfg_conditions(i).results_folder, 'dir')
% %                                 mkdir(cfg_conditions(i).results_folder)
% %                             end
%                         end
%                     end
%                 end
%             end
%         end
%     end
    
    %cfg_conditions = lfp_tfa_compare_conditions(states_lfp, lfp_tfa_cfg);
       
    % condition based TFS
    sites_pow = struct();
    session_pow = struct();
    session_pow.session = states_lfp(1).session;
    % loop through each site
    for i = 1:length(states_lfp)  
        % target for this site
        target = unique({states_lfp(i).target});
        site_conditions = lfp_tfa_compare_conditions(lfp_tfa_cfg, target);
        
        % folder to save sitewise results
        site_results_folder = fullfile(results_folder_psd, states_lfp(i).site_ID);
        if ~exist(site_results_folder, 'dir')
            mkdir(site_results_folder);
        end
        sites_pow(i).condition = struct();
        sites_pow(i).site_ID = states_lfp(i).site_ID;
        sites_pow(i).session = states_lfp(i).session;
        sites_pow(i).target = states_lfp(i).target;
        % flag to indicate if this site should be used for
        % averaging based on minimum no:of trials per condition
        sites_pow(i).use_for_avg = 1;
        
        for cn = 1:length(site_conditions)

            sites_pow(i).condition(cn).label = site_conditions(cn).label;
            sites_pow(i).condition(cn).cfg_cond = site_conditions(cn);

            % states to be analysed
        %     analyse_states = {states_lfp(1).states.name};

            % hand-space tuning of LFP
            hs_labels = unique({states_lfp(1).trials.hndspc_lbl});

            % current trial condition analysed
            cfg_condition = site_conditions(cn);


            % cell array to store time frequency average across sites
            %TFR_avg = cell(length(analyse_states),length(hs_labels));

            % loop through each site
            nsites = length(states_lfp); 
                 
            
            % consider site based on recorded hemispace
            sites_pow(i).condition(cn).label = site_conditions(cn).label;
            sites_pow(i).condition(cn).cfg_condition = site_conditions(cn);
            % make a struct for storing LFP power spectrum for all states
            sites_pow(i).condition(cn).hs_tuned_tfs = struct(); 
            sites_pow(i).condition(cn).ntrials = zeros(1,length(hs_labels));           

            for hs = 1:length(hs_labels)
                cond_trials = lfp_tfa_get_condition_trials(states_lfp(i), site_conditions(cn));

                cond_trials = cond_trials & ...
                    strcmp({states_lfp(i).trials.hndspc_lbl}, hs_labels(hs));
                sites_pow(i).condition(cn).ntrials(hs) = sum(cond_trials);

                fprintf('Condition %s - %s\n', site_conditions(cn).label, hs_labels{hs});
                fprintf('Total number of trials %g\n', sum(cond_trials));

                sites_pow(i).condition(cn).noisytrials(hs) = ...
                    sum(cond_trials & [states_lfp(i).trials.noisy]); 

                % consider only non noisy trials
                fprintf('Number of noisy trials %g\n', sum(cond_trials ...
                    & [states_lfp(i).trials.noisy]));
                cond_trials = cond_trials & ~[states_lfp(i).trials.noisy];
                % check if the site contains a specified minimum number
                % of trials for all conditions
                if sum(cond_trials) < lfp_tfa_cfg.mintrials_percondition
                    sites_pow(i).use_for_avg = 0;
                end


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
                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).mean = epoch_psd_mean;
                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).freq = epoch_tfs.psd_f;
                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).trials = find(cond_trials);
                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).hs_label = hs_labels(hs);
                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).epoch_name = lfp_tfa_cfg.epochs(ep).name;

                    end

                end

            end
%             else
%                 continue;
% 
%             end
            
            % TFR
            if ~isempty(sites_pow(i).condition(cn).hs_tuned_power)

                plottitle = ['Site ID: ', sites_pow(i).site_ID ...
                    ', Target = ' sites_pow(i).target ', '  ...
                '(Perturb ' num2str(site_conditions(cn).perturbation_group{1}) '), '];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                else
                    plottitle = [plottitle 'Choice trials'];
                end
                site_results_folder = fullfile(results_folder_psd, 'sites', sites_pow(i).site_ID);
                if ~exist(site_results_folder, 'dir')
                    mkdir(site_results_folder);
                end
                result_file = fullfile(site_results_folder, ...
                    ['LFP_Power_' states_lfp(i).site_ID '_' site_conditions(cn).label '.png']);
                lfp_tfa_plot_hs_tuned_psd(sites_pow(i).condition(cn).hs_tuned_power, ...
                    lfp_tfa_cfg, plottitle, result_file);
            end

        end
        % save mat file for site
        site_pow = sites_pow(i);
        save(fullfile(site_results_folder, ...
            ['LFP_Power_' site_pow.site_ID '.mat']), 'site_pow');

        session_pow.sites(i) = site_pow;
    end
        
    
    
    % Calculate average power spectrum across all sites
    session_avg = struct();
    % targets for this session
    targets = unique({states_lfp.target});
    for t = 1:length(targets)
        session_avg(t).target = targets{t};
        for cn = 1:length(site_conditions)
            %nsites = sum(contains({states_lfp.target}, site_conditions(cn).target));
            session_avg(t).condition(cn).hs_tuned_power = [];
            % variable to store no:of sites with trials satisfying this
            % condition
            isite = 0;            
            for i = 1:length(states_lfp)
                if ~strcmp(states_lfp(i).target, targets{t})
                    continue;
                end

                % calculate the average LFP power spectrum across sites for this condition 
                if ~isempty(sites_pow(i).condition(cn).hs_tuned_power) && ...
                        isfield(sites_pow(i).condition(cn).hs_tuned_power, 'mean')
                    isite = isite + 1;
                    % struct to store average evoked LFP across sites
                    %cond_based_psd(cn).hs_tuned_power = struct();%cell(length(analyse_states),length(hs_labels));            


                    for hs = 1:length(hs_labels)
                        for ep = 1:length(lfp_tfa_cfg.epochs)
                            if ~isempty(sites_pow(i).condition(cn).hs_tuned_power(ep, hs).mean)
                                
                                if isite == 1
                                    session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean = ...
                                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).mean ;
                                    session_avg(t).condition(cn).hs_tuned_power(ep, hs).freq = ...
                                        sites_pow(i).condition(cn).hs_tuned_power(ep, hs).freq;
                                else
                                    if ~isempty(session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean)
                                        nfreqs = length(session_avg(t).condition(cn).hs_tuned_power(ep, hs).freq);
                                        if nfreqs > length(sites_pow(i).condition(cn).hs_tuned_power(ep, hs).freq)
                                            nfreqs = length(sites_pow(i).condition(cn).hs_tuned_power(ep, hs).freq);
                                        end                               
                                        session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean = ...
                                            session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean(1:nfreqs) + ...
                                            sites_pow(i).condition(cn).hs_tuned_power(ep, hs).mean(1:nfreqs) ;
                                    else
                                        session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean = ...
                                            sites_pow(i).condition(cn).hs_tuned_power(ep, hs).mean ;
                                    end
                                end
                                session_avg(t).condition(cn).hs_tuned_power(ep, hs).hs_label = ...
                                    sites_pow(i).condition(cn).hs_tuned_power(ep, hs).hs_label;
                                session_avg(t).condition(cn).hs_tuned_power(ep, hs).epoch_name = ...
                                    sites_pow(i).condition(cn).hs_tuned_power(ep, hs).epoch_name;
                                session_avg(t).condition(cn).condition = site_conditions(cn);
                                session_avg(t).condition(cn).session = states_lfp(i).session;
                                session_avg(t).condition(cn).target = states_lfp(i).target;
                                session_avg(t).condition(cn).condition = site_conditions(cn);
                                session_avg(t).condition(cn).label = site_conditions(cn).label;
                            end
                        end
                    end
                end
            end
            
            % average TFR across sites for a session
            if isfield(session_avg(t).condition(cn).hs_tuned_power, 'mean')
                for hs = 1:length(hs_labels)
                    for ep = 1:length(lfp_tfa_cfg.epochs)
                        session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean = ...
                            session_avg(t).condition(cn).hs_tuned_power(ep, hs).mean / isite;
                        session_avg(t).condition(cn).hs_tuned_tfs(ep, hs).nsites = isite;
                    end
                end
            end 
            
            % plot average power spectrum across sites for this session
            if ~isempty(session_avg(t).condition(cn).hs_tuned_power)
                plottitle = ['Session: ', session_avg(t).condition(cn).session ...
                    ', Target = ' session_avg(t).condition(cn).target ', '  ...
                    'Perturb ' num2str(site_conditions(cn).perturbation_group{1}) ', '];
                if site_conditions(cn).choice == 0
                    plottitle = [plottitle 'Instructed trials'];
                else
                    plottitle = [plottitle 'Choice trials'];
                end
                results_file = fullfile(results_folder_psd, ...
                    ['LFP_Power_' session_avg(t).condition(cn).session '_' site_conditions(cn).label '.png']);
                lfp_tfa_plot_hs_tuned_psd(session_avg(t).condition(cn).hs_tuned_power, ...
                            lfp_tfa_cfg, plottitle, results_file);
            end
            
        end        
    end
    
    close all;
    
    session_pow.session_avg = session_avg;
    
    % save mat files
    save(fullfile(results_folder_psd, ['LFP_Power_' session_pow.session '.mat']), 'session_pow');
    % save settings file
    save(fullfile(results_folder_psd, 'lfp_tfa_settings.mat'), 'lfp_tfa_cfg');
end
        
