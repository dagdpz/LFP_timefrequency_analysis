function [ ft_data_sites, TFR_avg ] = computeAndPlotTFR( ft_data_sites, baseline ) 
%, choice, inactivation, analyse_states, all_states )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This script takes the data structure with site wise and trial wise 
    % arranged LFP data and creates a data structure which is compatible with 
    % the FT_DATATYPE_RAW for time frequency analysis of LFP
    % The data structure will be arranged as follows:
    % Sites -> Site (site ID, site name, states) -> State (state name,
    % handspace) -> handspace (label, trial (LFP), time)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %close all;
    % folder to save figures
    fig_folder = 'C:\Data\MIP_timefreq_analysis\Figures\';
    
    
    sessionName = ft_data_sites(1).session;
    fig_folder_tfr = fullfile(fig_folder, sessionName, date, 'handspace_LFP_timefreq_analysis');
    mkdir(fig_folder_tfr);
    

    % define the states to be analysed, should have same name as defined in
    % all_states data structure
    analyse_states = {ft_data_sites(1).states.name};
    % hand-space tuning of LFP
    hs_label = {ft_data_sites(1).states(1).hndspc.label};
    
%     % whether to consider choice trial or instructed for analysis (0 =
%     % instructed only, 2 = choice only, else both)
%     choice = 0;
%     % whether to choose control or inactivation trials
%     % 1 = only control trials, 0 = only inactivation trials
%     control = 0; 

    % cell array to store time frequency average across sites
    TFR_avg = cell(length(analyse_states),length(hs_label));
    
    % loop through each site
    nsites = length(ft_data_sites);

%     % create a data structure for time-freq analysis per site
%     ft_data_sites = struct();
% 

    
    %ft_data_allsites.nsites = ft_data_allsites.nsites + nsites;
%     for i = 1:length(session_lfp)
%         site = session_lfp(i);
%         ft_data_sites(i).session = site.session;
%         ft_data_sites(i).target = site.target;
%         ft_data_sites(i).siteID = site.site_ID;
%         ft_data_sites(i).states = struct();
%         ft_data_sites(i).TFR_hs = cell(length(analyse_states),length(hs_label));
%         for st = 1:length(analyse_states)
%             % store name of the state
%             ft_data_sites(i).states(st).name = analyse_states(st);%all_states([all_states.state_ID] == analyse_states{st}).state_name;
%             % struct to save hand-space tuned FT type LFP data for this state
%             ft_data_sites(i).states(st).hndspc = struct();
%             % initialize hand-space tuning structures
%             for hs = 1:length( hs_label)
%                 ft_data_sites(i).states(st).hndspc(hs).label = hs_label(hs);
%                 ft_data_sites(i).states(st).hndspc(hs).trial = {};
%                 ft_data_sites(i).states(st).hndspc(hs).time  = {};
%             end
%         end
%         % loop through each trial for this site
%         blocks = [];
%         ntrials = 0;
%         for t = 1:length(site.trials)
%             trial = site.trials(t);
%             % whether this trial should be considered
%             consider_trial = ~trial.noisy; % if trial is not noisy 
%             % consider a trial based on whether it is instructed or choice
%             % trial
%             if choice == 1 % user wants to select choice trials
%                 consider_trial = consider_trial && trial.choice_trial; % depending on whether trial is instructed or not
%             else % user wants to select instructed trials
%                 consider_trial = consider_trial && trial.choice_trial == 0;
%             end
%             % consider a trial based on whether it is control or
%             % inactivation trial
% %             if inactivation == 1 % user wants to select inactivation trials
% %                 consider_trial = consider_trial && trial.perturbation; 
% %             else
% %                 consider_trial = consider_trial && trial.perturbation == 0; 
% %             end
%             % if trial should be considered
%             if consider_trial%~trial.noisy && ~trial.choice_trial && 
%                 ntrials = ntrials + 1;
%                 blocks = [blocks, trial.block];
%                 hs_tuning = trial.hndspc_lbl;
%                 hs_idx = find(strcmp(hs_label, hs_tuning));
% 
%                 % loop through states
%                 for st = 1:length(trial.states)
%                     state = trial.states(st);
%                     st_idx = find(strcmp(analyse_states, state.name));
%                     ft_data_sites(i).states(st_idx).hndspc(hs_idx).trial = ...
%                         [ft_data_sites(i).states(st_idx).hndspc(hs_idx).trial; ...
%                         trial.lfp_data(state.start_s:state.end_s)];
%                     ft_data_sites(i).states(st_idx).hndspc(hs_idx).time = ...
%                         [ft_data_sites(i).states(st_idx).hndspc(hs_idx).time; ...
%                         trial.time(state.start_s:state.end_s) - trial.time(state.onset_s)];
%                     ft_data_sites(i).states(st_idx).hndspc(hs_idx).fsample = trial.fsample;
%                 end            
%             end
%         end
%         ft_data_sites(i).nblocks = length(unique(blocks));
%         ft_data_sites(i).ntrials = ntrials;
        
        % per-site time frequency analysis
        % make a folder for saving figures
        mkdir(fig_folder);
        
    for i = 1:length(ft_data_sites)
        figure; colormap jet;%
        % make a struct for concatenating TFR for all states
        concat_states_TFR = struct();
        % loop through each hand-space tuning
        for hs = 1:length(hs_label)
            % baseline for this handspace tuning
            data_lfp_baseline = ft_data_sites(i).baseline.hndspc(hs);
            % configuration for TFR of baseline
            if ~isempty(data_lfp_baseline.trial)
                nsample = length(data_lfp_baseline.time{1});
                fs = data_lfp_baseline.fsample;
                ts = 1/fs;
                tstart = data_lfp_baseline.time{1}(1);
                tend = data_lfp_baseline.time{1}(end);
                cfg              = [];
                cfg.output       = 'pow';
                cfg.method       = 'mtmconvol';
                cfg.taper        = 'hanning';
                cfg.pad          = 'nextpow2';
                cfg.foi          = fs/500:fs/500:fs/10;             % analysis 2 to 100 Hz in steps of 2 Hz
                cfg.t_ftimwin    = ones(length(cfg.foi),1).*200*ts;    % length of time window = 0.2 sec
                cfg.toi          = linspace(tstart, tend, nsample);% time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
                cfg.channel      = data_lfp_baseline.label;
                disp(['Processing ' cfg.channel]);
                TFR_baseline     = ft_freqanalysis(cfg, data_lfp_baseline);
                % get mean value of baseline power for each frequency
                baseline_mean = nanmean(TFR_baseline.powspctrm, 3);
            end
            % loop through each state
            for st = 1:length(ft_data_sites(i).states)
                % now this is a FT_DATATYPE_RAW on which FT time frequency
                % analzsis can be done
                data_lfp_st_hs = ft_data_sites(i).states(st).hndspc(hs);
                nsample = length(ft_data_sites(i).states(st).hndspc(hs).time{1});
                if ~isempty(data_lfp_st_hs.trial)
                    fs = data_lfp_st_hs.fsample;
                    ts = 1/fs;
                    tstart = data_lfp_st_hs.time{1}(1);
                    tend = data_lfp_st_hs.time{1}(end);
                    % configuration for TFR
                    cfg              = [];
                    cfg.output       = 'pow';
                    cfg.method       = 'mtmconvol';
                    cfg.taper        = 'hanning';
                    cfg.pad          = 'nextpow2';
                    cfg.foi          = fs/500:fs/500:fs/10;             % analysis 2 to 100 Hz in steps of 2 Hz
                    cfg.t_ftimwin    = ones(length(cfg.foi),1).*200*ts;    % length of time window = 0.2 sec
                    cfg.toi          = linspace(tstart, tend, nsample);% time window "slides" from -0.5 to 1.5 sec in steps of 0.05 sec (50 ms)
                    cfg.channel      = data_lfp_st_hs.label;
                    disp(['Processing ' cfg.channel]);
                    TFRhann_fixed = ft_freqanalysis(cfg, data_lfp_st_hs);    
                    % subtract baseline mean for each frequency
                    TFRhann_fixed.powspctrm = TFRhann_fixed.powspctrm - ...
                        repmat(baseline_mean, [1 1 size(TFRhann_fixed.powspctrm, 3)]);
                    % concatenate TFRs for each state
                    if st == 1
                        concat_states_TFR = TFRhann_fixed;
                    else
                        % concatenate times
                        % add an arbitrary large number for concatenated trials
                        concat_states_TFR.time = [concat_states_TFR.time, TFRhann_fixed.time];
                        concat_states_TFR.powspctrm = cat(3, concat_states_TFR.powspctrm, TFRhann_fixed.powspctrm);
                    end

                    % now save the TFR to corresponding hand space tuning
                    ft_data_sites(i).TFR{st, hs} = TFRhann_fixed;

                    % accumulate TFR per site for averaging across sites
                    if i == 1
                        TFR_avg{st, hs} = TFRhann_fixed;
                        TFR_avg{st, hs}.powspctrm = (1/nsites) * TFRhann_fixed.powspctrm;
                        TFR_avg{st, hs}.ntrials = length(data_lfp_st_hs.trial);
                    else
                        TFR_avg{st, hs}.powspctrm = TFR_avg{st, hs}.powspctrm + (1/nsites) * TFRhann_fixed.powspctrm;
                        TFR_avg{st, hs}.ntrials = TFR_avg{st, hs}.ntrials + length(data_lfp_st_hs.trial);
                    end
                    TFR_avg{st, hs}.label = strcat(analyse_states(st),hs_label(hs));
                else
                    ft_data_sites(i).TFR{st, hs} = [];
                end
            
            end 

            % Visualization
            if ~isempty(concat_states_TFR) > 0
                % normalize power spectrogram by z-score
%                 concat_states_TFR.powspctrm = (concat_states_TFR.powspctrm - ...
%                     nanmean(concat_states_TFR.powspctrm(:))) / nanstd(concat_states_TFR.powspctrm(:));
                time_to_event = concat_states_TFR.time;
                tshift = time_to_event(1);
                % first rearrange the time axis of TFR
                bin_t = concat_states_TFR.time(2) - concat_states_TFR.time(1);
                concat_states_TFR.time = 0:bin_t:bin_t*(length(concat_states_TFR.time)-1);
                % now find the state onsets in the new time axis
                rearr_state_onset_t = concat_states_TFR.time(time_to_event == 0);
                %600*ts%0:1:length(concat_states_TFR.time)-1;
                % rearrange the baseline
                baseline = [-0.4 -0.1];
                baseline_shift = baseline - tshift;

                cfg = [];
                cfg.baseline     = 'no'; %baseline_shift;                 % -400ms to -100ms before the onset of first state
                cfg.baselinetype = 'absolute';  
                cfg.maskstyle    = 'saturation';
                cfg.zlim         = [-1.5e-9 1.5e-9];
                cfg.masknans     = 'yes';
                cfg.channel  = concat_states_TFR.label;
                subplot(2,2,hs)
                ft_singleplotTFR(cfg, concat_states_TFR);
                % mark state onsets
                set(gca,'xtick',rearr_state_onset_t)
                set(gca,'xticklabels',strrep(analyse_states, '_', '\_'))
                for t = rearr_state_onset_t
                    line([t t], ylim);
                end
                xlabel('Time');
                ylabel('Frequency');
                title(sprintf('%s (ntrials = %g)', cfg.channel{1}, length(data_lfp_st_hs.trial)));
                line([0 0], ylim, 'color', 'k');
            end
        end
        % plot title
    %     sgtitle(sprintf('Session: %s, Target = %s, Site: %s (nblocks = %g)', ft_data_sites(i).session, ...
    %         ft_data_sites(i).target, ft_data_sites(i).siteID, ft_data_sites(i).nblocks));
        plottitle = ['Session: ' ft_data_sites(i).session ', Target = ' ft_data_sites(i).target ', Site ID: ' ft_data_sites(i).siteID ...
            '(nblocks = ' num2str(ft_data_sites(i).nblocks) ')'];
        ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
            , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
        saveas(gca, fullfile(fig_folder_tfr, [ft_data_sites(i).siteID '.png']));
    end

    % save variables
    save ft_data_sites ft_data_sites;
    %save TFR_avg TFR_avg;

    %%% Visualiyation of time frequency average across site %%%
    figure; colormap jet;
    for hs = 1:length(hs_label)
        for st = 1:length(analyse_states)
            % first concatenate TFR for different states for one hs tuning for
            % visualization
            if st == 1
                concat_TFR_avg = TFR_avg{st,hs};
                concat_TFR_avg.label = hs_label{hs};
            else
                concat_TFR_avg.time = [concat_TFR_avg.time TFR_avg{st,hs}.time];
                concat_TFR_avg.powspctrm = cat(3, concat_TFR_avg.powspctrm, TFR_avg{st,hs}.powspctrm);
            end

        end

        % Visualization

        time_to_event = concat_TFR_avg.time;
        tshift = time_to_event(1);
        % first rearrange the time axis of TFR
        bin_t = concat_TFR_avg.time(2) - concat_TFR_avg.time(1);
        concat_TFR_avg.time = 0:bin_t:bin_t*(length(concat_TFR_avg.time)-1);
        % now find the state onsets in the new time axis
        rearr_state_onset_t = concat_TFR_avg.time(time_to_event == 0);
        % rearrange the baseline
        baseline = [-0.4 -0.1];
        baseline_shift = baseline - tshift;

        cfg = [];
        cfg.baseline     = 'no';%baseline_shift;                 % -400ms to -100ms before the state onset
        cfg.baselinetype = 'absolute';  
        cfg.maskstyle    = 'saturation';
        cfg.zlim         = [-1.5e-9 1.5e-9];
        cfg.masknans     = 'yes';
        %concat_TFR_avg.label = TFR_avg{st,hs}.label{1};
        cfg.channel  = concat_TFR_avg.label;
        subplot(2,2,hs)
        ft_singleplotTFR(cfg, concat_TFR_avg);
        % mark state onsets
        set(gca,'xtick',rearr_state_onset_t)
        set(gca,'xticklabels',strrep(analyse_states, '_', '\_'))
        for t = rearr_state_onset_t
            line([t t], ylim);
        end
        xlabel('Time');
        ylabel('Frequency');
        title(sprintf('%s', cfg.channel{1}));
        line([0 0], ylim, 'color', 'k');
    end

    % plot title
    plottitle = ['Session: ' ft_data_sites(i).session '(nsites = ', ...
        num2str(length(ft_data_sites)), 'nblocks = ' num2str(ft_data_sites(i).nblocks) ')'];
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    % sgtitle(sprintf('Session: %s, Target = %s, (nsites = %g)', ft_data_sites(1).session, ...
    %     ft_data_sites(1).target(1:end-2), length(ft_data_sites)));
    saveas(gca, fullfile(fig_folder_tfr, [ft_data_sites(1).session '.png']));

end

