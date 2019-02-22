function [ evoked_lfp ] = evokedLFPAnalysis( ft_data_sites )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evoked LFP Analysis Script
    % This script calculates the mean and standard deviation of LFP data from
    % all trials per site and also average across sites
    % This script takes the structure created in Prepare_FT_Datastruct_persite
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %close all;
    % create a struct to store evoked LFP analysis results
    evoked_lfp = struct();

    % folder to save figures
    fig_folder = 'C:\Data\MIP_timefreq_analysis\Figures\';
    
    
    sessionName = ft_data_sites(1).session;
    fig_folder_evoked = fullfile(fig_folder, sessionName, date, 'LFP_evoked_analysis');
    mkdir(fig_folder_evoked);


    % cell array to store handspace tuned lfp data of all sites
    grand_lfp_arr = cell(1,length([ft_data_sites(1).states(1).hndspc.label]));
    % loop through each site
    for i = 1:length(ft_data_sites)
        figure; % to plot mean and std of LFP for each site
        % get site info
        evoked_lfp(i).siteID = ft_data_sites(i).siteID;
        evoked_lfp(i).target = ft_data_sites(i).target;
        evoked_lfp(i).session = ft_data_sites(i).session;
        evoked_lfp(i).nblocks = ft_data_sites(i).nblocks;
        % loop through each handspace label
        hs_labels = [ft_data_sites(1).states(1).hndspc.label];
        for hs = 1:length(hs_labels)
            evoked_lfp(i).label(hs).name = hs_labels(hs);
            % loop through each state
            analyse_states = {ft_data_sites(1).states  .name};
            lfp_arr = [];
            lfp_time = [];
            state_start_s = zeros(1,length(analyse_states)); % start sample of each state
            state_end_s = zeros(1,length(analyse_states)); % end sample of each state
            for st = 1:length(analyse_states)
                evoked_lfp(i).label(hs).state(st).name = analyse_states(st);
                state_start_s(st) = size(lfp_arr, 2);
                % concatenate all states
                lfp_arr = [lfp_arr cell2mat(ft_data_sites(i).states(st).hndspc(hs).trial)];
                lfp_time = [lfp_time ft_data_sites(i).states(st).hndspc(hs).time{1}];
                state_end_s(st) = size(lfp_arr, 2);
            end
            % accumulate sites
            grand_lfp_arr{hs} = [grand_lfp_arr{hs}; lfp_arr];
            state_onset_s = find(lfp_time == 0);
            evoked_lfp(i).label(hs).ntrials = size(lfp_arr, 1);
            % save state info in evoked LFP struct
            evoked_lfp(i).label(hs).state(st).start_s = state_start_s(st);
            evoked_lfp(i).label(hs).state(st).onset_s = state_onset_s(st);
            evoked_lfp(i).label(hs).state(st).end_s = state_end_s(st);
            % find mean and standard deviation
            lfp_mean = mean(lfp_arr, 1);
            evoked_lfp(i).label(hs).lfp_mean = lfp_mean;
            % plot mean
            subplot(2,2,hs)
            plot(lfp_mean)

            std_lfp = std(lfp_arr, 1);
            evoked_lfp(i).label(hs).std_LFP = std_lfp;
            % plot std
            hold on; plot(lfp_mean + std_lfp, 'r', 'linestyle', '--');
            plot(lfp_mean - std_lfp, 'r', 'linestyle', '--');

            % mark states
            for st = 1:length(analyse_states)
                line([state_onset_s(st) state_onset_s(st)], ylim, 'color', 'k');
                line([state_end_s(st) state_end_s(st)], ylim, 'color', 'k', 'linewidth',3);
                line([state_onset_s(st) state_onset_s(st)], ylim, 'color', 'k');
                line([state_end_s(st) state_end_s(st)], ylim, 'color', 'k', 'linewidth',3);
                % mark state onsets
                set(gca,'xtick',state_onset_s)
                set(gca,'xticklabels',strrep(analyse_states, '_', '\_'))
            end  
            title(hs_labels(hs));
        end
        plottitle = ['Session: ' evoked_lfp(i).session ', Target = ' evoked_lfp(i).target ', Site ID: ' evoked_lfp(i).siteID ...
            '(nblocks = ' num2str(evoked_lfp(i).nblocks) ')'];
        ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
            , 'EdgeColor', 'None', 'HorizontalAlignment', 'center');
        saveas (gca, fullfile(fig_folder_evoked, [evoked_lfp(i).siteID '_evoked_LFP.png']))
    end
    save evoked_LFP_analysis evoked_lfp;

    % Average across sites
    nblocks = 0;
    nsites = numel(evoked_lfp);
    for site = evoked_lfp
        nblocks = site.nblocks + nblocks;
    end
    figure;
    for hs = 1:length(site.label)
        lfp_mean = mean(grand_lfp_arr{hs}, 1);
        lfp_std = std(grand_lfp_arr{hs}, 1);
        % plot
        % plot mean
        subplot(2,2,hs)
        plot(lfp_mean)

        % plot std
        hold on; plot(lfp_mean + std_lfp, 'r', 'linestyle', '--');
        plot(lfp_mean - lfp_std, 'r', 'linestyle', '--');

        % mark states
        for st = 1:length(analyse_states)
            line([state_onset_s(st) state_onset_s(st)], ylim, 'color', 'k');
            line([state_end_s(st) state_end_s(st)], ylim, 'color', 'k', 'linewidth',3);
            line([state_onset_s(st) state_onset_s(st)], ylim, 'color', 'k');
            line([state_end_s(st) state_end_s(st)], ylim, 'color', 'k', 'linewidth',3);
            % mark state onsets
            set(gca,'xtick',state_onset_s)
            set(gca,'xticklabels',strrep(analyse_states, '_', '\_'))
        end  
        title(hs_labels(hs));
    end  
    plottitle = ['Session: ' evoked_lfp(1).session ' (nsites = ' num2str(nsites) ', nblocks = ' num2str(evoked_lfp(1).nblocks) ')'];
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'None', 'HorizontalAlignment', 'center');
    saveas (gca, fullfile(fig_folder_evoked, [evoked_lfp(i).session '_evoked_LFP.png']))
end

