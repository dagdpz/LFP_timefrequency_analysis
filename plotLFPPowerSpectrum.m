function [ lfp_spectrum ] = plotLFPPowerSpectrum( ft_data_sites )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plotLFPPowerSpectrum.m plots power spectrum for LFP per site and
    % average acorss sites for 4 conditions
    % This script calculates the mean and standard deviation of LFP data from
    % all trials per site and also average across sites
    % This script takes the structure created in Prepare_FT_Datastruct_persite
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    lfp_spectrum = struct();

    % folder to save figures
    fig_folder = 'C:\Data\MIP_timefreq_analysis\Figures\';
    
    
    sessionName = ft_data_sites(1).session;
    fig_folder_psd = fullfile(fig_folder, sessionName, date, 'LFP_powerspectrum');
    mkdir(fig_folder_psd);

    % cell array to store handspace tuned lfp data of all sites
    %session_lfp_spectrum = cell(1,length([ft_data_sites(1).states(1).hndspc.label]));
    % loop through each site
    for i = 1:length(ft_data_sites)
        figure; % to plot mean and std of LFP for each site
        % get site info
        lfp_spectrum(i).siteID = ft_data_sites(i).siteID;
        lfp_spectrum(i).target = ft_data_sites(i).target;
        lfp_spectrum(i).session = ft_data_sites(i).session;
        lfp_spectrum(i).nblocks = ft_data_sites(i).nblocks;
        % loop through each handspace label
        hs_labels = [ft_data_sites(1).states(1).hndspc.label];
        for hs = 1:length(hs_labels)
            lfp_spectrum(i).label(hs).name = hs_labels(hs);
            % loop through each state
            analyse_states = {ft_data_sites(1).states  .name};
            lfp_arr = [];
            lfp_time = [];
            state_end_s = zeros(1,length(analyse_states)); % end sample of each state
            for st = 1:length(analyse_states)
                lfp_spectrum(i).label(hs).state(st).name = analyse_states(st);
                % concatenate all trials
                lfp_arr = [lfp_arr cell2mat(ft_data_sites(i).states(st).hndspc(hs).trial)];
                lfp_time = [lfp_time ft_data_sites(i).states(st).hndspc(hs).time{1}];
                lfp_fs = ft_data_sites(i).states(st).hndspc(hs).fsample;
                state_end_s(st) = size(lfp_arr, 2);
            end
            % accumulate sites
            % grand_lfp_arr{hs} = [grand_lfp_arr{hs}; lfp_arr];
            lfp_spectrum(i).label(hs).ntrials = size(lfp_arr, 1);
            % array to store mean spectrum for each condition
            lfp_spectrum(i).label(hs).mean_spectrum = zeros(1, size(lfp_arr, 2)/2 + 1);
            % find mean lfp spectrum across all trials
            for t = 1:size(lfp_arr, 1)
                [psd, f] = powerspectrum(lfp_arr(t,:), lfp_fs);
                lfp_spectrum(i).label(hs).mean_spectrum = ...
                    lfp_spectrum(i).label(hs).mean_spectrum + psd;
            end
            lfp_spectrum(i).label(hs).mean_spectrum = ...
                lfp_spectrum(i).label(hs).mean_spectrum / lfp_spectrum(i).label(hs).ntrials;
            lfp_spectrum(i).label(hs).frequency = f;
            
            % plot mean power spectrum
            subplot(2,2,hs)
            plot((lfp_spectrum(i).label(hs).frequency), ...
                10*log10(lfp_spectrum(i).label(hs).mean_spectrum));
            xlabel('Frequency (Hz)');
            ylabel('Power (dB)');
            set(gca, 'XScale', 'log')
            set(gca, 'xlim', [1e-1, 1e3]);
            
            grid on;

            title(hs_labels(hs));
        end
        plottitle = ['Session: ' lfp_spectrum(i).session ', Target = ' lfp_spectrum(i).target ', Site ID: ' lfp_spectrum(i).siteID ...
            '(nblocks = ' num2str(lfp_spectrum(i).nblocks) ')'];
        ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
            , 'EdgeColor', 'None', 'HorizontalAlignment', 'center');
        saveas (gca, fullfile(fig_folder, [lfp_spectrum(i).siteID '_spectrum.png']))
    end
    save lfp_spectrum lfp_spectrum;

    % Average across sites
    figure;
    for site = lfp_spectrum
        for hs = 1:length(site.label)
            session_mean_psd = mean(site.label(hs).mean_spectrum, 1);
            frequency = site.label(hs).frequency;
            % plot
            % plot mean
            subplot(2,2,hs)
            plot((frequency), 10*log10(session_mean_psd));
            xlabel ('Frequency (Hz)');
            ylabel ('Power (dB)');
            set(gca, 'XScale', 'log')
            set(gca, 'xlim', [1e-1, 1e3]);
            
            grid on;
            
            title(hs_labels(hs));
        end  
    end
    plottitle = ['Session: ' lfp_spectrum(1).session ' (nsites = ' num2str(length(lfp_spectrum)) ', nblocks = ' num2str(lfp_spectrum(1).nblocks) ')'];
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'None', 'HorizontalAlignment', 'center');
    saveas (gca, fullfile(fig_folder_psd, [lfp_spectrum(i).session '_spectrum.png']))
end

