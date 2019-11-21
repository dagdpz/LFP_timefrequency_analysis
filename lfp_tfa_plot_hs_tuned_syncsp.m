function lfp_tfa_plot_hs_tuned_syncsp( avg_hs_tuned_syncsp, lfp_tfa_cfg, plottitle, results_file )
%lfp_tfa_plot_hs_tuned_syncsp  - Plots the LFP-LFP phase sync spectrum 
%averages for different hand-space conditions to be compared
%
% USAGE:
%   lfp_tfa_plot_hs_tuned_syncsp( avg_hs_tuned_syncsp, lfp_tfa_cfg, ...
%       plottitle, results_file )
%
% INPUTS:
%       avg_hs_tuned_syncsp     - average LFP-LFP sync spectrum for 
%       different hand-space conditions to be compared
%		lfp_tfa_cfg      - struct containing the required settings
%           Required Fields: see settings/lfp_tfa_settings_example
%               1. analyse_epochs - epochs to be analysed
%               2. compare.reach_hands          - hand labels to compare
%               3. compare.reach_spaces         - space labels to compare
%       plottitle       - title for the plot
%       results_file    - path to filename to store the resulting image
%
% REQUIRES:	export_fig
%
% See also settings/lfp_tfa_settings_example, 
% lfp_tfa_sitepair_averaged_syncspctrm, colormap, export_fig
% lfp_tfa_avg_sessions_syncspctrm, lfp_tfa_avg_sitepairs_syncspctrm
%
% Author(s):	S.Nair, DAG, DPZ
% URL:		http://www.dpz.eu/dag
%
% Change log:
% 2019-02-15:	Created function (Sarath Nair)
% 2019-03-05:	First Revision
% ...
% $Revision: 1.0 $  $Date: 2019-03-05 17:18:00 $

% ADDITIONAL INFO:
% ...
%%%%%%%%%%%%%%%%%%%%%%%%%[DAG mfile header version 1]%%%%%%%%%%%%%%%%%%%%%%%%%
    
    h = figure;
    set(h, 'position', [100, 100,900, 675]);
    
    %cm = colormap(othercolor('Cat_12', size(avg_hs_tuned_syncsp, 1)));
    % color scheme for plotting spectra for each epoch
    if isfield(lfp_tfa_cfg, 'epoch_colors') && ...
            length(lfp_tfa_cfg.epoch_colors) == size(avg_hs_tuned_syncsp, 1)
        cm = colormap(lfp_tfa_cfg.epoch_colors);
    else % default color scheme
        cm = colormap(othercolor('Cat_12', size(avg_hs_tuned_syncsp, 1)));
    end
    %cm = colormap(jet(size(avg_hs_tuned_syncsp, 1)));
    
    % number of subplots required
    nhandlabels = length(lfp_tfa_cfg.compare.reach_hands);
    nspacelabels = length(lfp_tfa_cfg.compare.reach_spaces);
       
    % loop through handspace
    for hs = 1:size(avg_hs_tuned_syncsp, 2)
        % check if no trials exist for this condition and HS
        if isfield(cat(2, avg_hs_tuned_syncsp(:, hs).ppc), 'ppcspctrm')
            ppc =  cat(2, avg_hs_tuned_syncsp(:, hs).ppc);
        else
            continue;
        end   
        if ~isempty(cat(2, ppc.ppcspctrm))
            subplot(nhandlabels, nspacelabels, hs)
            for ep = 1:size(avg_hs_tuned_syncsp, 1)
                % check if no trials exist for this epoch and HS
                if ~isempty(avg_hs_tuned_syncsp(ep, hs).ppc.ppcspctrm) && ...
                        ~isempty(avg_hs_tuned_syncsp(ep, hs).ppc.freq)


                    % now plot

                    hold on;

                    plot(1:numel(avg_hs_tuned_syncsp(ep, hs).ppc.freq), avg_hs_tuned_syncsp(ep, hs).ppc.ppcspctrm, ...
                        'Color', cm(ep,:), 'LineWidth', 2);                


                end
                
                % set y-axis (power) limits
                set(gca, 'ylim', [0 1]);
                
                %...
                    %round((avg_hs_tuned_syncsp(ep, hs).ppc.freq([1:8:numel(avg_hs_tuned_syncsp(ep, hs).ppc.freq)]))), ...
                    %'fontsize', 8);
%                 set(gca, 'xtick', (lfp_tfa_cfg.tfr.foi([1:8:numel(lfp_tfa_cfg.tfr.foi)])));
%                 set(gca, 'xticklabel', ...
%                     round((lfp_tfa_cfg.tfr.foi([1:8:numel(lfp_tfa_cfg.tfr.foi)]))), ...
%                     'fontsize', 8);
                set(gca, 'xlim', [1 numel(avg_hs_tuned_syncsp(ep, hs).ppc.freq)]);
                %set(gca, 'xlim', [lfp_tfa_cfg.tfr.foi(1) lfp_tfa_cfg.tfr.foi(end)]);
                set(gca, 'box', 'on');                
                ylabel('PPC');
                xlabel('Frequency (Hz)');
                subplottitle = avg_hs_tuned_syncsp(ep, hs).hs_label{1};
                if isfield(avg_hs_tuned_syncsp(1, hs), 'nsessions')
                    subplottitle = [subplottitle ' (nsessions = ' num2str(avg_hs_tuned_syncsp(1, hs).nsessions) ')'];
                end
                if isfield(avg_hs_tuned_syncsp(1, hs), 'trials')
                    subplottitle = [subplottitle ' (ntrials = ' num2str(length(avg_hs_tuned_syncsp(1, hs).trials)) ')'];
                end
                if isfield(avg_hs_tuned_syncsp(1, hs), 'nsites')
                    subplottitle = [subplottitle ' (nsites = ' num2str(avg_hs_tuned_syncsp(1, hs).nsites) ')'];
                end
            end
            title(subplottitle); 
            legend({lfp_tfa_cfg.analyse_epochs{:,2}});
            
            fbandstart = [2, 4, 8, 12, 18, 32, 80];
            fbandstart_idx = zeros(size(fbandstart));
            for f = fbandstart
                f_idx = find(abs(avg_hs_tuned_syncsp(ep, hs).ppc.freq - f) == ...
                    min(abs(avg_hs_tuned_syncsp(ep, hs).ppc.freq - f)), 1, 'first');
                line([f_idx f_idx], ylim,'color', 'k', 'linestyle', '--');
                fbandstart_idx(fbandstart == f) = f_idx;
            end

            % log y axis ticks
            set(gca, 'xtick', fbandstart_idx);%[1:8:numel(avg_hs_tuned_syncsp(ep, hs).ppc.freq)]);
            set(gca, 'xticklabel', fbandstart);
            
        end
    end
        
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    export_fig(h, results_file);

end

