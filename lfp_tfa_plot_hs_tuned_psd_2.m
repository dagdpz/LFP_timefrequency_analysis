function lfp_tfa_plot_hs_tuned_psd_2( avg_lfp_psd, lfp_tfa_cfg, plottitle, results_file )
%lfp_tfa_plot_hs_tuned_psd_2  - Plots the LFP power spectrum 
%averages for different hand-space conditions to be compared
%
% USAGE:
%   lfp_tfa_plot_hs_tuned_psd_2( avg_lfp_psd, lfp_tfa_cfg, plottitle, results_file )
%
% INPUTS:
%       avg_lfp_psd      - average LFP power spectrum for different
%       hand-space conditions to be compared
%		lfp_tfa_cfg      - struct containing the required settings
%           Required Fields: see settings/lfp_tfa_settings_example
%               1. compare.reach_hands  - reach hands to be included
%               2. compare.reach_spaces - reach spaces to be included
%               3. analyse_epochs - epochs to be analysed
%       plottitle       - title for the plot
%       results_file    - path to filename to store the resulting image
%
% REQUIRES:	export_fig
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_plot_site_powspctrm, 
% lfp_tfa_avg_pow_across_sites, lfp_tfa_avg_pow_across_sessions, 
% export_fig
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
    
    cm = colormap(othercolor('Cat_12', size(avg_lfp_psd, 1)));
    %cm = colormap(jet(size(avg_lfp_psd, 1)));
    
    % number of subplots required
    nhandlabels = length(lfp_tfa_cfg.compare.reach_hands);
    nspacelabels = length(lfp_tfa_cfg.compare.reach_spaces);
    
    % loop through handspace
    for hs = 1:size(avg_lfp_psd, 2)
        if ~isempty([avg_lfp_psd(:,hs).mean] )       
            subplot(nhandlabels,nspacelabels,hs)
            for ep = 1:size(avg_lfp_psd, 1)
                % check if no trials exist for this epoch and HS
                if ~isempty(avg_lfp_psd(ep, hs).mean)


                    % now plot

                    hold on;

                    plot(1:numel(avg_lfp_psd(ep, hs).freq), 10*log10(avg_lfp_psd(ep, hs).mean), ...
                        'Color', cm(ep,:), 'LineWidth', 2);                


                end
                
                % set y-axis (power) limits in dB
                set(gca, 'ylim', [-100 -40]);
                                
                set(gca, 'box', 'on');                
                ylabel('LFP Power (dB)');
                xlabel('Frequency (Hz)');
                subplottitle = avg_lfp_psd(ep, hs).hs_label{1};
                if isfield(avg_lfp_psd(1, hs), 'nsessions')
                    subplottitle = [subplottitle ' (nsessions = ' num2str(avg_lfp_psd(1, hs).nsessions) ')'];
                elseif isfield(avg_lfp_psd(1, hs), 'trials')
                    subplottitle = [subplottitle ' (ntrials = ' num2str(length(avg_lfp_psd(1, hs).trials)) ')'];
                elseif isfield(avg_lfp_psd(1, hs), 'nsites')
                    subplottitle = [subplottitle ' (nsites = ' num2str(avg_lfp_psd(1, hs).nsites) ')'];
                end
            end
            title(subplottitle); 
            legend({lfp_tfa_cfg.analyse_epochs{:,2}});
            
            % freq bands
            fbandstart = [2, 4, 8, 12, 18, 32, 80];
            fbandstart_idx = zeros(size(fbandstart));
            for f = fbandstart
                f_idx = find(abs(avg_lfp_psd(ep, hs).freq - f) == ...
                    min(abs(avg_lfp_psd(ep, hs).freq - f)), 1, 'first');
                line([f_idx f_idx], ylim,'color', 'k', 'linestyle', '--');
                fbandstart_idx(fbandstart == f) = f_idx;
            end
            % x axis ticks
            set(gca, 'xtick', fbandstart_idx);%[1:8:numel(avg_lfp_psd(ep, hs).freq)]);
            set(gca, 'xticklabel', fbandstart);%...
                %round((avg_lfp_psd(ep, hs).freq([1:8:numel(avg_lfp_psd(ep, hs).freq)]))), ...
                %'fontsize', 8);
            set(gca, 'xlim', [1 numel(avg_lfp_psd(ep, hs).freq)]);
                
        end
    end
        
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    export_fig(h, results_file);

end

