function lfp_tfa_plot_hs_tuned_csd( avg_hs_tuned_syncsp, lfp_tfa_cfg, plottitle, results_file )
%lfp_tfa_plot_hs_tuned_tfr  - Plots the LFP power spectrum 
%averages for different hand-space conditions to be compared
%
% USAGE:
%   lfp_tfa_plot_hs_tuned_psd( avg_lfp_psd, lfp_tfa_cfg, plottitle, results_file )
%
% INPUTS:
%       avg_lfp_psd      - average LFP power spectrum for different
%       hand-space conditions to be compared
%		lfp_tfa_cfg      - struct containing the required settings
%           Required Fields: see lfp_tfa_settings
%               1. analyse_epochs - epochs to be analysed
%       plottitle       - title for the plot
%       results_file    - path to filename to store the resulting image
%
% REQUIRES:	
%
% See also lfp_tfa_settings, lfp_tfa_plot_site_powspctrm, 
% lfp_tfa_avg_pow_across_sites, lfp_tfa_avg_pow_across_sessions
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
    
    cm = colormap(othercolor('BrBG4', size(avg_hs_tuned_syncsp, 1)));
       
    % loop through handspace
    for hs = 1:size(avg_hs_tuned_syncsp, 2)
        % check if no trials exist for this condition and HS
        if isfield(cat(2, avg_hs_tuned_syncsp(:, hs).csd), 'crsspctrm_abs_mean')
            csd =  cat(2, avg_hs_tuned_syncsp(:, hs).csd);
        else
            continue;
        end   
        if ~isempty(cat(2, csd.crsspctrm_abs_mean))
            subplot(2,2,hs)
            for ep = 1:size(avg_hs_tuned_syncsp, 1)
                % check if no trials exist for this epoch and HS
                if ~isempty(avg_hs_tuned_syncsp(ep, hs).csd.crsspctrm_abs_mean) && ...
                        ~isempty(avg_hs_tuned_syncsp(ep, hs).csd.freq)


                    % now plot

                    
                    semilogx(avg_hs_tuned_syncsp(ep, hs).csd.freq, ...
                        10*log10(avg_hs_tuned_syncsp(ep, hs).csd.crsspctrm_abs_mean), ...
                        'Color', cm(ep,:));                


                end
                % weird error: see lfp_tfa_plot_hs_tuned_syncsp
                hold on;

                % log y axis ticks
                set(gca, 'xtick', round(lfp_tfa_cfg.tfr.foi([1:8:numel(lfp_tfa_cfg.tfr.foi)])));
                set(gca, 'xticklabel', ...
                    round((lfp_tfa_cfg.tfr.foi([1:8:numel(lfp_tfa_cfg.tfr.foi)]))), ...
                    'fontsize', 8);
                set(gca, 'xlim', [lfp_tfa_cfg.tfr.foi(1) lfp_tfa_cfg.tfr.foi(end)]);
                % set y-axis (power) limits in dB
                set(gca, 'ylim', [-100 -40]);
                set(gca, 'box', 'on');                
                ylabel('Power (dB)');
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
        end
    end
        
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    export_fig(h, results_file);

end

