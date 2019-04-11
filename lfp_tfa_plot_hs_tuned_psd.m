function lfp_tfa_plot_hs_tuned_psd( avg_lfp_psd, lfp_tfa_cfg, plottitle, results_file )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    figure;
    cm = colormap('jet');
    cm(1,:,:) = [1,1,1];
    colormap(cm);

    % loop through handspace
    for hs = 1:size(avg_lfp_psd, 2)
        subplot(2,2,hs)
        for ep = 1:size(avg_lfp_psd, 1)
            % check if no trials exist for this epoch and HS
            if ~isempty(avg_lfp_psd(ep, hs).mean)
            
            
                % now plot
                
                hold on;

                plot(avg_lfp_psd(ep, hs).freq, 10*log10(avg_lfp_psd(ep, hs).mean));                
                               
            
            end
            % log y axis ticks
            set(gca, 'xtick', (avg_lfp_psd(ep, hs).freq([1:8:numel(avg_lfp_psd(ep, hs).freq)])));
            set(gca, 'xticklabel', ...
                round((avg_lfp_psd(ep, hs).freq([1:8:numel(avg_lfp_psd(ep, hs).freq)]))));
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
        legend({lfp_tfa_cfg.epochs.name});
    end
    
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    saveas(gca, results_file);

end

