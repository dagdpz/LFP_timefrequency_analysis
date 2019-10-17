function lfp_tfa_plot_Rpeak_ref_state_onsets( Rpeak_evoked_states, lfp_tfa_cfg, plottitle, results_file, varargin )
%lfp_tfa_plot_evoked_lfp  - Plots the LFP evoked response
%averages for different hand-space conditions to be compared
%
% USAGE:
%   lfp_tfa_plot_evoked_lfp( evoked_lfp, lfp_tfa_cfg, plottitle, results_file )
%
% INPUTS:
%       evoked_lfp       - average LFP power spectrum for different
%       hand-space conditions to be compared
%		lfp_tfa_cfg      - struct containing the required settings
%           Required Fields: see lfp_tfa_settings
%               1. 
%       plottitle        - title for the plot
%       results_file     - path to filename to store the resulting image
%
% REQUIRES:	
%
% See also lfp_tfa_settings, lfp_tfa_plot_site_evoked_LFP, 
% lfp_tfa_avg_evoked_LFP_across_sites, 
% lfp_tfa_avg_evoked_LFP_across_sessions
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
    
    % number of subplots required
    nstates = size(lfp_tfa_cfg.analyse_Rpeak_states, 1);
    nhandlabels = length(lfp_tfa_cfg.compare.reach_hands);
    nspacelabels = length(lfp_tfa_cfg.compare.reach_spaces);
    nhandspacelabels = nhandlabels * nspacelabels;
    
    % loop through handspace
    for hs = 1:size(Rpeak_evoked_states, 2)
        if ~isempty([Rpeak_evoked_states(:,hs).abs_timefromRpeak])
            for st = 1:size(Rpeak_evoked_states, 1)
                
%               % state onset probability vs. Rpeak phase
                timebinedges = Rpeak_evoked_states(st, hs).rel_timeprob.timebins * 2*pi;
                timebinmiddle = (timebinedges(1:end-1) + timebinedges(2:end))/2;
                mean_prob = nanmean(Rpeak_evoked_states(st, hs).rel_timeprob.prob, 1);
                null_hyp_prob = ones(size(mean_prob)) * (1/length(mean_prob));
                std_prob = nanstd(Rpeak_evoked_states(st, hs).rel_timeprob.prob, 1);
                
                % state onset probability vs relative time from Rpeak
                subplot(nhandspacelabels*2, nstates, (hs-1)*nstates*2 + st);
                % join the mean
                prob_plot = polar([timebinmiddle timebinmiddle(1)], ...
                    [mean_prob mean_prob(1)], 'r-');
                prob_plot.LineWidth = 2;
                % SN(14.10.2019) - not showing polar plot when hold on is written before
                % first polar plot
                hold on;
                polar([timebinmiddle timebinmiddle(1)], ...
                    [null_hyp_prob null_hyp_prob(1)], 'k-');
                view([90 -90]);
                %phase_plot.LineSize = 2;
                % standard deviation
                if size(Rpeak_evoked_states(st, hs).rel_timeprob.prob, 1) > 1
                    polar([timebinmiddle timebinmiddle(1)], ...
                        [mean_prob mean_prob(1)] + ...
                        [std_prob std_prob(1)], 'r--');
                    polar([timebinmiddle timebinmiddle(1)], ...
                        [mean_prob mean_prob(1)] - ...
                        [std_prob std_prob(1)], 'r--');
                end
                %set(gca, 'XTick', timebinedges);
                subplottitle = Rpeak_evoked_states(st, hs).state_name;
                if isfield(Rpeak_evoked_states(st, hs), 'ntrials') && ...
                        ~isempty(Rpeak_evoked_states(st, hs).ntrials)
                    subplottitle = [subplottitle, sprintf(' (%g)', ...
                        Rpeak_evoked_states(st, hs).ntrials)];
                end
                title(subplottitle); %box on;
                %xlabel('Rel. Time from Rpeak')
                ylabel('Rpeak phase[°]');
                
                % state onset probability vs absolute time around Rpeak
                subplot(nhandspacelabels*2, nstates, (hs-1)*nstates*2 + nstates + st);
                
                timebinedges = Rpeak_evoked_states(st, hs).abs_timeprob.timebins;
                timebinmiddle = (timebinedges(1:end-1) + timebinedges(2:end))/2;
                if size(Rpeak_evoked_states(st, hs).abs_timeprob.prob, 1) > 1
                    plot(timebinmiddle, Rpeak_evoked_states(st, hs).abs_timeprob.prob, ...
                        'Color', [0.6, 0.6, 0.6], 'LineStyle', '-', 'LineWidth', 0.5);
                end
                hold on;
                % join the mean
                plot(timebinmiddle, ...
                    nanmean(Rpeak_evoked_states(st, hs).abs_timeprob.prob, 1), 'r*-', 'LineWidth', 2);
                % standard deviation
                if size(Rpeak_evoked_states(st, hs).abs_timeprob.prob, 1) > 1
                    plot(timebinmiddle, ...
                        nanmean(Rpeak_evoked_states(st, hs).abs_timeprob.prob, 1) + ...
                        nanstd(Rpeak_evoked_states(st, hs).abs_timeprob.prob, 1), 'r--');
                    plot(timebinmiddle, ...
                        nanmean(Rpeak_evoked_states(st, hs).abs_timeprob.prob, 1) - ...
                        nanstd(Rpeak_evoked_states(st, hs).abs_timeprob.prob, 1), 'r--');
                end
                set(gca, 'XLim', lfp_tfa_cfg.analyse_Rpeak_states{st, 3});
                box on;
                % line at Rpeak onset
                line([0 0], ylim, 'Color', 'k', 'LineStyle', '--');
                %set(gca, 'XTick', timebinedges);
                %grid on;
                subplottitle = Rpeak_evoked_states(st, hs).state_name;
                if isfield(Rpeak_evoked_states(st, hs), 'ntrials') && ...
                        ~isempty(Rpeak_evoked_states(st, hs).ntrials)
                    subplottitle = [subplottitle, sprintf(' (%g)', ...
                        Rpeak_evoked_states(st, hs).ntrials)];
                end
                title(subplottitle); box on;
                xlabel('Time from Rpeak (s)')
                ylabel('P(onset)');
                
            end
            
        end
    end
    if isfield(Rpeak_evoked_states, 'nsessions')
        plottitle = [plottitle, ' (nsessions = ' num2str(Rpeak_evoked_states(1, hs).nsessions) ')'];
    end
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    export_fig(h, results_file);

end

