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
    
    yaxislabel = 'LFP amplitude';
    if nargin > 4
        settings = struct(varargin{:});
        yaxislabel = settings.ylabel;
    end
    

    % number of offset samples to divide between time windows
    noffset = 100;
    
    % number of subplots required
    nstates = size(lfp_tfa_cfg.analyse_Rpeak_states, 1);
    nhandlabels = length(lfp_tfa_cfg.compare.reach_hands);
    nspacelabels = length(lfp_tfa_cfg.compare.reach_spaces);
    nhandspacelabels = nhandlabels * nspacelabels;
    
    % loop through handspace
    for hs = 1:size(Rpeak_evoked_states, 2)
        if ~isempty([Rpeak_evoked_states(:,hs).abs_timefromRpeak])
%             % concatenate states
%             concat_states_lfp = struct();
%             concat_states_lfp.mean = [];
%             concat_states_lfp.std = [];
%             concat_states_lfp.time = [];
% 
%             state_info = struct();
            for st = 1:size(Rpeak_evoked_states, 1)
                
%                 state_info(st).onset_s = find(...
%                     Rpeak_evoked_state_onsets(st, hs).time <= 0, 1, 'last');
%                 state_info(st).onset_t = 0;
%                 state_info(st).start_s = 1;
%                 state_info(st).start_t = Rpeak_evoked_state_onsets(st, hs).time(1);
%                 state_info(st).finish_s = length(Rpeak_evoked_state_onsets(st, hs).time);
%                 state_info(st).finish_t = Rpeak_evoked_state_onsets(st, hs).time(end);                    
% 
%                 if st > 1
%                     state_info(st).start_s = length(concat_states_lfp.time) + ...
%                         state_info(st).start_s;
%                     state_info(st).finish_s = length(concat_states_lfp.time) + ...
%                         state_info(st).finish_s;
%                     state_info(st).onset_s = length(concat_states_lfp.time) + ...
%                         state_info(st).onset_s;
%                 end
% 
%                 % concatenate mean, std and time of evoked LFP for
%                 % different states
%                 concat_states_lfp.mean = [concat_states_lfp.mean, ...
%                     Rpeak_evoked_state_onsets(st, hs).mean, nan(size(Rpeak_evoked_state_onsets(st, hs).mean, 1), noffset)];
%                 concat_states_lfp.std = [concat_states_lfp.std, ...
%                     Rpeak_evoked_state_onsets(st, hs).std, nan(size(Rpeak_evoked_state_onsets(st, hs).std, 1), noffset)];
%                 concat_states_lfp.time = [concat_states_lfp.time, ...
%                     Rpeak_evoked_state_onsets(st, hs).time, nan(1, noffset)];
%                 concat_states_lfp.label = Rpeak_evoked_state_onsets(st, hs).hs_label;
%                 if isfield(Rpeak_evoked_state_onsets(st, hs), 'legend')
%                     concat_states_lfp.legend = Rpeak_evoked_state_onsets(st, hs).legend;
%                 end
%                 
%                 % now plot
                subplot(nhandspacelabels*2, nstates, (hs-1)*nstates*2 + st);
                
%                 histogram(Rpeak_evoked_states(st, hs).abs_timefromRpeak, 'BinWidth', 0.05, ...
%                     'Normalization', 'probability')
                timebinedges = Rpeak_evoked_states(st, hs).abs_timeprob.timebins;
                timebinmiddle = (timebinedges(1:end-1) + timebinedges(2:end))/2;
                plot(timebinmiddle, Rpeak_evoked_states(st, hs).abs_timeprob.prob, 'bs');
                hold on;
                % join the mean
                plot(timebinmiddle, ...
                    nanmean(Rpeak_evoked_states(st, hs).abs_timeprob.prob, 1), 'r*-');
                %set(gca, 'XTick', timebinedges);
                %grid on;
                subplottitle = Rpeak_evoked_states(st, hs).state_name;
                if isfield(Rpeak_evoked_states(st, hs), 'ntrials') && ...
                        ~isempty(Rpeak_evoked_states(st, hs).ntrials)
                    subplottitle = [subplottitle, sprintf(' (%g)', ...
                        Rpeak_evoked_states(st, hs).ntrials)];
                end
                title(subplottitle);
                xlabel('Time from Rpeak (s)')
                ylabel('P(onset)');
                
                subplot(nhandspacelabels*2, nstates, (hs-1)*nstates*2 + nstates + st);
                
%                 histogram(Rpeak_evoked_states(st, hs).rel_timefromRpeak, 'BinWidth', 0.1, ...
%                     'Normalization', 'probability')
                timebinedges = Rpeak_evoked_states(st, hs).rel_timeprob.timebins;
                timebinmiddle = (timebinedges(1:end-1) + timebinedges(2:end))/2;
                plot(timebinmiddle, Rpeak_evoked_states(st, hs).rel_timeprob.prob, 'bs');
                hold on;
                % join the mean
                plot(timebinmiddle, ...
                    nanmean(Rpeak_evoked_states(st, hs).rel_timeprob.prob, 1), 'r*-');
                %set(gca, 'XTick', timebinedges);
                subplottitle = Rpeak_evoked_states(st, hs).state_name;
                if isfield(Rpeak_evoked_states(st, hs), 'ntrials') && ...
                        ~isempty(Rpeak_evoked_states(st, hs).ntrials)
                    subplottitle = [subplottitle, sprintf(' (%g)', ...
                        Rpeak_evoked_states(st, hs).ntrials)];
                end
                title(subplottitle);
                xlabel('Rel. Time from Rpeak')
                ylabel('P(onset)');

            end
            
            %state_onsets = find(concat_states_lfp.time(1:end-1) .* ...
            %    concat_states_lfp.time(2:end) <= 0);
%             state_onsets = find(concat_states_lfp.time == 0);
%             state_samples = sort([state_info.start_s, state_info.onset_s, ...
%                 state_info.finish_s]);
% 
%             % now plot
%             subplot(nhandspacelabels, nstates, (hs-1)*nstates + st)
%             hold on;
%             colors = ['b', 'r', 'g', 'y', 'm', 'c', 'k'];
%             for i = 1:size(concat_states_lfp.mean, 1)
%                 plot(concat_states_lfp.mean(i, :), colors(i));
%             end
%             if isfield(concat_states_lfp, 'legend')
%                 legend(concat_states_lfp.legend);
%             end
%             for i = 1:size(concat_states_lfp.mean, 1)
%                 plot(concat_states_lfp.mean(i, :) + concat_states_lfp.std(i, :), [colors(mod(i, length(colors))) ':']);
%                 plot(concat_states_lfp.mean(i, :) - concat_states_lfp.std(i, :), [colors(mod(i, length(colors))) ':']);
%             end
%             % mark state onsets
%             %if isfield(evoked_lfp, 'state_name')
%             set(gca,'xtick',state_samples)
%             for so = state_onsets
%                 line([so so], ylim, 'color', 'k'); 
%                 if isfield(Rpeak_evoked(state_onsets == so, hs), 'state_name') && ...
%                         ~isempty(Rpeak_evoked(state_onsets == so, hs).state_name)
%                     state_name = Rpeak_evoked(state_onsets == so, hs).state_name;
%                     ypos = ylim;
%                     ypos = ypos(1) + (ypos(2) - ypos(1))*0.2;
%                     text(so+1, ypos, state_name, 'fontsize', 8);
%                 end
%             end
%             %end
%             set(gca,'xticklabels', round(concat_states_lfp.time(state_samples), 2), 'fontsize', 8)
%             set(gca, 'xticklabelrotation', 45);
%             xlabel('Time(s)');
%             ylabel(yaxislabel);
%             
%             subplottitle = [concat_states_lfp.label{1}];
%             if isfield(Rpeak_evoked(1, hs), 'nsessions')
%                 subplottitle = [subplottitle ' (nsessions = ' num2str(Rpeak_evoked(1, hs).nsessions) ')'];
%             elseif isfield(Rpeak_evoked(1, hs), 'nsites')
%                 subplottitle = [subplottitle ' (nsites = ' num2str(Rpeak_evoked(1, hs).nsites) ')'];
%             elseif isfield(Rpeak_evoked(1, hs), 'trials')
%                 subplottitle = [subplottitle ' (ntrials = ' ...
%                     num2str(length(Rpeak_evoked(1, hs).trials)) ')'];            
%             end
%             title(subplottitle);
        end
    end
    if isfield(Rpeak_evoked_states, 'nsessions')
        plottitle = [plottitle, ' (nsessions = ' num2str(Rpeak_evoked_states(1, hs).nsessions) ')'];
    end
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    export_fig(h, results_file);

end

