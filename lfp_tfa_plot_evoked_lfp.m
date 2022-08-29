function lfp_tfa_plot_evoked_lfp( evoked_lfp, lfp_tfa_cfg, plottitle, results_file )
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
%           Required Fields: see settings/lfp_tfa_settings_example
%               1. compare.reach_hands  - reach hands to be included
%               2. compare.reach_spaces - reach spaces to be included
%       plottitle        - title for the plot
%       results_file     - path to filename to store the resulting image
%
% REQUIRES:	export_fig
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_plot_site_evoked_LFP, 
% lfp_tfa_avg_evoked_LFP_across_sites,
% lfp_tfa_avg_evoked_LFP_across_sessions, external/export_fig/export_fig.m
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

    % number of offset samples to divide between time windows
    noffset = 100;
    
    % number of subplots required
    nhandlabels = length(lfp_tfa_cfg.compare.reach_hands);
    nspacelabels = length(lfp_tfa_cfg.compare.reach_spaces);
    
    % loop through handspace
    for hs = 1:size(evoked_lfp, 2)
        if ~isempty([evoked_lfp(:,hs).lfp])% &&  ~isempty([evoked_lfp(:,hs).std])
            % concatenate states
            concat_states_lfp = struct();
            concat_states_lfp.mean = [];
            concat_states_lfp.error = [];
            concat_states_lfp.time = [];

            state_info = struct();
            for st = 1:size(evoked_lfp, 1)
                
                state_info(st).onset_s = find(...
                    evoked_lfp(st, hs).time <= 0, 1, 'last');
                state_info(st).onset_t = 0;
                state_info(st).start_s = 1;
                state_info(st).start_t = evoked_lfp(st, hs).time(1);
                state_info(st).finish_s = length(evoked_lfp(st, hs).time);
                state_info(st).finish_t = evoked_lfp(st, hs).time(end);                    

                if st > 1
                    state_info(st).start_s = length(concat_states_lfp.time) + ...
                        state_info(st).start_s;
                    state_info(st).finish_s = length(concat_states_lfp.time) + ...
                        state_info(st).finish_s;
                    state_info(st).onset_s = length(concat_states_lfp.time) + ...
                        state_info(st).onset_s;
                end

                % concatenate mean, std and time of evoked LFP for
                % different states
                [state_mean, state_error] = ...
                    lfp_tfa_compute_statistics(evoked_lfp(st, hs).lfp, lfp_tfa_cfg.error_measure);
                concat_states_lfp.mean = [concat_states_lfp.mean, ...
                    state_mean, nan(1, noffset)];
                concat_states_lfp.error = [concat_states_lfp.error, ...
                    state_error, nan(2, noffset)];
                concat_states_lfp.time = [concat_states_lfp.time, ...
                    evoked_lfp(st, hs).time, nan(1, noffset)];
                concat_states_lfp.label = evoked_lfp(st, hs).hs_label;         

                

            end
            
            %state_onsets = find(concat_states_lfp.time(1:end-1) .* ...
            %    concat_states_lfp.time(2:end) <= 0);
            state_onsets = find(concat_states_lfp.time == 0);
            state_samples = sort([state_info.start_s, state_info.onset_s, ...
                state_info.finish_s]);

            % now plot
            subplot(nhandlabels, nspacelabels, hs)
            hold on;
            plot(concat_states_lfp.mean, 'b');
%             plot(concat_states_lfp.error', 'b--');
            %plot(concat_states_lfp.mean - concat_states_lfp.std, 'r--');
            % mark state onsets
            %if isfield(evoked_lfp, 'state_name')
            set(gca,'xtick',state_samples)
            for so = state_onsets
                line([so so], ylim, 'color', 'k'); 
                if isfield(evoked_lfp(state_onsets == so, hs), 'state_name') && ...
                        ~isempty(evoked_lfp(state_onsets == so, hs).state_name)
                    state_name = evoked_lfp(state_onsets == so, hs).state_name;
                    ypos = ylim;
                    ypos = ypos(1) + (ypos(2) - ypos(1))*0.2;
                    text(so+1, ypos, state_name, 'fontsize', 8);
                end
            end
            %end
            set(gca,'xticklabels', round(concat_states_lfp.time(state_samples), 2), 'fontsize', 8)
            set(gca, 'xticklabelrotation', 45);
            xlabel('Time(s)');
            ylabel('LFP amplitude');
            subplottitle = [concat_states_lfp.label{1}];
            if isfield(evoked_lfp(1, hs), 'nsessions')
                subplottitle = [subplottitle ' (nsessions = ' num2str(evoked_lfp(1, hs).nsessions) ')'];
            elseif isfield(evoked_lfp(1, hs), 'nsites')
                subplottitle = [subplottitle ' (nsites = ' num2str(evoked_lfp(1, hs).nsites) ')'];
            elseif isfield(evoked_lfp(1, hs), 'trials')
                subplottitle = [subplottitle ' (ntrials = ' ...
                    num2str(length(evoked_lfp(1, hs).trials)) ')'];            
            end
            title(subplottitle);
        end
    end
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    
    fig_formats = {'png'}; %default
    if isfield(lfp_tfa_cfg, 'save_fig_format') && ~isempty(lfp_tfa_cfg.save_fig_format)
        fig_formats = lfp_tfa_cfg.save_fig_format;
    end
    for fmt = fig_formats
        export_fig(h, results_file, ['-' fmt{:}]);
    end

end

