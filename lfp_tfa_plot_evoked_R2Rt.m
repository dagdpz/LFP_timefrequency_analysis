function lfp_tfa_plot_evoked_R2Rt( evoked_R2Rt, lfp_tfa_cfg, plottitle, results_file, varargin )
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


    % defaults
    if lfp_tfa_cfg.normalize_R2Rt
        yaxislabel = 'Relative R2R time';
    else
        yaxislabel = 'R2R time (s)';
    end
    err = 'stdev'; % standard error
    
    % get settings
    if nargin > 4
        settings = struct(varargin{:});
        if isfield(settings, 'ylabel')
            yaxislabel = settings.ylabel;
        end
        if isfield(settings, 'err')
            err = settings.err;
        end
        
    end
    
    h1 = figure;
    
    %set(h1, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    set(h1, 'position', [100, 100,900, 675]);
    
    % number of offset samples to divide between time windows
    noffset = 150;
    
    % number of subplots required
    nhandlabels = length(lfp_tfa_cfg.compare.reach_hands);
    nspacelabels = length(lfp_tfa_cfg.compare.reach_spaces);
    
    % loop through handspace
    for hs = 1:size(evoked_R2Rt, 2)
        if ~isempty([evoked_R2Rt(:,hs).mean]) &&  ~isempty([evoked_R2Rt(:,hs).std])
            % concatenate states
            concat_states_R2Rt = struct();
            concat_states_R2Rt.trial = {};
            concat_states_R2Rt.mean = [];
            concat_states_R2Rt.err = [];
            concat_states_R2Rt.time = [];

            state_info = struct();
            for st = 1:size(evoked_R2Rt, 1)
                
                state_info(st).onset_s = find(...
                    evoked_R2Rt(st, hs).time <= 0, 1, 'last');
                state_info(st).onset_t = 0;
                state_info(st).start_s = 1;
                state_info(st).start_t = evoked_R2Rt(st, hs).time(1);
                state_info(st).finish_s = length(evoked_R2Rt(st, hs).time);
                state_info(st).finish_t = evoked_R2Rt(st, hs).time(end);                    

                if st > 1
                    state_info(st).start_s = length(concat_states_R2Rt.time) + ...
                        state_info(st).start_s;
                    state_info(st).finish_s = length(concat_states_R2Rt.time) + ...
                        state_info(st).finish_s;
                    state_info(st).onset_s = length(concat_states_R2Rt.time) + ...
                        state_info(st).onset_s;
                end

                % concatenate mean, std and time of evoked LFP for
                % different states
                if isfield(evoked_R2Rt(1, hs), 'ntrials') && ...
                     ~isempty(evoked_R2Rt(1, hs).ntrials)
                    state_evoked_R2Rt = evoked_R2Rt(st, hs).ecg_b2bt;
                    if lfp_tfa_cfg.normalize_R2Rt
                           state_evoked_R2Rt = state_evoked_R2Rt ./ ...
                               repmat(evoked_R2Rt(st, hs).trials_mean', ...
                               [1 size(state_evoked_R2Rt, 2)]);
                    end
                    concat_states_R2Rt.trial = [concat_states_R2Rt.trial, ...
                    horzcat(state_evoked_R2Rt, nan(size(evoked_R2Rt(st, hs).ecg_b2bt, 1), noffset))];
                end
                
                state_evoked_mean = evoked_R2Rt(st, hs).mean;
                state_evoked_err = evoked_R2Rt(st, hs).std;
                if strcmp(err, 'stderr')
                    state_evoked_err = evoked_R2Rt(st, hs).std ...
                        / sqrt(size(evoked_R2Rt(st, hs).ecg_b2bt, 1));
                elseif strcmp(err, 'stdev')
                    state_evoked_err = evoked_R2Rt(st, hs).std;
                end                
                
                concat_states_R2Rt.mean = [concat_states_R2Rt.mean, ...
                    state_evoked_mean, nan(size(evoked_R2Rt(st, hs).mean, 1), noffset)];
                concat_states_R2Rt.err = [concat_states_R2Rt.err, ...
                    state_evoked_err, nan(size(evoked_R2Rt(st, hs).std, 1), noffset)];
                concat_states_R2Rt.time = [concat_states_R2Rt.time, ...
                    evoked_R2Rt(st, hs).time, nan(1, noffset)];
                concat_states_R2Rt.label = evoked_R2Rt(st, hs).hs_label;
                if isfield(evoked_R2Rt(st, hs), 'legend')
                    concat_states_R2Rt.legend = evoked_R2Rt(st, hs).legend;
                end

            end
            
            %state_onsets = find(concat_states_lfp.time(1:end-1) .* ...
            %    concat_states_lfp.time(2:end) <= 0);
            state_onsets = find(concat_states_R2Rt.time == 0);
            state_samples = sort([state_info.start_s, state_info.onset_s, ...
                state_info.finish_s]);

            % now plot
            subplot(nhandlabels*2, nspacelabels, hs)
            hold on;
            colors = ['b'; 'r'; 'g'; 'y'; 'm'; 'c'; 'k'];
            if isfield(evoked_R2Rt(1, hs), 'color') && ...
                    ~isempty(evoked_R2Rt(1, hs).color)
                colors = evoked_R2Rt(1, hs).color;
            end
            % plot individual trials
            if isfield(evoked_R2Rt(1, hs), 'ntrials') && ~isempty(evoked_R2Rt(1, hs).ntrials)
                xx = [0];
                for s = 1:length(concat_states_R2Rt.trial)
                    xx = linspace(1, size(concat_states_R2Rt.trial{s}, 2), ...
                        size(concat_states_R2Rt.trial{s}, 2)) + xx(end);
                    plot(xx, concat_states_R2Rt.trial{s}, 'Color', [0.6, 0.6, 0.6])
                end
            end
            for i = 1:size(concat_states_R2Rt.mean, 1)
                plot(concat_states_R2Rt.mean(i, :), 'Color', colors(i, :), 'LineWidth', 2);
            end
            if isfield(concat_states_R2Rt, 'legend')
                legend(concat_states_R2Rt.legend);
            end
            for i = 1:size(concat_states_R2Rt.mean, 1)
                plot(concat_states_R2Rt.mean(i, :) + concat_states_R2Rt.err(i, :), ...
                    '--', 'Color', colors(i, :), 'LineWidth', 1);
                plot(concat_states_R2Rt.mean(i, :) - concat_states_R2Rt.err(i, :), ...
                    '--', 'Color', colors(i, :), 'LineWidth', 1);
            end
            % mark state onsets
            %if isfield(evoked_lfp, 'state_name')
            set(gca,'xtick', unique(state_samples))
            set(gca, 'YLim', ylim);
            for so = state_onsets
                line([so so], ylim, 'color', 'k'); 
                if isfield(evoked_R2Rt(state_onsets == so, hs), 'state_name') && ...
                        ~isempty(evoked_R2Rt(state_onsets == so, hs).state_name)
                    state_name = evoked_R2Rt(state_onsets == so, hs).state_name;
                    plottxt = state_name;
                    if isfield(evoked_R2Rt(state_onsets == so, hs), 'ntrials') && ...
                        ~isempty(evoked_R2Rt(state_onsets == so, hs).ntrials)
                        ntrials = (evoked_R2Rt(state_onsets == so, hs).ntrials);
                        plottxt = sprintf('%s \n(%g)', plottxt, ntrials);
                    end
                    ypos = ylim;
                    ypos = ypos(1) + (ypos(2) - ypos(1))*0.2;
                    text(so+1, ypos, plottxt, 'fontsize', 10);
                end
            end
            %end
            set(gca,'xticklabels', round(concat_states_R2Rt.time(unique(state_samples)), 2), 'fontsize', 10)
            set(gca, 'xticklabelrotation', 45);
            xlabel('Time(s)');
            ylabel(yaxislabel);
            
            subplottitle = [concat_states_R2Rt.label{1}];
            if isfield(evoked_R2Rt(1, hs), 'nsessions')
                subplottitle = [subplottitle ' (nsessions = ' num2str(evoked_R2Rt(1, hs).nsessions) ')'];
%             elseif isfield(evoked_R2Rt(1, hs), 'ntrials') && ...
%                     ~isempty(evoked_R2Rt(1, hs).ntrials)
%                 subplottitle = [subplottitle ' (ntrials = ' ...
%                     num2str(evoked_R2Rt(1, hs).ntrials) ')'];            
            end
            title(subplottitle);
        end
    end
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    export_fig(h1, [results_file '.png']);
    print(h1, '-depsc', [results_file '.ai']);

end
