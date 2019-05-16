function lfp_tfa_plot_hs_tuned_tfr_multiple_img( avg_tfr, lfp_tfa_cfg, plottitle, results_file, varargin )
%lfp_tfa_plot_hs_tuned_tfr  - Plots the LFP time frequency spectrogram 
%averages for different hand-space conditions to be compared
%
% USAGE:
%   lfp_tfa_plot_hs_tuned_tfr( avg_tfr, lfp_tfa_cfg, plottitle, results_file )
%   lfp_tfa_plot_hs_tuned_tfr( avg_tfr, lfp_tfa_cfg, plottitle, results_file, 'bluewhitered' )
%
%
% INPUTS:
%       avg_tfr         - average LFP time frequency response for different
%       hand-space conditions to be compared
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields: see lfp_tfa_settings
%               1. baseline_method             - method used for baseline
%               normalization
%               2. compare.reach_hands          - hand labels to compare
%               3. compare.reach_spaces         - space labels to compare
%       plottitle       - title for the plot
%       results_file    - path to filename to store the resulting image
%       varargin        - colormap to be used (default = 'jet', can be any 
%       standard colormap additionally supported is 'bluewhitered')
%
% REQUIRES:	bluewhitered
%
% See also lfp_tfa_settings, lfp_tfa_plot_site_average_tfr, 
% lfp_tfa_avg_tfr_across_sessions, lfp_tfa_avg_tfr_across_sites, 
% bluewhitered, lfp_tfa_define_settings
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
    
    figure;
    hold on
    
    % colorbar title
    if strcmp(lfp_tfa_cfg.baseline_method, 'zscore')
        cbtitle = 'Z-score';
    elseif strcmp(lfp_tfa_cfg.baseline_method, 'division')
        cbtitle = 'P / \mu';
    elseif strcmp(lfp_tfa_cfg.baseline_method, 'subraction')
        cbtitle = 'P - \mu';
    elseif strcmp(lfp_tfa_cfg.baseline_method, 'relchange')
        cbtitle = '(P - \mu) / \mu';
    end
    
    % offset samples
    noffset = 100;
    % number of subplots required
    nhandlabels = length(lfp_tfa_cfg.compare.reach_hands);
    nspacelabels = length(lfp_tfa_cfg.compare.reach_spaces);
    
    % loop through handspace
    for hs = 1:size(avg_tfr, 2)
        % check if no trials exist for this condition and HS
        if ~isempty(cat(3, avg_tfr(:, hs).powspctrm))
            % concatenate tfs for different state windows for plotting
            concat_states_tfs = struct();
            concat_states_tfs.powspctrm = [];
            concat_states_tfs.state_time = [];
            concat_states_tfs.freq = avg_tfr(1, hs).freq;
            concat_states_tfs.label = avg_tfr(1, hs).hs_label;
            
            subplot(nhandlabels, nspacelabels, hs);
            hold on;

            state_info = struct();
            for st = 1:size(avg_tfr, 1)

                concat_states_tfs.powspctrm = cat(3, concat_states_tfs.powspctrm, ...
                    avg_tfr(st, hs).powspctrm);
                concat_states_tfs.state_time = [concat_states_tfs.state_time, ...
                    avg_tfr(st, hs).time];
                % append nans to separate the states
                if st < size(avg_tfr, 1)
                    concat_states_tfs.powspctrm = cat(3, concat_states_tfs.powspctrm, ...
                        nan(1, length(concat_states_tfs.freq), 100/25));
                    concat_states_tfs.state_time = [concat_states_tfs.state_time, ...
                        nan(1, 100/25)];
                end
                
                % state timing information
                % state onset sample number
                state_info(st).onset_s = find(...
                    avg_tfr(st, hs).time <= 0, 1, 'last'); 
                % state onset time
                state_info(st).onset_t = 0; 
                % start start sample
                state_info(st).start_s = 1;
                % state start time
                state_info(st).start_t = avg_tfr(st, hs).time(1);
                % start finish sample
                state_info(st).finish_s = length(avg_tfr(st, hs).time);
                % start end sample
                state_info(st).finish_t = avg_tfr(st, hs).time(end);                    
                
                % state onset, start and finish samples for further states
                % offset from previous state window
                if st > 1
                    state_info(st).start_s = length(avg_tfr(st-1, hs).time) + ...
                        state_info(st).start_s + (st-1)*(100/25);
                    state_info(st).finish_s = length(avg_tfr(st-1, hs).time) + ...
                        state_info(st).finish_s + (st-1)*(100/25);
                    state_info(st).onset_s = length(avg_tfr(st-1, hs).time) + ...
                        state_info(st).onset_s + (st-1)*(100/25);
                end
                
                imagesc(linspace(state_info(st).start_s, state_info(st).finish_s, ...
                    state_info(st).finish_s - state_info(st).start_s + 1), ...
                    linspace(1, length(avg_tfr(st, hs).freq), length(avg_tfr(st, hs).freq)), ...
                    squeeze(avg_tfr(st, hs).powspctrm) ,[-1 1]);
                
                % horizontal lines to separate frequency bands
                for f = [2, 4, 8, 12, 18, 32, 80]
                    f_idx = find(abs(avg_tfr(st, hs).freq - f) == ...
                        min(abs(avg_tfr(st, hs).freq - f)), 1, 'first');
                    line([state_info(st).start_s state_info(st).finish_s], ...
                        [f_idx f_idx], 'color', 'k', 'linestyle', '--');
                end


            end
            concat_states_tfs.time = 1:1:size(concat_states_tfs.powspctrm, 3);
            state_onsets = find(concat_states_tfs.state_time(1:end-1) .* ...
                concat_states_tfs.state_time(2:end) <= 0);
            state_samples = sort([state_info.start_s, state_info.onset_s, ...
                state_info.finish_s]);

            % now plot
            
            %subplot(nhandlabels, nspacelabels, hs)
            %imagesc(concat_states_tfs.time, [1:numel(concat_states_tfs.freq)], squeeze(concat_states_tfs.powspctrm), [-1 1]);
            axis xy, cb = colorbar;
            set(get(cb,'title'),'string', cbtitle, 'fontsize',6);
            set(gca,'TickDir','out')
            % log y axis ticks
            set(gca, 'ytick', ([1:8:numel(concat_states_tfs.freq)]));
            set(gca, 'yticklabel', ...
                round(concat_states_tfs.freq([1:8:numel(concat_states_tfs.freq)])));
            % add 0.5 at end since the time value is the center of the bin
            % add 0 at beginning to make x-axis visible
            set(gca, 'ylim', [0 numel(avg_tfr(st, hs).freq) + 0.5]);
            % mark state onsets
            set(gca,'xtick',state_samples)
            for so = state_onsets
                line([so so], ylim, 'color', 'k'); 
                state_name = avg_tfr(state_onsets == so, hs).state_name;
                text(so+1, 10, state_name, 'fontsize', 8);
            end
            set(gca,'xticklabels', round(concat_states_tfs.state_time(state_samples), 1), 'fontsize', 7)
            set(gca, 'xticklabelrotation', 45)
            % add 0.5 since the time value is the center of the bin
            % add 0 at the beginning to make the y-axis visible
            set(gca, 'xlim', [0 state_samples(end) + 0.5]);
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
                        
            subplottitle = concat_states_tfs.label{1};
            if isfield(avg_tfr(1, hs), 'nsessions')
                subplottitle = [subplottitle ' (nsessions = ' num2str(avg_tfr(1, hs).nsessions) ')'];
            end
            if isfield(avg_tfr(1, hs), 'nsites')
                subplottitle = [subplottitle ' (nsites = ' num2str(avg_tfr(1, hs).nsites) ')'];
            end
            if isfield(avg_tfr(1, hs), 'ntrials') && ~isempty(avg_tfr(1, hs).ntrials)
                subplottitle = [subplottitle ' (ntrials = ' num2str(avg_tfr(1, hs).ntrials) ')'];            
            end
            title(subplottitle);
            %line([0 0], ylim, 'color', 'k');
            
        end
    end
    
    % plot title
    plottitle = strrep(plottitle, '_', '\_');
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', plottitle...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    
    % define colormap
    cm = colormap('jet'); % default
    if nargin > 4
        cm = colormap(varargin{1});
        colorbar;
    end
    % white separation between two state windows - commented since there is
    % a new method for separating windows - test and remove
    %cm(1,:,:) = [1,1,1];
    colormap(cm);
    
    saveas(gca, results_file);

end

