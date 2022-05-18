function lfp_tfa_plot_hs_tuned_tfr_band_average_combined( avg_tfr, lfp_tfa_cfg, plottitle, results_file, trial_type,  varargin )
%lfp_tfa_plot_hs_tuned_tfr_multiple_img  - Plots the LFP time frequency spectrogram
%averages for different hand-space conditions to be compared
%
% USAGE:
%   lfp_tfa_plot_hs_tuned_tfr_multiple_img( avg_tfr, lfp_tfa_cfg, plottitle, results_file )
%   lfp_tfa_plot_hs_tuned_tfr_multiple_img( avg_tfr, lfp_tfa_cfg, plottitle, results_file, cm )
%   lfp_tfa_plot_hs_tuned_tfr_multiple_img( avg_tfr, lfp_tfa_cfg, plottitle, results_file, cm, plot_significant )
%
%
% INPUTS:
%       avg_tfr         - average LFP time frequency response for different
%       hand-space conditions to be compared
%		lfp_tfa_cfg     - struct containing the required settings
%           Required Fields: see settings/lfp_tfa_settings_example
%               1. baseline_method             - method used for baseline
%               normalization
%               2. compare.reach_hands          - hand labels to compare
%               3. compare.reach_spaces         - space labels to compare
%       plottitle       - title for the plot
%       results_file    - path to filename to store the resulting image
%       varargin        - colormap to be used (default = 'jet', can be any
%                       standard colormap additionally supported is 'bluewhitered')
%                       - flag to indicate if only significant difference
%                       bins should be plotted
%
% REQUIRES:	bluewhitered, export_fig
%
% See also settings/lfp_tfa_settings_example, lfp_tfa_plot_site_average_tfr,
% lfp_tfa_avg_tfr_across_sessions, lfp_tfa_avg_tfr_across_sites,
% bluewhitered, colormap, lfp_tfa_compute_difference_condition_tfr
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
for nband = 1:length(fieldnames(lfp_tfa_cfg.band))
    band_labels = fieldnames(lfp_tfa_cfg.band);
    results_file_tmp = [];
    
    h = figure;
    set(h, 'position', [100, 100,900, 675]);
    hold on
    
    % plot_significant = 0;
    % if nargin > 5
    %     plot_significant = varargin{2};
    % end
    
    % colorbar title
    if strcmp(lfp_tfa_cfg.baseline_method, 'zscore')
        cbtitle = 'Z-score';
        imscale = [-1, 1];
    elseif strcmp(lfp_tfa_cfg.baseline_method, 'division')
        cbtitle = 'P / \mu';
        imscale = [0, 2];
    elseif strcmp(lfp_tfa_cfg.baseline_method, 'subtraction')
        cbtitle = 'P - \mu';
        imscale = [0 1e-8];
    elseif strcmp(lfp_tfa_cfg.baseline_method, 'relchange')
        cbtitle = '(P - \mu) / \mu';
        imscale = [-1, 1];
    end
    
    % offset samples
    noffset = 100;
    
    % number of subplots required
    nhandlabels = length(lfp_tfa_cfg.compare.reach_hands);
    nspacelabels = length(lfp_tfa_cfg.compare.reach_spaces);
    if trial_type == 1
        cond_to_plot = [1,3];
    elseif trial_type == 2
        cond_to_plot = [2,4];
    end
    concat_states_lfp = struct();
    for cnd = 1:length(cond_to_plot)
        
        % loop through handspace
        for hs = 1:size(avg_tfr(cond_to_plot(cnd)).hs_tuned_band, 2)
            
            
            if isfield(cat(2, avg_tfr(cond_to_plot(cnd)).hs_tuned_band.freq), (band_labels{nband}))
                freq =  cat(2, avg_tfr(cond_to_plot(cnd)).hs_tuned_band.freq);
            else
                continue;
            end
            
            % concatenate tfs for different state windows for plotting
            concat_states_tfs = struct();
            concat_states_tfs(cnd).powspctrm = [];
            concat_states_tfs(cnd).state_time = [];
            concat_states_tfs(cnd).freq = avg_tfr(cond_to_plot(cnd)).hs_tuned_band(cnd,hs).freq.freq;
            concat_states_tfs(cnd).label = avg_tfr(cond_to_plot(cnd)).hs_tuned_band(cnd,hs).hs_label;
            
            subplot(nhandlabels, nspacelabels, hs);
            hold on;
            
            state_info = struct();
            for st = 1:size(avg_tfr(cond_to_plot(cnd)).hs_tuned_band, 1)
                
                % state timing information
                % state onset sample number
                state_info(st).onset_s = find(...
                    avg_tfr(cond_to_plot(cnd)).hs_tuned_band(st,hs).freq.time <= 0, 1, 'last');
                % state onset time
                state_info(st).onset_t = 0;
                % start start sample
                state_info(st).start_s = 1;
                % state start time
                state_info(st).start_t = avg_tfr(cond_to_plot(cnd)).hs_tuned_band(st,hs).freq.time(1);
                % start finish sample
                state_info(st).finish_s = length(avg_tfr(cond_to_plot(cnd)).hs_tuned_band(st,hs).freq.time);
                % start end sample
                state_info(st).finish_t =avg_tfr(cond_to_plot(cnd)).hs_tuned_band(st,hs).freq.time(end);
                
                % state onset, start and finish samples for further states
                % offset from previous state window
                if st > 1
                    state_info(st).start_s = length(concat_states_tfs(cnd).state_time) + ...
                        state_info(st).start_s;
                    state_info(st).finish_s = length(concat_states_tfs(cnd).state_time) + ...
                        state_info(st).finish_s;
                    state_info(st).onset_s = length(concat_states_tfs(cnd).state_time) + ...
                        state_info(st).onset_s;
                end
                
                concat_states_tfs(cnd).powspctrm = cat(3, concat_states_tfs(cnd).powspctrm, ...
                    avg_tfr(cond_to_plot(cnd)).hs_tuned_band(st,hs).freq.(band_labels{nband}), ...
                    nan(size(avg_tfr(cond_to_plot(cnd)).hs_tuned_band(st,hs).freq.(band_labels{nband}), 1), 1, 100/25));
                concat_states_tfs(cnd).state_time = [concat_states_tfs(cnd).state_time, ...
                    avg_tfr(cond_to_plot(cnd)).hs_tuned_band(st,hs).freq.time, nan(1, 100/25)];
                
                
                state_lfp_powspctrm = nanmean(avg_tfr(cond_to_plot(cnd)).hs_tuned_band(st,hs).freq.(band_labels{nband}), 1);
                
                stderr = nanstd(avg_tfr(cond_to_plot(cnd)).hs_tuned_band(st,hs).freq.(band_labels{nband}), 1, 1)...
                    / sqrt(size(avg_tfr(cond_to_plot(cnd)).hs_tuned_band(st,hs).freq.(band_labels{nband}),1));
                state_lfp_powspctrm_error = [state_lfp_powspctrm + stderr...
                    state_lfp_powspctrm - stderr];
                if cnd == 1
                    plot(linspace(state_info(st).start_s, state_info(st).finish_s, ...
                        state_info(st).finish_s - state_info(st).start_s + 1),squeeze(state_lfp_powspctrm),'b')
                    plot(linspace(state_info(st).start_s, state_info(st).finish_s, ...
                        state_info(st).finish_s - state_info(st).start_s + 1),squeeze(state_lfp_powspctrm_error),'b--')
                    
                else
                    plot(linspace(state_info(st).start_s, state_info(st).finish_s, ...
                        state_info(st).finish_s - state_info(st).start_s + 1),squeeze(state_lfp_powspctrm),'r')
                    plot(linspace(state_info(st).start_s, state_info(st).finish_s, ...
                        state_info(st).finish_s - state_info(st).start_s + 1),squeeze(state_lfp_powspctrm_error),'r--')
                end
            end
            concat_states_tfs(cnd).time = 1:1:size(concat_states_tfs(cnd).powspctrm, 3);
            state_onsets = find(concat_states_tfs(cnd).state_time == 0);
            state_samples = sort([state_info.start_s, state_info.onset_s, ...
                state_info.finish_s]);
            
            % now plot
            
            %             set(get(cb,'title'),'string', cbtitle, 'fontsize',8);
            %             set(gca,'TickDir','out')
            
            % mark state onsets
            if cnd == 1
                
                
                
                set(gca,'xtick',state_samples)
                for so = state_onsets
                    line([so so], ylim, 'color', 'k');
                    if isfield(avg_tfr(cond_to_plot(cnd)).hs_tuned_band(state_onsets == so, hs), 'state_name') && ...
                            ~isempty(avg_tfr(cond_to_plot(cnd)).hs_tuned_band(state_onsets == so, hs).state_name)
                        state_name = avg_tfr(cond_to_plot(cnd)).hs_tuned_band(state_onsets == so, hs).state_name;
                        text(so+1, 10, state_name, 'fontsize', 8);
                    end
                end
                set(gca,'xticklabels', round(concat_states_tfs(cnd).state_time(state_samples), 1), 'fontsize', 8)
                set(gca, 'xticklabelrotation', 45)
                % add 0.5 since the time value is the center of the bin
                % add 0 at the beginning to make the y-axis visible
                set(gca, 'xlim', [0 state_samples(end) + 0.5]);
                xlabel('Time (s)');
                ylabel(cbtitle);
                
                subplottitle = concat_states_tfs(cnd).label{1};
                if isfield(avg_tfr(1, hs), 'nsessions')
                subplottitle = [subplottitle ' (nsessions = ' num2str(avg_tfr(1, hs).nsessions) ')'];
            elseif isfield(avg_tfr(cnd).hs_tuned_band(st,hs), 'nsites')
                subplottitle = [subplottitle ' (nsites = ' num2str(avg_tfr(cnd).hs_tuned_band(st,hs).nsites) ')'];
            elseif isfield(avg_tfr(1, hs), 'ntrials') && ~isempty(avg_tfr(1, hs).ntrials)
                subplottitle = [subplottitle ' (ntrials = ' num2str(avg_tfr(1, hs).ntrials) ')'];
            end
                
%                  if isfield(avg_tfr.hs_tuned_band(cnd,hs), 'nsessions')
%                     subplottitle = [subplottitle ' (nsessions = ' num2str(avg_tfr.hs_tuned_band(cnd,hs).nsessions) ')'];
%                 elseif isfield(avg_tfr.hs_tuned_band(cnd,hs), 'nsites')
%                     subplottitle = [subplottitle ' (nsites = ' num2str(avg_tfr.hs_tuned_band(cnd,hs).nsites) ')'];
%                 elseif isfield(avg_tfr.hs_tuned_band(cnd,hs), 'ntrials') && ~isempty(avg_tfr.hs_tuned_band(cnd,hs).ntrials)
%                     subplottitle = [subplottitle ' (ntrials = ' num2str(avg_tfr.hs_tuned_band(cnd,hs).ntrials) ')'];
%                 end
                title(subplottitle);
                %line([0 0], ylim, 'color', 'k');
                
                %change aspect ration if only 2 conditions
                %             if size(avg_tfr, 2) < 3
                %                 set(gca,'DataAspectRatio', [1 0.6 1]);
                %             end
                
            end
        end
    end
    % plot title
    plottitle = [plottitle ' ' band_labels{nband}];
    plottitle = strrep(plottitle, '_', '\_');
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', plottitle...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    
    
    
    
    
    %export_fig(h, results_file);
    results_file_tmp = [results_file ' ' band_labels{nband}];
    fig_formats = {'png'}; %default
    if isfield(lfp_tfa_cfg, 'save_fig_format') && ~isempty(lfp_tfa_cfg.save_fig_format)
        fig_formats = lfp_tfa_cfg.save_fig_format;
    end
    for fmt = fig_formats
        export_fig(h,  results_file_tmp, ['-' fmt{:}]);
    end
end


end


