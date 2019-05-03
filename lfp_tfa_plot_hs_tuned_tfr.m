function lfp_tfa_plot_hs_tuned_tfr( avg_tfr, lfp_tfa_cfg, plottitle, results_file, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    figure;    
    
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
    % loop through handspace
    for hs = 1:size(avg_tfr, 2)
        % check if no trials exist for this condition and HS
        if ~isempty(cat(3, avg_tfr(:, hs).powspctrm))
            % concatenate states
            concat_states_tfs = struct();
            concat_states_tfs.powspctrm = [];
            concat_states_tfs.state_time = [];
            concat_states_tfs.freq = avg_tfr(1, hs).freq;
            concat_states_tfs.label = avg_tfr(1, hs).hs_label;

            state_info = struct();
            for st = 1:size(avg_tfr, 1)


                concat_states_tfs.powspctrm = cat(3, concat_states_tfs.powspctrm, ...
                    avg_tfr(st, hs).powspctrm);
                concat_states_tfs.state_time = [concat_states_tfs.state_time, ...
                    avg_tfr(st, hs).time];
                % append nans to separate the states
                if st < size(avg_tfr, 1)
                    concat_states_tfs.powspctrm = cat(3, concat_states_tfs.powspctrm, ...
                        nan(1, length(concat_states_tfs.freq), 300/25));
                    concat_states_tfs.state_time = [concat_states_tfs.state_time, ...
                        nan(1, 300/25)];
                end

                state_info(st).onset_s = find(...
                    avg_tfr(st, hs).time <= 0, 1, 'last');
                state_info(st).onset_t = 0;
                state_info(st).start_s = 1;
                state_info(st).start_t = avg_tfr(st, hs).time(1);
                state_info(st).finish_s = length(avg_tfr(st, hs).time);
                state_info(st).finish_t = avg_tfr(st, hs).time(end);                    

                if st > 1
                    state_info(st).start_s = length(avg_tfr(st-1, hs).time) + ...
                        state_info(st).start_s + (st-1)*(300/25);
                    state_info(st).finish_s = length(avg_tfr(st-1, hs).time) + ...
                        state_info(st).finish_s + (st-1)*(300/25);
                    state_info(st).onset_s = length(avg_tfr(st-1, hs).time) + ...
                        state_info(st).onset_s + (st-1)*(300/25);
                end


            end
            concat_states_tfs.time = 1:1:size(concat_states_tfs.powspctrm, 3);
            state_onsets = find(concat_states_tfs.state_time(1:end-1) .* ...
                concat_states_tfs.state_time(2:end) <= 0);
            state_samples = sort([state_info.start_s, state_info.onset_s, ...
                state_info.finish_s]);

            % now plot
            nhandlabels = length(lfp_tfa_cfg.compare.reach_hands);
            nspacelabels = length(lfp_tfa_cfg.compare.reach_spaces);
            subplot(nhandlabels, nspacelabels, hs)
            cfg = [];
%                 cfg.baseline     = 'no'; %baseline_shift;                 % -400ms to -100ms before the onset of first state
%                 cfg.maskstyle    = 'saturation';
%                 cfg.interactive  = 'no';
%                 cfg.zlim         = [-1 1];
            cfg.channel  = concat_states_tfs.label;
%                 subplot(2,2,hs)
            %ft_singleplotTFR(cfg, concat_states_tfs);
            %Z = zeros(size(squeeze(concat_states_tfs.powspctrm)));
            imagesc(concat_states_tfs.time, [1:numel(concat_states_tfs.freq)], squeeze(concat_states_tfs.powspctrm), [-1 1]);
            axis xy, cb = colorbar;
            set(get(cb,'title'),'string', cbtitle, 'fontsize',6);
            % log y axis ticks
            set(gca, 'ytick', ([1:8:numel(concat_states_tfs.freq)]));
            set(gca, 'yticklabel', ...
                round(concat_states_tfs.freq([1:8:numel(concat_states_tfs.freq)])));
            % mark state onsets
            set(gca,'xtick',state_samples)
%                 xticklabels = [];
            for so = state_onsets
                line([so so], ylim, 'color', 'k'); 
                state_name = avg_tfr(state_onsets == so, hs).state_name;
                text(so, 10, state_name);
            end
            set(gca,'xticklabels', round(concat_states_tfs.state_time(state_samples), 1), 'fontsize', 8)
            xlabel('Time (s)');
            ylabel('Frequency (Hz)');
            subplottitle = concat_states_tfs.label{1};
            if isfield(avg_tfr(1, hs), 'nsessions')
                subplottitle = [subplottitle ' (nsessions = ' num2str(avg_tfr(1, hs).nsessions) ')'];
            elseif isfield(avg_tfr(1, hs), 'nsites')
                subplottitle = [subplottitle ' (nsites = ' num2str(avg_tfr(1, hs).nsites) ')'];
            elseif isfield(avg_tfr(1, hs), 'ntrials') && ~isempty(avg_tfr(1, hs).ntrials)
                subplottitle = [subplottitle ' (ntrials = ' num2str(avg_tfr(1, hs).ntrials) ')'];            
            end
            title(subplottitle);
            line([0 0], ylim, 'color', 'k');
            % horizontal lines to separate frequency bands
            for f = [2, 4, 8, 12, 18, 32, 80]
                f_idx = find(abs(concat_states_tfs.freq - f) == ...
                    min(abs(concat_states_tfs.freq - f)), 1, 'first');
                line(xlim, [f_idx f_idx], 'color', 'k', 'linestyle', '--');
            end
        end
    end
    
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    
    % define colormap
    cm = colormap('jet'); % default
    if nargin > 4
        cm = colormap(varargin{1});
        colorbar;
        %end
    end
    
    %cm = colormap('jet'); 
    cm(1,:,:) = [1,1,1];
    colormap(cm);
    
    saveas(gca, results_file);

end

