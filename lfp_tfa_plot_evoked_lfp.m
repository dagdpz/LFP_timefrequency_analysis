function lfp_tfa_plot_evoked_lfp( evoked_lfp, lfp_tfa_cfg, plottitle, results_file )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    figure;

    % loop through handspace
    for hs = 1:size(evoked_lfp, 2)
        if ~isempty([evoked_lfp(:,hs).mean]) &&  ~isempty([evoked_lfp(:,hs).std])
            % concatenate states
            concat_states_lfp = struct();
            concat_states_lfp.mean = [];
            concat_states_lfp.std = [];
            concat_states_lfp.time = [];

            state_info = struct();
            for st = 1:size(evoked_lfp, 1)

                % concatenate mean, std and time of evoked LFP for
                % different states
                concat_states_lfp.mean = [concat_states_lfp.mean, ...
                    evoked_lfp(st, hs).mean, nan(1, 150)];
                concat_states_lfp.std = [concat_states_lfp.std, ...
                    evoked_lfp(st, hs).std, nan(1, 150)];
                concat_states_lfp.time = [concat_states_lfp.time, ...
                    evoked_lfp(st, hs).time, nan(1, 150)];
                concat_states_lfp.label = evoked_lfp(st, hs).hs_label;                

                state_info(st).onset_s = find(...
                    evoked_lfp(st, hs).time <= 0, 1, 'last');
                state_info(st).onset_t = 0;
                state_info(st).start_s = 1;
                state_info(st).start_t = evoked_lfp(st, hs).time(1);
                state_info(st).finish_s = length(evoked_lfp(st, hs).time);
                state_info(st).finish_t = evoked_lfp(st, hs).time(end);                    

                if st > 1
                    state_info(st).start_s = length(evoked_lfp(st-1, hs).time) + ...
                        state_info(st).start_s + (st-1)*150;
                    state_info(st).finish_s = length(evoked_lfp(st-1, hs).time) + ...
                        state_info(st).finish_s + (st-1)*150;
                    state_info(st).onset_s = length(evoked_lfp(st-1, hs).time) + ...
                        state_info(st).onset_s + (st-1)*150;
                end

            end
            %lfp_time = 1:1:length(concat_states_lfp.time);
            state_onsets = find(concat_states_lfp.time(1:end-1) .* ...
                concat_states_lfp.time(2:end) <= 0);
            state_samples = sort([state_info.start_s, state_info.onset_s, ...
                state_info.finish_s]);

            % now plot
            subplot(2,2,hs)
            hold on;
            plot(concat_states_lfp.mean, 'b');
            plot(concat_states_lfp.mean + concat_states_lfp.std, 'r--');
            plot(concat_states_lfp.mean - concat_states_lfp.std, 'r--');
            % mark state onsets
            set(gca,'xtick',state_samples)
%                 xticklabels = [];
            for so = state_onsets
                line([so so], ylim); 
                state = evoked_lfp(state_onsets == so, hs).state;
                state_name = lfp_tfa_cfg.all_states(...
                    [lfp_tfa_cfg.all_states.state_ID] == state).state_name;
                ypos = ylim;
                ypos = ypos(1) + (ypos(2) - ypos(1))*0.2;
                text(so, ypos, state_name);
            end
            set(gca,'xticklabels', round(concat_states_lfp.time(state_samples), 2), 'fontsize', 8)
            xlabel('Time(s)');
            ylabel('LFP amplitude');
            subplottitle = [concat_states_lfp.label{1}];
            if isfield(evoked_lfp(1, hs), 'trials')
                subplottitle = [subplottitle ' (ntrials = ' ...
                    num2str(length(evoked_lfp(1, hs).trials)) ')'];
            elseif isfield(evoked_lfp(1, hs), 'nsites')
                subplottitle = [subplottitle ' (nsites = ' num2str(evoked_lfp(1, hs).nsites) ')'];
            end
            title(subplottitle);
            %line([0 0], ylim, 'color', 'k');
        end
    end
    ann = annotation('textbox', [0 0.9 1 0.1], 'String', strrep(plottitle, '_', '\_')...
        , 'EdgeColor', 'none', 'HorizontalAlignment', 'center');
    saveas(gca, results_file);

end

