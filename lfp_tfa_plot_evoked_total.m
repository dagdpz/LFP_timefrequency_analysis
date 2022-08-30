function lfp_tfa_plot_evoked_total(sessions_avg, lfp_tfa_cfg)
%lfp_tfa_plot_evoked_total plots the average amplitude of each of the four
%conditions (specific for Keren's analysis) - CS/IS, Choice 0/1, hopefully
%with significance also plotted

% USAGE:
%       lfp_tfa_plot_evoked_total(sessions_avg, lfp_tfa_cfg, plottitle, results_file)
% INPUTS:
%       sessions_avg    - a structure containing lfp amplitude data, averaged
%       across sessions. (lfp_avoked.sessions_avg)
%       lfp_tfa_cgf   - struct containing the required settings
%       plottile - title of the figure
%       results_file - file dir to where the figure will be saved

plottitle = 'CS vs IS, Instr and Choice';
results_file = [lfp_tfa_cfg.analyse_lfp_folder filesep plottitle];
h = figure('Name', plottitle);
set(h, 'position', [100, 100,1000, 1000]);
nstates = size(lfp_tfa_cfg.analyse_states, 1);
ncond = length(lfp_tfa_cfg.conditions);
i = 1;
for cn = 1:ncond
    for st = 1:nstates
        subplot(ncond, nstates, i)
        for hs = 1:length(lfp_tfa_cfg.compare.reach_spaces)
            mean_plot = sessions_avg(1).condition(cn).avg_across_sessions(st,hs).mean;
            time_plot = sessions_avg(1).condition(cn).avg_across_sessions(st,hs).time;
            error_plot = sessions_avg(1).condition(cn).avg_across_sessions(st,hs).error;
            p = plot(time_plot, mean_plot);
            hold on
            x = [time_plot, fliplr(time_plot)];
            bw = [error_plot(1,:), fliplr(error_plot(2,:))];
            col = get(p, 'Color');
            fill(x, bw, col, 'FaceAlpha', .1, 'LineStyle', 'none')
            lgnd{i,hs} = sprintf('"%s %s"', sessions_avg(1).condition(cn).label,...
                cell2mat(lfp_tfa_cfg.compare.reach_spaces(hs)));
        end
        %h X and Y are for plotting significant results
        h_X = sessions_avg(1).difference(cn).hs_tuned_evoked(st,1).lfp.stat_test.h;
        h_X = time_plot(h_X==1);
        h_Y = ones(numel(h_X),1)*max(mean_plot);
        scatter(h_X, h_Y,5,'.')
        legend(lgnd{i,1}, 'sterr', lgnd{i,2}, 'sterr', 'Location', 'Best')
        state_name = cell2mat(lfp_tfa_cfg.analyse_states(st,3));
        line([0 0],ylim)
        title(sprintf('%s, %s', state_name, lfp_tfa_cfg.conditions(cn).label))
        ylabel('LFP amplitude')
        xlabel('time')
        hold off
        i = i+1;
    end
    
end


fig_formats = {'png'}; %default
if isfield(lfp_tfa_cfg, 'save_fig_format') && ~isempty(lfp_tfa_cfg.save_fig_format)
    fig_formats = lfp_tfa_cfg.save_fig_format;
end
for fmt = fig_formats
    export_fig(h, results_file, ['-' fmt{:}]);
end
%% this part is very specific for Keren's analysis
plottitle = 'CS-IS, Instructed vs Choice';
results_file = [lfp_tfa_cfg.analyse_lfp_folder filesep plottitle];

g = figure('Name', plottitle);
set(g, 'position', [100, 100,1000, 1000]);
nhs = 2; %ipsi and contra
nchoice = 2; %instr and choice
i = 1;
for hs = 1:nhs
    subplot(nhs, 1, i)
    for ch = 1:nchoice
        mean_plot = sessions_avg(1).condition(ch).avg_across_sessions(hs).mean;
        time_plot = sessions_avg(1).condition(ch).avg_across_sessions(hs).time;
        error_plot = sessions_avg(1).condition(ch).avg_across_sessions(hs).error;
        p = plot(time_plot, mean_plot);
        hold on
        x = [time_plot, fliplr(time_plot)];
        bw = [error_plot(1,:), fliplr(error_plot(2,:))];
        col = get(p, 'Color');
        fill(x, bw, col, 'FaceAlpha', .1, 'LineStyle', 'none')
        lgnd{i,ch} = sprintf('"%s"', lfp_tfa_cfg.conditions(ch).label);
    end
    %h X and Y are for plotting significant results
    h_X = sessions_avg(1).difference(3).hs_tuned_evoked(hs).lfp.stat_test.h;
    h_X = time_plot(h_X==1);
    h_Y = ones(numel(h_X),1)*max(mean_plot);
    scatter(h_X, h_Y,5,'.')
    legend(lgnd{i,1}, 'sterr', lgnd{i,2}, 'sterr', 'Location', 'Best')
    line([0 0],ylim)
    title(sessions_avg(1).condition(1).avg_across_sessions(hs).hs_label)
    ylabel('LFP amplitude')
    xlabel('time')
    hold off
    i = i+1;
end


fig_formats = {'png'}; %default
if isfield(lfp_tfa_cfg, 'save_fig_format') && ~isempty(lfp_tfa_cfg.save_fig_format)
    fig_formats = lfp_tfa_cfg.save_fig_format;
end
for fmt = fig_formats
    export_fig(g, results_file, ['-' fmt{:}]);
end
end