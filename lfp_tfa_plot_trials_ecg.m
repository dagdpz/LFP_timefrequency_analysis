function [ h ] = lfp_tfa_plot_trials_ecg( trial_ecg )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    h = figure(1); clf; hold on;

%             h1 = uicontrol('Position', [5 5 200 40], 'String', 'Continue', ...
%                       'Callback', 'uiresume(gcbf)');
%             h2 = uicontrol('Position', [300 5 200 40], 'String', 'Continue', ...
%                       'Callback', 'plottrials = 0; close;');
    % plot raw ECG
    ax1 = subplot(311);
    plot(trial_ecg.time, trial_ecg.ecg_data);
    box on;
    set(gca, 'Xlim', [trial_ecg.time(1), ...
        trial_ecg.time(end)]);
    set(gca, 'Ylim', ax1.YLim);

    % plot spikes
    subplot(312);
    plot(trial_ecg.time, trial_ecg.ECG_spikes);
    box on;
    set(gca, 'Xlim', [trial_ecg.time(1), trial_ecg.time(end)]);

    % plot ECG R2Rt valid
    subplot(313);
    plot(trial_ecg.time, trial_ecg.ECG_b2btime);
    box on;
    set(gca, 'Xlim', [trial_ecg.time(1), trial_ecg.time(end)]);

    % title
    subplot(311);
    title();
    % mark state onsets
    for state = (trial_ecg.states)
        if ~isnan(state.onset_t) || ~isempty( state.onset_t)
            line([state.onset_t state.onset_t], ax1.YLim, ...
                'Color', 'k', 'LineStyle', '--');
            text(double(state.onset_t), ax1.YLim(1) + ((ax1.YLim(2) - ax1.YLim(1))*(find([trial_ecg.states.id] == state.id)/length(trial_ecg.states))), num2str(state.id));
        end
    end

%             uiwait(gcf);

    if t == 1
        disp('Press any key to continue! To abort, press Ctrl+C');
    end
    pause;

end

