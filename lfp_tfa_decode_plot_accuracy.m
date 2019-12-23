function lfp_tfa_decode_plot_accuracy( lfp_decode, figtitle, results_folder )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

h = figure('name', figtitle);
set(h, 'position', [100, 100,900, 675]);
noffset = 2;
hold on;
set(gca, 'YLim', [0, 1])

timebin_samples = [0];
timebin_values = [];
event_onset_samples = zeros(1, length(lfp_decode.timebins));
wnd_start_samples = zeros(1, length(lfp_decode.timebins));
wnd_end_samples = zeros(1, length(lfp_decode.timebins));

for ep = 1:length(lfp_decode.timebins)
    event_onset_sample = find(abs(lfp_decode.timebins{ep}) == ...
        min(abs(lfp_decode.timebins{ep})), 1, 'last');
    event_onset_samples(ep) = timebin_samples(end) + event_onset_sample;
    wnd_start_samples(ep) = timebin_samples(end) + 1;
    wnd_end_samples(ep) = timebin_samples(end) + length(lfp_decode.timebins{ep});
    timebin_samples = ...
        timebin_samples(end) + (1:(length(lfp_decode.timebins{ep}) + noffset));
    timebin_values = [timebin_values, lfp_decode.timebins{ep}, nan(1, noffset)];
    shadedErrorBar(timebin_samples(1:length(lfp_decode.timebins{ep})), ...
        nanmean(lfp_decode.test_accuracy{ep}, 2), ...
        nanstd(lfp_decode.test_accuracy{ep}, 0, 2), 'b');

    line([event_onset_samples(ep) event_onset_samples(ep)], ylim, 'color', 'k');
    text(event_onset_samples(ep) + 0.5, 0.5, lfp_decode.epoch_name{ep});

end


title(figtitle);
xlabel('Time (s)');
ylabel('Accuracy');
xticks = [wnd_start_samples; event_onset_samples; wnd_end_samples];
xticks = unique(xticks(:));
set(gca, 'xtick', xticks)
set(gca, 'xticklabels', round(timebin_values(xticks), 1))
set(gca, 'xticklabelrotation', 45);

% save figure
results_file = fullfile(results_folder, figtitle);
try
    export_fig(h, results_file, '-png');
    export_fig(h, results_file, '-pdf');
catch e
    warning('Cannot save figures. Reason: %s\n\n', e.message());
end

end

